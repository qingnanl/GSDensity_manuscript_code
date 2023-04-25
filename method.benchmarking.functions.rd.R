#library(GSVA)
# library(singleCellTK)
# source("sergio_analysis_functions.R")
library(BiocParallel)
library(GSVA)
library(Seurat)
# library(SeuratData)
library(dplyr)
library(ggrepel)
library(ggplot2)
# to be bench marked
library(gsdensity)
library(CelliD)
# library(escape)
library(VAM)
library(DescTools)
library(GSEABase)
library(AUCell)
library(nnet)
library(reshape2)
library(RColorBrewer)

igraph::igraph.options(sparsematrices = FALSE)
# get a function for calculate AUC
# get a function for calculate AUC
calc.auc <- function(input, 
                     pos.id, 
                     method.id,
                     sample.points = c(seq(from = 0, to = 0.2, by = 0.025), seq(from = 0.3, to = 1, by = 0.1))){
  total.pos <- sum(input[, 2] == pos.id)
  input.order <- input[order(input[, 1], decreasing = T), ]
  output.df <- as.data.frame(sample.points)
  output.df$positive.capture <- 0
  output.df$method <- method.id
  for (i in 2:length(sample.points)){
    pos.count <- sum(input.order[1:floor(nrow(input.order)*sample.points[i]), 2] == pos.id)
    output.df[i, "positive.capture"] <- pos.count/total.pos
  }
  return(output.df)
}

# get cell-level activity scores
# gsd
get.activity.list <- function(glist, pos.id, object, el, cells, cells_rankings){
  lst <- list()
  print(glist)
  gsd <- run.rwr(el = el, gene_set = glist, cells = cells)
  gsd.res <- as.data.frame(gsd)
  gsd.res$cell.type <- object$cell.type
  lst[['gsdensity']] <- gsd.res
  # gsea
  #gsea <- enrichIt(obj = object, 
  #                 gene.sets = list(glist = glist), groups = ncol(object))
  mt <- match(CharacterList(list(glist)), rownames(object))
  mapdgenesets <- as.list(mt[!is.na(mt)])
  
  a <- GSVA:::.gsva(expr = object@assays$RNA@data, gset.idx.list = mapdgenesets, abs.ranking = F,tau = 0.25, rnaseq = T,
                    method = "ssgsea", kcdf = "Poisson")
  gsea.res <- as.data.frame(t(a))
  rownames(gsea.res) <- colnames(object)
  colnames(gsea.res) <- 'gsea'
  gsea.res$cell.type <- object$cell.type
  lst[['gsea']] <- gsea.res
  
  b <- GSVA:::.gsva(expr = object@assays$RNA@data, gset.idx.list = mapdgenesets, method = "gsva", rnaseq = F, 
                    kernel = F, parallel.sz = 40, BPPARAM = SnowParam(), verbose = T)
  gsva.res <- as.data.frame(t(b))
  rownames(gsva.res) <- colnames(object)
  colnames(gsva.res) <- 'gsva'
  gsva.res$cell.type <- object$cell.type
  lst[['gsva']] <- gsva.res
  # cellid
  
  HGT_glist <- RunCellHGT(object, pathways = list(glist = glist), dims = 1:30, minSize = 5)
  cellid.res <- as.data.frame(HGT_glist)
  colnames(cellid.res) <- 'cellid'
  cellid.res$cell.type <- object$cell.type
  lst[['cellid']] <- cellid.res
  #
  
  geneSets <- list(geneSet1 = glist)
  cells_AUC <- AUCell_calcAUC(geneSets, rankings = cells_rankings)
  cells_AUC_df <- as.data.frame(t(cells_AUC@assays@data@listData$AUC))
  aucell.res <- as.data.frame(cells_AUC_df)
  colnames(aucell.res) <- 'aucell'
  aucell.res$cell.type <- object$cell.type
  lst[['aucell']] <- aucell.res
  
  coll <- list()
  coll[[1]] <- intersect(rownames(object), glist)
  names(coll)[1] <- "VAM"
  gene.set.collection = createGeneSetCollection(gene.ids=rownames(object),
                                                gene.set.collection=coll)
  obj <- vamForSeurat(seurat.data=object,
                      gene.set.collection=gene.set.collection,
                      center=F, gamma=T, sample.cov=F, return.dist=F)
  vam.res <- as.data.frame(t(obj@assays$VAMcdf@data))
  vam.res$cell.type <- object$cell.type
  lst[['vam']] <- vam.res
  
  lst1 <- list()
  for (method in names(lst)){
    lst1[[method]] <- calc.auc(input = lst[[method]], method.id = method, pos.id = pos.id)
  }
  return(list(score.list = lst, auc.list = lst1))
}
benchmark.AUC <- function(object, clist, methods, markers, el, cells, cells_rankings){
  AUC.collection <- list()
  for (i in 1:length(clist)){
    # for cells
    glist <- markers[[clist[i]]]
    # print(glist)
    ga <- get.activity.list(glist = glist, pos.id = clist[i], object = object, 
                            el = el, cells = cells, cells_rankings = cells_rankings)
    AUC.list <- c()
    for (method in names(ga$auc.list)){
      AUC.list <- c(AUC.list, AUC(x = ga$auc.list[[method]][, 'sample.points'],
                                  y = ga$auc.list[[method]][, 'positive.capture']))
    }
    AUC.df <- data.frame(AUC.score = AUC.list)
    AUC.df$methods <- methods# c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam")
    AUC.df$cell.type <- clist[i]
    AUC.collection[[clist[i]]] <- AUC.df
  }
  AUC.col.plot <- do.call(rbind, AUC.collection)
  return(AUC.col.plot)
}
get.gsd.activity.list <- function(glist, pos.id, object, el, cells){
  lst <- list()
#  print(glist)
  gsd <- run.rwr.list(el = el, gene_set = glist, cells = cells)
  gsd.res <- as.data.frame(gsd)
  gsd.res$cell.type <- object$cell.type
  n.cell.types <- length(glist)
  gsd.auc <- lapply(1:n.cell.types,
                   FUN = function(x){calc.auc(input = gsd.res[, c(pos.id[x], "cell.type")], 
                                              method.id = "gsdensity", 
                                              pos.id = pos.id[x])})
  names(gsd.auc) <- names(glist)
  #lst[['gsdensity']] <- gsd.res
  #lst1 <- list()
  #for (method in names(lst)){
  #  lst1[[method]] <- calc.auc(input = lst[[method]], method.id = method, pos.id = pos.id)
  #}
  return(gsd.auc)
  # return(list(score.list = lst, auc.list = lst1))
}
benchmark.AUC.2 <- function(object, clist, methods = "gsdensity", markers, el, cells){
 # AUC.collection <- list()
  glist <- markers
  ga <- get.gsd.activity.list(glist = glist, pos.id = clist, object = object, 
                              el = el, cells = cells)
  AUC.list <- c()
  for (cell in names(ga)){
    AUC.list <- c(AUC.list, AUC(x = ga[[cell]][, 'sample.points'],
                                y = ga[[cell]][, 'positive.capture']))
  }
  AUC.df <- data.frame(AUC.score = AUC.list)
  # AUC.df$methods <- methods# c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam")
  AUC.df$cell.type <- names(ga)
  return(AUC.df)
#  AUC.col.plot <- do.call(rbind, AUC.collection)
#  return(AUC.col.plot)
}


benchmark.ACC <- function(object, clist, methods, markers, el, cells, cells_rankings){
  score.collection <- list()
  for (i in 1:length(clist)){
    # for cells
    glist <- markers[[clist[i]]]
    # print(glist)
    ga <- get.activity.list(glist = glist, pos.id = clist[i], object = object, 
                            el = el, cells = cells, cells_rankings = cells_rankings)
    score.df <- do.call(cbind, ga$score.list)[, c(1, 3, 5, 7, 9, 11, 12)]
    colnames(score.df) <- c(methods, "cell.type")
    score.collection[[clist[i]]] <- score.df
  }
  method.vec <- list()
  cell.vec <- list()
  for (method in methods){
    for (cell in clist){
      cell.vec[[cell]] <-  score.collection[[cell]][, method]
    }
    cell.df <- do.call(cbind, cell.vec)
    method.vec[[method]] <- as.data.frame(cell.df)
  }
  
  for (method in methods){
    method.vec[[method]]$max <- colnames(method.vec[[method]])[apply(method.vec[[method]], 1, which.is.max)]
    method.vec[[method]]$cell.type <- score.collection[[clist[1]]][, "cell.type"]
  }
  
  # method.output <- do.call(cbind, method.vec)
  
  ACC.list <- c()
  for (method in methods){
    ACC.list <- c(ACC.list, ACC.calc(ACC.df = method.vec[[method]]))
  }
  ACC.df <- data.frame(ACC.score = ACC.list)
  ACC.df$methods <- methods
  return(list(ACC.df = ACC.df, method.raw = method.vec))
}

plot.AUC <- function(AUC.df){
  p <- ggplot(data = AUC.col.plot, aes(x = methods, y = AUC.score, fill = methods))+
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_brewer(palette="BuPu") +
    geom_point(data = AUC.col.plot, aes(x = methods, y = AUC.score, color = cell.type), 
               shape = 2) + 
    theme_classic() + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15))
  return(p)
}



ACC.calc <- function(ACC.df){
  ACC.rate <- sum(ACC.df$max == ACC.df$cell.type)/nrow(ACC.df)
  return(ACC.rate)
}

benchmark.real <- function(fname, k.range = 5, fixed.ds = 10){
  sr.fname <- paste0(fname, ".rds")
  sr.markers.fname <- paste0(fname, ".markers.rds")
  # obtain seurat object
  sr <- preprocess_sr(sr.fname)
  # print(sr)
  # pre-compute
  # gsdensity
  ce <- compute.mca(object = sr)
  # print(head(ce))
  cells <- colnames(sr)
  el <- compute.nn.edges(coembed = ce, nn.use = 300)
  # print(head(el))
  # cellid
  sr <- RunMCA(sr)
  # aucell
  cells_rankings <- AUCell_buildRankings(sr@assays$RNA@data, plotStats = F, splitByBlocks = TRUE)
  
  
  # obtain gene list sum
  #sum <- preprocess_regulation(design = design, version.n = version.n, n.tight.strong.targets.r = n.tight.strong.targets.r, 
  #                             n.tight.medium.targets.r = n.tight.medium.targets.r, n.tight.weak.targets.r = n.tight.weak.targets.r, 
  #                             n.multi.targets.r = n.multi.targets.r, n.rest.r = n.rest.r, n.cell.types = n.cell.types, n.genes = n.genes)
  
  # print(str(sum))
  #tar.type <- c("strong", "medium", "weak")
  marker.list <- readRDS(sr.markers.fname)
  clist <- names(marker.list)
  print(clist)
  n.cell.types <- length(clist)
  # AUC.tar <- list()
  AUC.batch <- list()
  AUC.res <- benchmark.AUC(object = sr, clist = clist, 
                           methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam"), 
                           markers = marker.list, el = el, cells = cells, cells_rankings = cells_rankings)
  AUC.res$design <- "full.marker.set"
  # AUC.res$drop.out <- drop.out
  AUC.res$batch <- 1
  # AUC.res$target.type <- tar.type.name
  AUC.batch[["full.marker.set"]] <- AUC.res
  for (k in 1:k.range){
    ds.markers <- list()
    for (q in 1:n.cell.types){
      ds.markers[[q]] <- ds.glist(glist = marker.list[[q]], all.genes = rownames(sr), mode = "ratio", ratio = 0.7, s = k)
    }
    names(ds.markers) <- clist
    
    AUC.res <- benchmark.AUC(object = sr, clist = clist, 
                             methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam"), 
                             markers = ds.markers, el = el, cells = cells, cells_rankings = cells_rankings)
    AUC.res$design <- "70_percent.ds"
    # AUC.res$drop.out <- drop.out
    AUC.res$batch <- k
    # AUC.res$target.type <- tar.type.name
    AUC.batch[[paste0("70_percent.ds.", k)]] <- AUC.res
  }
  for (k in 1:k.range){
    ds.markers <- list()
    for (q in 1:n.cell.types){
      ds.markers[[q]] <- ds.glist(glist = marker.list[[q]], all.genes = rownames(sr), mode = "ratio", ratio = 0.3, s = k)
    }
    names(ds.markers) <- clist
    
    AUC.res <- benchmark.AUC(object = sr, clist = clist, 
                             methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam"), 
                             markers = ds.markers, el = el, cells = cells, cells_rankings = cells_rankings)
    AUC.res$design <- "30_percent.ds"
    # AUC.res$drop.out <- drop.out
    AUC.res$batch <- k
    # AUC.res$target.type <- tar.type.name
    AUC.batch[[paste0("30_percent.ds.", k)]] <- AUC.res
  }
  for (k in 1:k.range){
    ds.markers <- list()
    for (q in 1:n.cell.types){
      ds.markers[[q]] <- ds.glist.mix(glist = marker.list[[q]], all.genes = rownames(sr), ratio = 0.3, s = k, mix.ratio = 1)
    }
    names(ds.markers) <- clist
    
    AUC.res <- benchmark.AUC(object = sr, clist = clist, 
                             methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam"), 
                             markers = ds.markers, el = el, cells = cells, cells_rankings = cells_rankings)
    AUC.res$design <- "30_percent.ds.mix1"
    # AUC.res$drop.out <- drop.out
    AUC.res$batch <- k
    # AUC.res$target.type <- tar.type.name
    AUC.batch[[paste0("30_percent.ds.mix1.", k)]] <- AUC.res
  }
  for (k in 1:k.range){
    ds.markers <- list()
    for (q in 1:n.cell.types){
      ds.markers[[q]] <- ds.glist.mix(glist = marker.list[[q]], all.genes = rownames(sr), ratio = 0.3, s = k, mix.ratio = 3)
    }
    names(ds.markers) <- clist
    
    AUC.res <- benchmark.AUC(object = sr, clist = clist, 
                             methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam"), 
                             markers = ds.markers, el = el, cells = cells, cells_rankings = cells_rankings)
    AUC.res$design <- "30_percent.ds.mix3"
    # AUC.res$drop.out <- drop.out
    AUC.res$batch <- k
    # AUC.res$target.type <- tar.type.name
    AUC.batch[[paste0("30_percent.ds.mix3.", k)]] <- AUC.res
  }
  for (k in 1:k.range){
    ds.markers <- list()
    for (q in 1:n.cell.types){
      ds.markers[[q]] <- ds.glist.mix(glist = marker.list[[q]], all.genes = rownames(sr), ratio = 0.3, s = k, mix.ratio = 5)
    }
    names(ds.markers) <- clist
    
    AUC.res <- benchmark.AUC(object = sr, clist = clist, 
                             methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam"), 
                             markers = ds.markers, el = el, cells = cells, cells_rankings = cells_rankings)
    AUC.res$design <- "30_percent.ds.mix5"
    # AUC.res$drop.out <- drop.out
    AUC.res$batch <- k
    # AUC.res$target.type <- tar.type.name
    AUC.batch[[paste0("30_percent.ds.mix5.", k)]] <- AUC.res
  }
  AUC.tar <- as.data.frame(do.call(rbind, AUC.batch))
  
  ACC.batch <- list()
  # ACC.tar <- list()
  method.batch <- list()
  for (k in 1:k.range){
    ds.markers <- list()
    for (q in 1:n.cell.types){
      ds.markers[[q]] <- ds.glist(glist = marker.list[[q]], all.genes = rownames(sr), fixed = fixed.ds, s = k)
    }
    names(ds.markers) <- clist
    # print(ds.markers)
    ACC.full.res <- benchmark.ACC(object = sr, clist = clist, 
                                  methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam"), 
                                  markers = ds.markers, el = el, cells = cells, cells_rankings = cells_rankings)
    ACC.res <- ACC.full.res[["ACC.df"]]
    ACC.res$design <- fname
    #ACC.res$drop.out <- drop.out
    ACC.res$batch <- k
    # ACC.res$target.type <- tar.type.name
    ACC.batch[[k]] <- ACC.res
    method.batch[[k]] <- ACC.full.res[["method.raw"]]
  }
  ACC.tar <- as.data.frame(do.call(rbind, ACC.batch))
  # AUC.tar[[p]] <- as.data.frame(do.call(rbind, AUC.batch))
  
  data.sum <- list(object = sr, ACC = ACC.tar, AUC = AUC.tar, method.raw = method.batch)
  saveRDS(data.sum$object, paste0(fname, ".sr.rds"))
  saveRDS(data.sum[2:4], paste0(fname, ".info.rds"))
  return(data.sum)
}

plot.ACC <- function(ACC.df){
  p <- ggplot(data = ACC.df, aes(x = methods, y = ACC.score, fill = methods))+
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_brewer(palette="BuPu") +
    geom_jitter(data = ACC.df, aes(x = methods, y = ACC.score, alpha = 0.5), 
                shape = 2, size = 0.5) + 
    theme_classic() + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15))
  return(p)
}



benchmark.sergio <- function(design, version.n, drop.out, n.cell.types, n.cells.per.type, 
                             n.tight.strong.targets.r = 0.1, 
                             n.tight.medium.targets.r = 0.1, n.tight.weak.targets.r = 0.1, 
                             n.multi.targets.r = 0.3, n.rest.r = 0.395, n.genes, k.range = 5, fixed.ds = 10){
  # obtain seurat object
  sr <- preprocess_sr(design = design, version.n = version.n, drop.out = drop.out, 
                      n.cell.types = n.cell.types, n.cells.per.type = n.cells.per.type)
  # print(sr)
  # pre-compute
  # gsdensity
  ce <- compute.mca(object = sr)
  # print(head(ce))
  cells <- colnames(sr)
  el <- compute.nn.edges(coembed = ce, nn.use = 300)
  # print(head(el))
  # cellid
  sr <- RunMCA(sr)
  # aucell
  cells_rankings <- AUCell_buildRankings(sr@assays$RNA@data, plotStats = F, splitByBlocks = TRUE)
  
  clist <- paste0("cell", 1:n.cell.types)
  # print(clist)
  # obtain gene list sum
  sum <- preprocess_regulation(design = design, version.n = version.n, n.tight.strong.targets.r = n.tight.strong.targets.r, 
                               n.tight.medium.targets.r = n.tight.medium.targets.r, n.tight.weak.targets.r = n.tight.weak.targets.r, 
                               n.multi.targets.r = n.multi.targets.r, n.rest.r = n.rest.r, n.cell.types = n.cell.types, n.genes = n.genes)
  
  # print(str(sum))
  tar.type <- c("strong", "medium", "weak")
  
  ACC.tar <- list()
  AUC.tar <- list()
  for (p in 1:3){
    # markers <- sum$target.list[c(p+0, p+3, p+6)]
    markers <- sum$target.list[seq(0, 3*(n.cell.types - 1), by = 3) + p]
    tar.type.name <- tar.type[p]
    names(markers) <- clist
    # print(markers)
    ACC.batch <- list()
    AUC.batch <- list()
    for (k in 1:k.range){
      ds.markers <- list()
     # ds.markers[[1]] <- ds.glist(glist = markers[[1]], all.genes = rownames(sr), fixed = fixed.ds, s = k)
     # ds.markers[[2]] <- ds.glist(glist = markers[[2]], all.genes = rownames(sr), fixed = fixed.ds, s = k)
     # ds.markers[[3]] <- ds.glist(glist = markers[[3]], all.genes = rownames(sr), fixed = fixed.ds, s = k)
      for (q in 1:n.cell.types){
        ds.markers[[q]] <- ds.glist(glist = markers[[q]], all.genes = rownames(sr), fixed = fixed.ds, s = k)
      }
      names(ds.markers) <- clist
      # print(ds.markers)
      ACC.res <- benchmark.ACC(object = sr, clist = clist, 
                               methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam"), 
                               markers = ds.markers, el = el, cells = cells, cells_rankings = cells_rankings)
      ACC.res$design.n <- paste0(design, ".", version.n)
      ACC.res$drop.out <- drop.out
      ACC.res$batch <- k
      ACC.res$target.type <- tar.type.name
      ACC.batch[[k]] <- ACC.res
      
      AUC.res <- benchmark.AUC(object = sr, clist = clist, 
                               methods = c("gsdensity", "gsea", "gsva","cellid", "aucell", "vam"), 
                               markers = ds.markers, el = el, cells = cells, cells_rankings = cells_rankings)
      AUC.res$design.n <- paste0(design, ".", version.n)
      AUC.res$drop.out <- drop.out
      AUC.res$batch <- k
      AUC.res$target.type <- tar.type.name
      AUC.batch[[k]] <- AUC.res
    }
    ACC.tar[[p]] <- as.data.frame(do.call(rbind, ACC.batch))
    AUC.tar[[p]] <- as.data.frame(do.call(rbind, AUC.batch))
  }
  data.sum <- list(object = sr, ACC = ACC.tar, AUC = AUC.tar, regulatory_summary = sum)
  saveRDS(data.sum$object, paste0(design, ".", version.n, ".", drop.out, ".sr.rds"))
  saveRDS(data.sum[2:4], paste0(design, ".", version.n, ".", drop.out, ".info.rds"))
  return(data.sum)
}

preprocess_sr <- function(sr.fname){
  sr <- readRDS(sr.fname)
  sr <- sr %>%
    NormalizeData() %>% 
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors() %>%
    RunUMAP(dims = 1:10)
  return(sr)
}

ds.glist <- function(glist, all.genes, ratio, fixed, mode = "fixed", bottom = 5, s){
  set.seed(s)
  glist <- intersect(glist, all.genes)
  if (mode == "fixed"){
    sample.size <- fixed
  }
  if (mode == "ratio"){
    sample.size <- max(floor(length(glist) * ratio), bottom)
  }
  idx <- sample(x = 1:length(glist), size = sample.size, replace = F)
  out.glist <- glist[idx]
  return(out.glist)
}

ds.glist.mix <- function(glist, all.genes, ratio, bottom = 5, mix.ratio = 1, s){
  set.seed(s)
  glist <- intersect(glist, all.genes)
  sample.size <- max(floor(length(glist) * ratio), bottom)
  idx <- sample(x = 1:length(glist), size = sample.size, replace = F)
  out.glist <- glist[idx]
  bg.genes <- setdiff(all.genes, glist)
  sample.bg.size <- floor(sample.size * mix.ratio)
  bg.idx <- sample(x = 1:length(bg.genes), size = sample.bg.size, replace = F)
  bg.sample <- bg.genes[bg.idx]
  out.glist <- c(out.glist, bg.sample)
  return(out.glist)
}


get.prediction <- function(info){
  output.list <- list(gsdensity = info[["method.raw"]][[1]][["gsdensity"]]$max, 
                      gsea = info[["method.raw"]][[1]][["gsea"]]$max,
                      gsva = info[["method.raw"]][[1]][["gsva"]]$max,
                      cellid = info[["method.raw"]][[1]][["cellid"]]$max,
                      aucell = info[["method.raw"]][[1]][["aucell"]]$max,
                      vam = info[["method.raw"]][[1]][["vam"]]$max)
  return(as.data.frame(do.call(cbind, output.list)))
}

plot.pred.umap <- function(object, fname, info){
  cell.type.prediction <- get.prediction(info)
  object@meta.data <- cbind(object@meta.data, cell.type.prediction)
  pdf(paste0(fname, ".pred.umap.pdf"), width = 12, height = 10)
  p <- DimPlot(object, group.by = c("gsdensity", "gsea", "aucell", "cellid", "gsva", "vam","cell.type"))
  print(p)
  dev.off()
}

print.auc <- function(info, fname){
  auc <- info$AUC
  designs <- unique(auc$design)
  for (des in designs){
    pdf(paste0(fname, ".", des, ".AUC.pdf"), width = 8, height = 6)
    p <- plot.AUC(auc[auc$design == des, ])
    print(p)
    dev.off()
  }
}

print.acc <- function(info, fname){
  acc <- info$ACC
  pdf(paste0(fname, ".ACC.pdf"), width = 8, height = 6)
  p <- plot.ACC(acc)
  print(p)
  dev.off()
}

print.auc.2 <- function(info, fname){
  auc <- info$AUC
  pdf(paste0(fname, ".sum.AUC.pdf"), width = 8, height = 6)
  p <- plot.AUC.2(auc)
  print(p)
  dev.off()
}

plot.AUC.2 <- function(AUC.df){
  AUC.df$design <- factor(AUC.df$design, levels=c("full.marker.set", "70_percent.ds", "30_percent.ds", 
                                                  "30_percent.ds.mix1", "30_percent.ds.mix3", "30_percent.ds.mix5"))
  p <- ggplot(data = AUC.df, aes(x = methods, y = AUC.score, fill = design))+
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_brewer(palette="BuPu") +
    geom_jitter(data = AUC.df, aes(x = methods, y = AUC.score, alpha = 0.5), 
                shape = 2, size = 0.5) + 
    theme_classic() + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15))
  return(p)
}

output.auc <- function(info){
  auc <- info$AUC
  designs <- unique(auc$design)
  auc.ls <- list()
  for (des in designs){
    dat <- auc[auc$design == des, ] %>%
      group_by(methods) %>%
      summarise_at(vars(AUC.score), list(name = mean)) %>% 
      as.data.frame()
    dat <- as.data.frame(dat$name)
    colnames(dat) <- des
    auc.ls[[des]] <- dat
  }
  out <- do.call(cbind, auc.ls)
  rownames(out) <- c("aucell", "cellid", "gsdensity", "gsea", "gsva", "vam")
  return(out)
}

output.acc <- function(info, design){
  acc <- info$ACC
  dat <- acc %>%
    group_by(methods) %>%
    summarise_at(vars(ACC.score), list(name = median)) %>% 
    as.data.frame()
  dat <- as.data.frame(dat$name)
  colnames(dat) <- "accuracy"
  rownames(dat) <- c("aucell", "cellid", "gsdensity", "gsea", "gsva", "vam")
  dat$design <- design
  dat$method <- rownames(dat)
  return(dat)
}

make.bubble.auc <- function(info, fname){
  plt.df <- melt(output.auc(info))
  plt.df$method <- rep(c("aucell", "cellid", "gsdensity", "gsea", "gsva", "vam"), 6)
  
  p <- ggplot(plt.df, aes(x = variable, y = method, label = value)) + 
    geom_point(aes(size = value, fill = value), shape = 21) +
    theme_classic() + 
    scale_x_discrete(position = "top") +
    labs( x= "", y = "", size = "AUC score", fill = "AUC score") + 
    scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), limits=c(0.5, 1)) + 
    scale_size_continuous(range = c(0, 10)) + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 0, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2))
  pdf(paste0(fname, ".bubble.AUC.pdf"), width = 7, height = 5)
  print(p)
  dev.off()
}

acc.nbrs.real <- function(fname, k.range = 5, fixed.ds = 10){
  sr.fname <- paste0(fname, ".sr.rds")
  sr <- readRDS(sr.fname)
  sr.markers.fname <- paste0(fname, ".markers.rds")
  # obtain seurat object
  # sr <- preprocess_sr(sr.fname)
  # print(sr)
  # pre-compute
  # gsdensity
  ce <- compute.mca(object = sr)
  # print(head(ce))
  cells <- colnames(sr)
  marker.list <- readRDS(sr.markers.fname)
  clist <- names(marker.list)
  print(clist)
  n.cell.types <- length(clist)
  #  ACC.batch <- list()
  # ACC.tar <- list()
  ACC.nbrs.batch <- list()
  AUC.nbrs.batch <- list()
  for (n.nbrs in c(20, 50, 100, 200, 300, 400, 500, 600)){ #
    n.nbrs.label <- paste0("nn.", n.nbrs)
    ACC.nbrs.batch[[n.nbrs.label]] <- c() 
    el <- compute.nn.edges(coembed = ce, nn.use = n.nbrs)
    for (k in 1:k.range){
      ds.markers <- list()
      for (q in 1:n.cell.types){
        ds.markers[[q]] <- ds.glist(glist = marker.list[[q]], all.genes = rownames(sr), fixed = fixed.ds, s = k)
      }
      names(ds.markers) <- clist
      gsd <- run.rwr.list(el = el, gene_set = ds.markers, cells = cells)
      gsd.res <- as.data.frame(gsd)
      gsd.res$max <- colnames(gsd.res)[apply(gsd.res, 1, which.is.max)]
      gsd.res$cell.type <- sr$cell.type
      ACC.nbrs.batch[[n.nbrs.label]] <- c(ACC.nbrs.batch[[n.nbrs.label]], ACC.calc(ACC.df = gsd.res))
    }
    # print(ds.markers)
    ACC.nbrs.batch[[n.nbrs.label]] <- as.data.frame(ACC.nbrs.batch[[n.nbrs.label]])
    
    #   ACC.res <- ACC.full.res[["ACC.df"]]
    ACC.nbrs.batch[[n.nbrs.label]]$design <- fname
    #ACC.res$drop.out <- drop.out
    #    nbrs.batch[[n.nbrs.label]]$batch <- k
    ACC.nbrs.batch[[n.nbrs.label]]$n.nbrs <- n.nbrs.label
    # ACC.res$target.type <- tar.type.name
    #   ACC.batch[[k]] <- ACC.res
    #   method.batch[[k]] <- ACC.full.res[["method.raw"]]
    # AUC
    AUC.batch <- list()
    AUC.res <- benchmark.AUC.2(object = sr, clist = clist, 
                               methods = c("gsdensity"), 
                               markers = marker.list, el = el, cells = cells)
    AUC.res$design <- "full.marker.set"
    # AUC.res$drop.out <- drop.out
    AUC.res$batch <- 1
    AUC.res$n.nbrs <- n.nbrs.label
    # AUC.res$target.type <- tar.type.name
    AUC.nbrs.batch[[paste0("full.marker.set.", n.nbrs.label)]] <- AUC.res
    for (k in 1:k.range){
      ds.markers <- list()
      for (q in 1:n.cell.types){
        ds.markers[[q]] <- ds.glist.mix(glist = marker.list[[q]], all.genes = rownames(sr), ratio = 0.3, s = k, mix.ratio = 1)
      }
      names(ds.markers) <- clist
      AUC.res <- benchmark.AUC.2(object = sr, clist = clist, 
                                 methods = c("gsdensity"), 
                                 markers = ds.markers, el = el, cells = cells)
      AUC.res$design <- "30_percent.ds.mix1"
      # AUC.res$drop.out <- drop.out
      AUC.res$batch <- k
      AUC.res$n.nbrs <- n.nbrs.label
      # AUC.res$target.type <- tar.type.name
      AUC.nbrs.batch[[paste0("30_percent.ds.mix1.", k, ".", n.nbrs.label)]] <- AUC.res
    }
    for (k in 1:k.range){
      ds.markers <- list()
      for (q in 1:n.cell.types){
        ds.markers[[q]] <- ds.glist.mix(glist = marker.list[[q]], all.genes = rownames(sr), ratio = 0.3, s = k, mix.ratio = 5)
      }
      names(ds.markers) <- clist
      AUC.res <- benchmark.AUC.2(object = sr, clist = clist, 
                                 methods = c("gsdensity"), 
                                 markers = ds.markers, el = el, cells = cells)
      AUC.res$design <- "30_percent.ds.mix5"
      # AUC.res$drop.out <- drop.out
      AUC.res$batch <- k
      AUC.res$n.nbrs <- n.nbrs.label
      # AUC.res$target.type <- tar.type.name
      AUC.nbrs.batch[[paste0("30_percent.ds.mix5.", k, ".", n.nbrs.label)]] <- AUC.res
    }
  }
  ACC.tar <- as.data.frame(do.call(rbind, ACC.nbrs.batch))
  AUC.tar <- as.data.frame(do.call(rbind, AUC.nbrs.batch))
  return(list(ACC.tar, AUC.tar))
  #   el <- compute.nn.edges(coembed = ce, nn.use = 300)
}


plot.AUC.3 <- function(AUC.df, fname){
  AUC.df$n.nbrs <- factor(AUC.df$n.nbrs, levels=c("nn.20", "nn.50", "nn.100", "nn.200", 
                                                  "nn.300", "nn.400", "nn.500", "nn.600"))
  AUC.df.1 <- AUC.df[AUC.df$design == "full.marker.set", ]
  p <- ggplot(data = AUC.df.1, aes(x = n.nbrs, y = AUC.score))+
    geom_boxplot(outlier.shape = NA, fill = "lightgrey") + 
    #   scale_fill_brewer(palette="BuPu") +
    coord_cartesian(ylim=c(0,1)) +
    geom_jitter(data = AUC.df.1, aes(x = n.nbrs, y = AUC.score, alpha = 0.5), 
                shape = 2, size = 0.5) + 
    theme_classic() + 
    NoLegend() + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15))
  pdf(paste0(fname, ".full.marker.set.scan.nn.auc.pdf"), width = 5, height = 3.5)
  print(p)
  dev.off()
  AUC.df.2 <- AUC.df[AUC.df$design == "30_percent.ds.mix1", ]
  p <- ggplot(data = AUC.df.2, aes(x = n.nbrs, y = AUC.score))+
    geom_boxplot(outlier.shape = NA, fill = "lightgrey") + 
    #   scale_fill_brewer(palette="BuPu") +
    coord_cartesian(ylim=c(0,1)) +
    geom_jitter(data = AUC.df.2, aes(x = n.nbrs, y = AUC.score, alpha = 0.5), 
                shape = 2, size = 0.5) + 
    theme_classic() + 
    NoLegend() + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15))
  pdf(paste0(fname, ".30_percent.ds.mix1.scan.nn.auc.pdf"), width = 5, height = 3.5)
  print(p)
  dev.off()
  AUC.df.3 <- AUC.df[AUC.df$design == "30_percent.ds.mix5", ]
  p <- ggplot(data = AUC.df.3, aes(x = n.nbrs, y = AUC.score))+
    geom_boxplot(outlier.shape = NA, fill = "lightgrey") + 
    #   scale_fill_brewer(palette="BuPu") +
    coord_cartesian(ylim=c(0,1)) +
    geom_jitter(data = AUC.df.3, aes(x = n.nbrs, y = AUC.score, alpha = 0.5), 
                shape = 2, size = 0.5) + 
    theme_classic() + 
    NoLegend() + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15))
  pdf(paste0(fname, ".30_percent.ds.mix5.scan.nn.auc.pdf"), width = 5, height = 3.5)
  print(p)
  dev.off()
}
plot.ACC.3 <- function(ACC.df, fname){
  ACC.df$n.nbrs <- factor(ACC.df$n.nbrs, levels=c("nn.20", "nn.50", "nn.100", "nn.200", 
                                                  "nn.300", "nn.400", "nn.500", "nn.600"))
  colnames(ACC.df)[1] <- "ACC.score"
  p <- ggplot(data = ACC.df, aes(x = n.nbrs, y = ACC.score))+
    geom_boxplot(outlier.shape = NA, fill = "lightgrey") + 
    #    scale_fill_brewer(palette="BuPu") +
    coord_cartesian(ylim=c(0,1)) +
    geom_jitter(data = ACC.df, aes(x = n.nbrs, y = ACC.score, alpha = 0.5), 
                shape = 2, size = 1.5) + 
    theme_classic() + 
    NoLegend() + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15))
  pdf(paste0(fname, ".scan.nn.acc.pdf"), width = 5, height = 3.5)
  print(p)
  dev.off()
  #  return(p)
}