library(GSVA)
# library(singleCellTK)
# source("sergio_analysis_functions.R")
library(BiocParallel)
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
  
  b <- GSVA:::.gsva(expr = object@assays$RNA@data, gset.idx.list = mapdgenesets, method = "gsva", parallel.sz = 6)
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
  
  ACC.list <- c()
  for (method in methods){
    ACC.list <- c(ACC.list, ACC.calc(ACC.df = method.vec[[method]]))
  }
  ACC.df <- data.frame(ACC.score = ACC.list)
  ACC.df$methods <- methods
  return(ACC.df)
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

ACC.calc <- function(ACC.df){
  ACC.rate <- sum(ACC.df$max == ACC.df$cell.type)/nrow(ACC.df)
  return(ACC.rate)
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
    markers <- sum$target.list[seq(0, 3*(n.cell.types - 1), by = 3) + p]
    tar.type.name <- tar.type[p]
    names(markers) <- clist
    # print(markers)
    ACC.batch <- list()
    AUC.batch <- list()
    for (k in 1:k.range){
      ds.markers <- list()
      for (q in 1:n.cell.types){
        ds.markers[[q]] <- ds.glist(glist = markers[[q]], all.genes = rownames(sr), fixed = fixed.ds, s = k)
      }
      #ds.markers[[1]] <- ds.glist(glist = markers[[1]], all.genes = rownames(sr), fixed = 20, s = k)
      #ds.markers[[2]] <- ds.glist(glist = markers[[2]], all.genes = rownames(sr), fixed = 20, s = k)
      #ds.markers[[3]] <- ds.glist(glist = markers[[3]], all.genes = rownames(sr), fixed = 20, s = k)
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


print.auc <- function(design, drop.out){
  info1 <- readRDS(paste0(design, ".1.", drop.out, ".info.rds"))
  info2 <- readRDS(paste0(design, ".2.", drop.out, ".info.rds"))
  info3 <- readRDS(paste0(design, ".3.", drop.out, ".info.rds"))
  
  df1 <- do.call(rbind, info1$AUC)
  df2 <- do.call(rbind, info2$AUC)
  df3 <- do.call(rbind, info3$AUC)
  
  plt.df <- do.call(rbind, list(df1, df2, df3))
  plt.df$target.type <- factor(plt.df$target.type , levels=c("strong", "medium", "weak"))
  
  pdf(paste0(design, ".", drop.out, ".AUC.pdf"), width = 8, height = 6)
  p <- plot.AUC(plt.df)
  print(p)
  dev.off()
}

plot.AUC <- function(AUC.df){
  p <- ggplot(data = AUC.df, aes(x = methods, y = AUC.score, fill = target.type))+
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_brewer(palette="BuPu") +
    geom_jitter(data = AUC.df, aes(x = methods, y = AUC.score, alpha = 0.1), #, color = cell.type
                shape = 2, size = 0.5) + 
    theme_classic() + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 1, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15))
  return(p)
}


output.auc <- function(design, drop.out){
  info1 <- readRDS(paste0(design, ".1.", drop.out, ".info.rds"))
  info2 <- readRDS(paste0(design, ".2.", drop.out, ".info.rds"))
  info3 <- readRDS(paste0(design, ".3.", drop.out, ".info.rds"))
  
  df1 <- do.call(rbind, info1$AUC)
  df2 <- do.call(rbind, info2$AUC)
  df3 <- do.call(rbind, info3$AUC)
  
  plt.df <- do.call(rbind, list(df1, df2, df3))
  plt.df$target.type <- factor(plt.df$target.type , levels=c("strong", "medium", "weak"))
  out.df <- plt.df %>%
    group_by(methods, target.type) %>%
    summarise_at(vars(AUC.score), list(name = mean)) %>% 
    as.data.frame()
  out.df$dp <- drop.out
  out.df$target.dp <- paste0(out.df$dp, "_", out.df$target.type)
  return(out.df)
}

output.acc <- function(design, drop.out){
  info1 <- readRDS(paste0(design, ".1.", drop.out, ".info.rds"))
  info2 <- readRDS(paste0(design, ".2.", drop.out, ".info.rds"))
  info3 <- readRDS(paste0(design, ".3.", drop.out, ".info.rds"))
  
  df1 <- do.call(rbind, info1$ACC)
  df2 <- do.call(rbind, info2$ACC)
  df3 <- do.call(rbind, info3$ACC)
  
  plt.df <- do.call(rbind, list(df1, df2, df3))
  plt.df$target.type <- factor(plt.df$target.type , levels=c("strong", "medium", "weak"))
  out.df <- plt.df %>%
    group_by(methods, target.type) %>%
    summarise_at(vars(ACC.score), list(name = mean)) %>% 
    as.data.frame()
  out.df$dp <- drop.out
  out.df$target.dp <- paste0(out.df$dp, "_", out.df$target.type)
  return(out.df)
}
make.bubble.auc2 <- function(design){
  auc.df <- do.call(rbind, list(output.auc(design = design, drop.out = "dyn.dp1"), 
                                output.auc(design = design, drop.out = "dyn.dp2"), 
                                output.auc(design = design, drop.out = "dyn.dp3")))
  auc.df$target.dp <- factor(auc.df$target.dp, levels=c("dyn.dp1_strong", "dyn.dp1_medium", "dyn.dp1_weak", 
                                                        "dyn.dp2_strong", "dyn.dp2_medium", "dyn.dp2_weak",
                                                        "dyn.dp3_strong", "dyn.dp3_medium", "dyn.dp3_weak"))
  p <- ggplot(auc.df, aes(x = target.dp, y = methods, label = name)) + 
    geom_point(aes(size = name, fill = name), shape = 21) +
    theme_classic() + 
    scale_x_discrete(position = "top") +
    labs( x= "", y = "", size = "AUC score", fill = "AUC score") + 
    scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), limits=c(0, 1)) + 
    guides(fill = guide_colorbar(order = 2),size = guide_legend(order = 1)) + 
    scale_size_continuous(range = c(0, 10)) + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 0, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2))
  pdf(paste0("sum.", design, ".bubble.AUC.pdf"), width = 7, height = 5)
  print(p)
  dev.off()
}

make.bubble.auc <- function(design){
  auc.df <- do.call(rbind, list(output.auc(design = design, drop.out = "dp1"), 
                                output.auc(design = design, drop.out = "dp2"), 
                                output.auc(design = design, drop.out = "dp3")))
  auc.df$target.dp <- factor(auc.df$target.dp, levels=c("dp1_strong", "dp1_medium", "dp1_weak", 
                                                        "dp2_strong", "dp2_medium", "dp2_weak",
                                                        "dp3_strong", "dp3_medium", "dp3_weak"))
  p <- ggplot(auc.df, aes(x = target.dp, y = methods, label = name)) + 
    geom_point(aes(size = name, fill = name), shape = 21) +
    theme_classic() + 
    scale_x_discrete(position = "top") +
    labs( x= "", y = "", size = "AUC score", fill = "AUC score") + 
    scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), limits=c(0, 1)) + 
    guides(fill = guide_colorbar(order = 2),size = guide_legend(order = 1)) + 
    scale_size_continuous(range = c(0, 10)) + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 0, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2))
  pdf(paste0("sum.", design, ".bubble.AUC.pdf"), width = 7, height = 5)
  print(p)
  dev.off()
}
make.bubble.acc2 <- function(design){
  acc.df <- do.call(rbind, list(output.acc(design = design, drop.out = "dyn.dp1"), 
                                output.acc(design = design, drop.out = "dyn.dp2"), 
                                output.acc(design = design, drop.out = "dyn.dp3")))
  acc.df$target.dp <- factor(acc.df$target.dp, levels=c("dyn.dp1_strong", "dyn.dp1_medium", "dyn.dp1_weak", 
                                                        "dyn.dp2_strong", "dyn.dp2_medium", "dyn.dp2_weak",
                                                        "dyn.dp3_strong", "dyn.dp3_medium", "dyn.dp3_weak"))
  p <- ggplot(acc.df, aes(x = target.dp, y = methods, label = name)) + 
    geom_point(aes(size = name, fill = name), shape = 21) +
    theme_classic() + 
    scale_x_discrete(position = "top") +
    labs( x= "", y = "", size = "ACC score", fill = "ACC score") + 
    scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), limits=c(0, 1)) + 
    scale_size_continuous(range = c(0, 10)) + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 0, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2))
  pdf(paste0("sum.", design, ".bubble.ACC.pdf"), width = 7, height = 5)
  print(p)
  dev.off()
}

make.bubble.acc <- function(design){
  acc.df <- do.call(rbind, list(output.acc(design = design, drop.out = "dp1"), 
                                output.acc(design = design, drop.out = "dp2"), 
                                output.acc(design = design, drop.out = "dp3")))
  acc.df$target.dp <- factor(acc.df$target.dp, levels=c("dp1_strong", "dp1_medium", "dp1_weak", 
                                                        "dp2_strong", "dp2_medium", "dp2_weak",
                                                        "dp3_strong", "dp3_medium", "dp3_weak"))
  p <- ggplot(acc.df, aes(x = target.dp, y = methods, label = name)) + 
    geom_point(aes(size = name, fill = name), shape = 21) +
    theme_classic() + 
    scale_x_discrete(position = "top") +
    labs( x= "", y = "", size = "ACC score", fill = "ACC score") + 
    scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), limits=c(0, 1)) + 
    scale_size_continuous(range = c(0, 10)) + 
    theme(axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 0, size = 15),
          axis.text.y = element_text(color = "black", size = 15),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size = 15),
          panel.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2))
  pdf(paste0("sum.", design, ".bubble.ACC.pdf"), width = 7, height = 5)
  print(p)
  dev.off()
}