library(Seurat)
library(gsdensity)
library(ggplot2) # for plotting
library(ggrepel)
library(reshape2)
library(msigdbr) # for gathering gene sets
# library(copykat)
library(dplyr)
library(future) # for parallel computing
library(future.apply) # for parallel computing
plan(multiprocess, workers = 20)

# 1. get gene-sets
mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
mdb_c5_bp <- mdb_c5[mdb_c5$gs_subcat == "GO:BP", ]
mdb_h <- msigdbr(species = "Homo sapiens", category = "H")
mdb_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
mdb_c2_cp <- mdb_c2[grep("CP", mdb_c2$gs_subcat), ]
gl.dataframe <- Reduce(rbind, list(mdb_c5_bp, mdb_h, mdb_c2_cp))

gene.set.list <- list()
for (gene.set.name in unique(gl.dataframe$gs_name)){
  gene.set.list[[gene.set.name]] <- gl.dataframe[gl.dataframe$gs_name %in% gene.set.name, ]$gene_symbol
}
print(length(gene.set.list))
saveRDS(gene.set.list, "gene.set.list.rds")

genes <- sapply(gene.set.list, 
                function(x)paste(x, collapse = ", "))
gene.set.list.df <- cbind(gene.set = names(gene.set.list), genes = genes)
rownames(gene.set.list.df) <- 1:nrow(gene.set.list.df)
write.csv(gene.set.list.df, "gene.set.list.df.csv")

# 
gene.set.list <- readRDS("gene.set.list.rds")

# get data
info <- read.csv("sample_info.csv")
info
for (i in 1:nrow(info)){
        sample.name <- info[i, "sample_name"]
        img <- Read10X_Image(image.dir = info[i, "imgpwd"], 
                             image.name = "tissue_lowres_image.png", filter.matrix = TRUE)
        data <- Load10X_Spatial(data.dir = info[i, "pwd"], 
                                info[i, "matrix_name"], 
                                assay = "RNA", 
                                slice = "slice1", 
                                filter.matrix = TRUE, 
                                to.upper = FALSE, 
                                image = img)
        data <- SCTransform(data, assay = "RNA", verbose = FALSE) %>%
                RunPCA(assay = "SCT", verbose = FALSE) %>%
                FindNeighbors(dims = 1:30) %>%
                FindClusters() %>%
                RunUMAP(dims = 1:30)
        ce <- compute.mca(object = data)
        res <- compute.kld(coembed = ce, 
                   genes.use = intersect(rownames(data), rownames(ce)), 
                   n.grids = 100, 
                   gene.set.list = gene.set.list,
                   gene.set.cutoff = 10,
                   n.times = 100)
        write.csv(res, paste0(sample.name, "_res.csv"))
        gene.set.deviated <- res[res$p.adj < 0.01, ]$gene.set
        cells <- colnames(data)
        el <- compute.nn.edges(coembed = ce, nn.use = 300)
        cv.df <- run.rwr.list(el = el, 
                              gene_set_list = gene.set.list[gene.set.deviated], 
                              cells = cells)
        cl.df <- compute.cell.label.df(cv.df)
        positive.count <- apply(cl.df, MARGIN = 2, 
                                FUN = function(x) {length(x[x == "positive"])})
        gene.set.deviated.2 <- names(positive.count[positive.count > 50])
        coords.df <- data@images$slice1@coordinates[, c("imagerow", "imagecol")]
        spatial.klds <- compute.spatial.kld.df(spatial.coords = coords.df, 
                                       weight_df = cv.df[, gene.set.deviated.2], 
                                       n = 10, n.times = 50)
        jsd.df <- compute.spec(cell_df = cv.df[, gene.set.deviated.2], 
                       metadata = data@meta.data, 
                       cell_group = "seurat_clusters")
        write.csv(spatial.klds, paste0(sample.name, ".spatial.klds.csv"))
        data@meta.data <- cbind(data@meta.data, cv.df[rownames(data@meta.data), ])
        saveRDS(cv.df, paste0(sample.name, ".cv.df.rds"))
        saveRDS(cl.df, paste0(sample.name, ".cl.df.rds"))
        saveRDS(jsd.df, paste0(sample.name, ".jsd.df.rds"))
        saveRDS(data, paste0(sample.name, ".data.rds"))
}
