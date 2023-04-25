library(dplyr)
library(monocle3)
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(msigdbr)
library(gsdensity)
library(TSrepr)
library(data.table)
library(cluster)
library(clusterCrit)
library(ape)

data <- Read10X_h5(filename = "./GSM5277844_E17_5_filtered_feature_bc_matrix.h5")
data <- CreateSeuratObject(counts = data, project='E17.5', min.cells = 3, min.features = 200)

data <- NormalizeData(data) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunUMAP(dims = 1:20)
Idents(data) <- "seurat_clusters"

pdf("data.umap.labeled.pdf", width = 5, height = 4)
DimPlot(data, label = T)
dev.off()

FeaturePlot(data, c("Eomes", "Neurod2"))

subdata <- subset(data, idents = c("0", "1", "2", "3", "6", "7", "9", "10", "11"))

pdf("subdata.umap.labeled.pdf", width = 5, height = 4)
DimPlot(subdata, label = T) + xlim(-10, 10) + ylim(-10, 10)
dev.off()

pdf("subdata.features.pdf", width = 10, height = 8)
FeaturePlot(subdata, c("Sox2", "Eomes", "Neurod2", "Fezf2"))
dev.off()

subdata <- NormalizeData(subdata)

#######################
# kegg gene set
mdb <- msigdbr(species = "Mus musculus", category = "C2")

# If we just want to do biological process:
mdb_kegg <- mdb[mdb$gs_subcat %in% c("CP:KEGG", "CP:BIOCARTA"), ]


# convert msigdbr gene sets to a list good for the input
gene.set.list <- list()
for (gene.set.name in unique(mdb_kegg$gs_name)){
  gene.set.list[[gene.set.name]] <- mdb_kegg[mdb_kegg$gs_name %in% gene.set.name, ]$gene_symbol
}

ce <- compute.mca(object = subdata)
res <- compute.kld(coembed = ce, 
                   genes.use = intersect(rownames(ce), rownames(subdata)), 
                   n.grids = 100, 
                   gene.set.list = gene.set.list,
                   gene.set.cutoff = 3,
                   n.times = 100)
gene.set.deviated <- res[res$p.adj < 0.05, ]$gene.set
cells <- colnames(subdata)
el <- compute.nn.edges(coembed = ce, nn.use = 300)
cv.df <- run.rwr.list(el = el, gene_set_list = gene.set.list[gene.set.deviated], cells = cells)
cl.df <- compute.cell.label.df(cv.df)
positive.count <- apply(cl.df, MARGIN = 2, FUN = function(x) {length(x[x == "positive"])})
gene.set.deviated.2 <- names(positive.count[positive.count > 100])


sr.mon <- as.cell_data_set(subdata)
sr.mon <- cluster_cells(cds = sr.mon, reduction_method = "UMAP")
sr.mon <- learn_graph(sr.mon, use_partition = TRUE)

roots <- colnames(subdata)[subdata$seurat_clusters == "6"]
sr.mon <- order_cells(sr.mon, reduction_method = "UMAP", root_cells = roots)

# plot trajectories colored by pseudotime

pdf("pseudotime.umap.pdf", width = 5, height = 4)
plot_cells(
  cds = sr.mon,
  color_cells_by = "pseudotime",
  show_trajectory_graph = F
) + xlim(-10, 10) + ylim(-10, 10)
dev.off()


subdata@meta.data <- cbind(subdata@meta.data, cv.df[rownames(subdata@meta.data), ])




colnames(sr.mon)
sr.mon@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]

pt.info <- as.data.frame(cbind(cells = colnames(sr.mon), 
                 pseudotime = sr.mon@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]))

pt.info$pseudotime <- as.numeric(pt.info$pseudotime)

# pt.info <- pt.info[order(pt.info$pseudotime, decreasing = F), ]
# cell.orders <- pt.info$cells

identical(rownames(pt.info), rownames(cv.df))
cv.df$time <- pt.info$pseudotime

library(plyr)
cv.df.bin <- t(ddply(cv.df, .(cut(cv.df$time, 20)), colwise(mean)))[2:98, ]
cv.df.bin <- apply(cv.df.bin, 2, as.numeric)
cv.df.bin <- scale(cv.df.bin)
# rownames(cv.df.bin) <- colnames(cv.df)[1:97]

data_gam_1 <- cv.df.bin

clusterings <- lapply(c(2:15), function(x)
  pam(data_gam_1, x))

DB_values <- sapply(seq_along(clusterings), function(x) 
  intCriteria(data_gam, as.integer(clusterings[[x]]$clustering),
              c("Davies_Bouldin")))
ggplot(data.table(Clusters = 2:15, DBindex = unlist(DB_values)),
       aes(Clusters, DBindex)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw()


# prepare data for plotting
data_plot <- data.table(melt(data.table(class = as.factor(clusterings[[7]]$clustering),
                                        data_gam_1)))
data_plot[, Time := rep(1:ncol(data_gam_1), each = nrow(data_gam_1))]
data_plot[, ID := rep(1:nrow(data_gam_1), ncol(data_gam_1))]

# prepare medoids
centers <- data.table(melt(clusterings[[7]]$medoids))
setnames(centers, c("Var1", "Var2"), c("class", "Time"))
centers[, ID := class]

# plot the results
pdf("trajectory.clustering.pdf", width = 6, height = 6)
ggplot(data_plot, aes(Time, value, group = ID)) +
  facet_wrap(~class, ncol = 2, scales = "free_y") +
  geom_line(color = "grey10", alpha = 0.25) +
  geom_line(data = centers, aes(Time, value),
            color = "firebrick1", alpha = 0.80, size = 1.2) +
  labs(x = "Partition of pseudotime", y = "Pathway relevance (scaled)") +
  theme_bw()
dev.off()


cluster.info <- as.data.frame(cbind(gene.sets = colnames(cv.df)[1:97], 
                      clustering = clusterings[[7]]$clustering))


saveRDS(ce, "ce.rds")
saveRDS(cluster.info, "cluster.info.rds")
saveRDS(cv.df, "cv.df.rds")
saveRDS(data, "data.rds")
saveRDS(data_gam, "data_gam.rds")
saveRDS(res, "res.rds")
saveRDS(sr.mon, "sr.mon.rds")
saveRDS(subdata, "subdata.rds")
