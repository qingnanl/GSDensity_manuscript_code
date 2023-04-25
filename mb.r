library(gsdensity)
library(ggplot2) # for plotting
library(ggrepel)
library(reshape2)
library(msigdbr) # for gathering gene sets
library(Seurat)
library(SeuratData)
library(future) # for parallel computing
library(future.apply) # for parallel computing
library(matrixStats)

mdb_c5 <- msigdbr(species = "Mus musculus", category = "C5")

# If we just want to do biological process:
mdb_c5_bp <- mdb_c5[mdb_c5$gs_subcat == "GO:BP", ]

# convert msigdbr gene sets to a list good for the input
gene.set.list <- list()
for (gene.set.name in unique(mdb_c5_bp$gs_name)){
  gene.set.list[[gene.set.name]] <- mdb_c5_bp[mdb_c5_bp$gs_name %in% gene.set.name, ]$gene_symbol
}
genes <- sapply(gene.set.list, 
                function(x)paste(x, collapse = ", "))
gene.set.list.df <- cbind(gene.set = names(gene.set.list), genes = genes)
rownames(gene.set.list.df) <- 1:nrow(gene.set.list.df)
write.csv(gene.set.list.df, "gene.set.info.csv")

brain <- LoadData("stxBrain", type = "anterior1")
print(brain)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
# Dimensionality reduction, clustering, and visualization

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# plot the clustering information
pdf("brain.umap.pdf", width = 5, height = 5)
DimPlot(brain, reduction = "umap", label = TRUE)
dev.off()

pdf("spatial.brain.overview.pdf", width = 6, height = 5)
SpatialDimPlot(brain, label = TRUE, label.size = 5)
dev.off()

ce <- compute.mca(object = brain)
res <- compute.kld(coembed = ce, 
                   genes.use = rownames(brain), 
                   n.grids = 100, 
                   gene.set.list = gene.set.list,
                   gene.set.cutoff = 20,
                   n.times = 100)

# we will then focus on th deviated gene sets; here we set a more stringent alpha level cutoff
gene.set.deviated <- res[res$p.adj < 0.01, ]$gene.set
# length(gene.set.deviated)

# compute a nearest neighbor graph (edge list) in the MCA space
cells <- colnames(brain)
el <- compute.nn.edges(coembed = ce, nn.use = 300)

# We then compute the relevance between each cell and the deviated gene sets

cv.df <- run.rwr.list(el = el, gene_set_list = gene.set.list[gene.set.deviated], cells = cells)

cl.df <- compute.cell.label.df(cv.df)


# An optional filtering step: we want to only keep the terms with certain numbers of positive cells; here we use 100

positive.count <- apply(cl.df, MARGIN = 2, FUN = function(x) {length(x[x == "positive"])})
gene.set.deviated.2 <- names(positive.count[positive.count > 100])
coords.df <- brain@images$anterior1@coordinates[, c("imagerow", "imagecol")]
head(coords.df)
spatial.klds <- compute.spatial.kld.df(spatial.coords = coords.df, 
                                       weight_df = cv.df[, gene.set.deviated.2], 
                                       n = 10)

jsd.df <- compute.spec(cell_df = cv.df[, gene.set.deviated.2], 
                       metadata = brain@meta.data, # each row is a cell; columns include partition information
                       cell_group = "seurat_clusters" # 'cell_group' should use a column name in the metadata as the input
)

saveRDS(cv.df, "cv.df.rds")
saveRDS(cl.df, "cl.df.rds")
saveRDS(ce, "ce.rds")
saveRDS(spatial.klds, "spatial.klds.rds")
saveRDS(res, "res.rds")
saveRDS(jsd.df, "jsd.df.rds")
saveRDS(brain, "brain.rds")

identical(rownames(jsd.df), names(spatial.klds))

jsd.df <- as.data.frame(jsd.df)
jsd.df$max.specificity <- rowMaxs(as.matrix(jsd.df))

spat.clu <- data.frame(cluster_specificity = jsd.df$max.specificity, 
                       spatial_relevance = spatial.klds$spatial.kld, 
                       gene.set = rownames(spatial.klds))

hl1 <- spat.clu[spat.clu$cluster_specificity > quantile(spat.clu$cluster_specificity, 0.95), ]
hl2 <- spat.clu[spat.clu$cluster_specificity < quantile(spat.clu$cluster_specificity, 0.5) & 
                  spat.clu$spatial_relevance > quantile(spat.clu$spatial_relevance, 0.8),  ]

write.csv(hl1, "hl1.csv")
write.csv(hl2, "hl2.csv")

pdf("spatial.cluster.scatter.pdf", width = 4, height = 4)
ggplot(spat.clu, aes(x=cluster_specificity, y=spatial_relevance)) + 
  geom_point(shape = 1, alpha = 0.7) + 
  # scale_y_continuous(trans='log10') + 
  geom_point(shape = 1, data = hl1, color = 'red', alpha = 0.7) + 
  geom_point(shape = 1, data = hl2, color = 'darkblue', alpha = 0.7) + 
  geom_segment(aes(x = 0.12, xend = quantile(spat.clu$cluster_specificity, 0.5), 
                   y = quantile(spat.clu$spatial_relevance, 0.8), 
                   yend = quantile(spat.clu$spatial_relevance, 0.8)),
             linetype="dashed", color = "blue", size = 0.2) + 
  geom_vline(xintercept=quantile(spat.clu$cluster_specificity, 0.5), 
             linetype="dashed", color = "blue", size = 0.2) + 
  geom_vline(xintercept=quantile(spat.clu$cluster_specificity, 0.95), 
             linetype="dashed", color = "red", size = 0.2) + 
  #geom_text_repel(aes(label = ifelse(term %in% highlighted$term, term, "")), 
  #                size = 2) + 
  xlab("Maximum cluster-wise specificity") + 
  ylab("Spatial relevance of gene set") + 
  theme_classic() + 
  theme(axis.title=element_text(size=14,face="bold"))
  
dev.off()
