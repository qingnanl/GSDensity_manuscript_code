library(Seurat)
library(copykat)
library(dplyr)
library(gsdensity)
library(ggplot2) # for plotting
library(reshape2)
library(msigdbr)
library(pheatmap)
library(ggrepel)
tnbc1 <- read.delim("GSM4476486_combined_UMIcount_CellTypes_TNBC1.txt")

meta <- as.data.frame(t(tnbc1[1:2, ]))

tnbc1.mat <- tnbc1[3:nrow(tnbc1), ]
tnbc1.mat <- as.matrix(tnbc1.mat)
tnbc1.mat <- matrix(as.numeric(tnbc1.mat),    # Convert to numeric matrix
                  ncol = ncol(tnbc1.mat))
rownames(tnbc1.mat) <- rownames(tnbc1)[3:nrow(tnbc1)]
colnames(tnbc1.mat) <- colnames(tnbc1)
tnbc1.obj <- CreateSeuratObject(counts = tnbc1.mat)

exp.rawdata <- as.matrix(tnbc1.obj@assays$RNA@counts)
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, 
                        sam.name="test", distance="euclidean", norm.cell.names="",output.seg="FLASE", 
                        plot.genes="TRUE", genome="hg20",n.cores=1)
pred.test <- data.frame(copykat.test$prediction)
# pred.test <- pred.test[pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
CNA.test <- data.frame(copykat.test$CNAmat)


my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

#jpeg("heatmap.tnbc1.pred.jpeg", width = 6, height = 4)
#heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
#          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
#          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
#          keysize=1, density.info="none", trace="none",
#          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
#          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
#
#legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
#dev.off()

tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =1, method = "euclidean"), method = "ward.D2")
hc.umap <- cutree(hcc,4)

rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:4])
subpop <- rbPal6(4)[as.numeric(factor(hc.umap))]
cells <- rbind(subpop,subpop)


jpeg("heatmap.tnbc1.pred.jpeg", width = 6000, height = 3000, res = 300)
heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=as.matrix(chr1[, 1]),RowSideColors=t(as.matrix(cells[1, ])),Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", c("c1","c2", "c3", "c4"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:4], cex=0.9, bty='n')
dev.off()

tnbc1.tu <- tnbc1.mat[, tumor.cells]
tnbc1.tu.obj <- CreateSeuratObject(tnbc1.tu)
tnbc1.tu.obj$cluster <- paste0("cluster_", as.numeric(factor(hc.umap)))
tnbc1.tu.obj <- tnbc1.tu.obj %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:10)


pdf("cluster.tnbc1.pdf", width = 5, height = 3)
DimPlot(tnbc1.tu.obj, group.by = "cluster", cols = rbPal6(4))
dev.off()

saveRDS(tnbc1.tu.obj, "tnbc1.tu.obj.rds")

tnbc1.tu.obj <- readRDS("tnbc1.tu.obj.rds")
ce <- compute.mca(object = tnbc1.tu.obj)
cells <- colnames(tnbc1.tu.obj)
el <- compute.nn.edges(coembed = ce, nn.use = 100)

mdb_h <- msigdbr(species = "Homo sapiens", category = "H")

gene.set.list <- list()
for (gene.set.name in unique(mdb_h$gs_name)){
  gene.set.list[[gene.set.name]] <- mdb_h[mdb_h$gs_name %in% gene.set.name, ]$gene_symbol
}
genes <- sapply(gene.set.list, 
                function(x)paste(x, collapse = ", "))
gene.set.list.df <- cbind(gene.set = names(gene.set.list), genes = genes)
rownames(gene.set.list.df) <- 1:nrow(gene.set.list.df)
write.csv(gene.set.list.df, "hallmarkr.list.info.csv")


res <- compute.kld(coembed = ce, genes.use = intersect(rownames(ce), rownames(tnbc1.tu.obj)), 
                   gene.set.list = gene.set.list )

cv.df <- run.rwr.list(el = el, gene_set_list = gene.set.list, cells = cells)
cl <- compute.cell.label.df(cv.df)

tnbc1.tu.obj@meta.data <- cbind(tnbc1.tu.obj@meta.data, cl[colnames(tnbc1.tu.obj), ])

DimPlot(tnbc1.tu.obj,
        group.by = c("HALLMARK_G2M_CHECKPOINT"),
        cols = c("grey", "blue"),
        raster = F)

library(cluster)

dist.matrix <- dist(x = Embeddings(object = tnbc1.tu.obj[["pca"]])[, 1:10])
sil.store <- c()
n.clusters <- c()
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2)){
  tnbc1.tu.obj <- FindClusters(object = tnbc1.tu.obj, resolution = res)
  clusters <- tnbc1.tu.obj$seurat_clusters
  n.clusters <- c(n.clusters, length(table(clusters)))
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  sil.store <- c(sil.store, mean(sil[, 3]))
}

sil.df <- data.frame(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2), 
                     silhouette = sil.store, 
                     n.clusters = n.clusters)
p <- ggplot(data = sil.df, aes(x = resolution, y = silhouette, label = n.clusters)) +        
  geom_text_repel() + 
  geom_point(alpha = 0.5) + 
  theme_classic() + 
  xlab("Resolution") + ylab("Silhouette score") + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )
p
pdf("scan.resolution.tnbc1.pdf", width = 4, height = 3)
p
dev.off()
