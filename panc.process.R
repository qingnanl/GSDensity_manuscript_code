library(Seurat)
library(ggplot2)
library(ggrepel)
library(UpSetR)
library(copykat)
library(gsdensity)
library(dior)

pc <- readRDS("visium_prostate_cancer.data.rds")
pacc <- readRDS("visium_prostate_acinar_cell_carcinoma.data.rds")
oc <- readRDS("visium_ovarian_cancer.data.rds")
ic <- readRDS("visium_intestinal_cancer.data.rds")
cc <- readRDS("visium_cervical_cancer.data.rds")
bc <- readRDS("visium_breast_cancer.data.rds")

runcopykat <- function(object, sample.name, nc){
  exp.rawdata <- as.matrix(object@assays$RNA@counts)
  copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", 
                          ngene.chr=5, win.size=25, KS.cut=0.1, sam.name= sample.name, 
                          distance="euclidean", norm.cell.names=nc,
                          output.seg="FLASE", plot.genes="FALSE", genome="hg20",n.cores=8)
  pred.test <- data.frame(copykat.test$prediction)
  # pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
  CNA.test <- data.frame(copykat.test$CNAmat)
  return(list(CNA.test, pred.test))
}
pc.nc <- colnames(pc)[pc$seurat_clusters == 8]
pc.cna <- runcopykat(object = pc, nc = pc.nc, sample.name = "PC")
saveRDS(pc.cna, "pc.cna.rds")

bc.nc <- colnames(bc)[bc$seurat_clusters == 9]
bc.cna <- runcopykat(object = bc, nc = bc.nc, sample.name = "BC")
saveRDS(bc.cna, "bc.cna.rds")

cc.nc <- colnames(cc)[cc$seurat_clusters == 11]
cc.cna <- runcopykat(object = cc, nc = cc.nc, sample.name = "CC")
saveRDS(cc.cna, "cc.cna.rds")

ic.nc <- colnames(ic)[ic$seurat_clusters == 4]
ic.cna <- runcopykat(object = ic, nc = ic.nc, sample.name = "IC")
saveRDS(ic.cna, "ic.cna.rds")

oc.nc <- colnames(oc)[oc$seurat_clusters == 1]
oc.cna <- runcopykat(object = oc, nc = oc.nc, sample.name = "OC")
saveRDS(oc.cna, "oc.cna.rds")

pacc.nc <- colnames(pacc)[pacc$seurat_clusters == 1]
pacc.cna <- runcopykat(object = pacc, nc = pacc.nc, sample.name = "PACC")
saveRDS(pacc.cna, "pacc.cna.rds")

pc$tumor.prediction <- ifelse(pc$prediction == "aneuploid", "tumor", "normal")
saveRDS(pc, "visium_prostate_cancer.data.rds")
bc$tumor.prediction <- ifelse(bc$prediction == "aneuploid", "tumor", "normal")
saveRDS(bc, "visium_breast_cancer.data.rds")
cc$tumor.prediction <- ifelse(cc$prediction == "aneuploid", "tumor", "normal")
saveRDS(cc, "visium_cervical_cancer.data.rds")
ic$tumor.prediction <- ifelse(ic$prediction == "aneuploid", "tumor", "normal")
saveRDS(ic, "visium_intestinal_cancer.data.rds")
oc$tumor.prediction <- ifelse(oc$prediction == "aneuploid", "tumor", "normal")
saveRDS(oc, "visium_ovarian_cancer.data.rds")
pacc$tumor.prediction <- ifelse(pacc$prediction == "aneuploid", "tumor", "normal")
saveRDS(pacc, "visium_prostate_acinar_cell_carcinoma.data.rds")


fl <- c("visium_prostate_cancer", "visium_prostate_acinar_cell_carcinoma", 
        "visium_ovarian_cancer", "visium_intestinal_cancer", "visium_cervical_cancer",
        "visium_breast_cancer")
basic.plot <- function(object, output.name){
  p1 <- DimPlot(object, reduction = "umap", label = TRUE)
  Idents(object) <- "seurat_clusters"
  p2 <- SpatialDimPlot(object, label = TRUE, label.size = 3)
  Idents(object) <- "tumor.prediction"
  p3 <- SpatialDimPlot(object, cells.highlight = colnames(object)[object$tumor.prediction == "tumor"]) + NoLegend()
  pdf(paste0("./plots/", output.name, ".overview.pdf"), width = 15, height = 5)
  print(p1 + p2 + p3)
  dev.off()
}



###########################
basic.plot(object = pc, output.name = "visium_prostate_cancer")
basic.plot(object = bc, output.name = "visium_breast_cancer")
basic.plot(object = cc, output.name = "visium_cervical_cancer")
basic.plot(object = ic, output.name = "visium_intestinal_cancer")
basic.plot(object = oc, output.name = "visium_ovarian_cancer")
basic.plot(object = pacc, output.name = "visium_prostate_acinar_cell_carcinoma_cancer")


plot.overlay.density <- function(object, term, output.name){
  scf <- object@images[["slice1"]]@scale.factors[["lowres"]]
  object.coords <- object@images$slice1@coordinates[, c("imagerow", "imagecol")]
  object.coords$imagerow <- object.coords$imagerow * scf
  object.coords$imagecol <- object.coords$imagecol * scf
  y.max <- max(object.coords$imagerow)
  y.min <- min(object.coords$imagerow)
  
  
  object.tumor.coords <- object@images$slice1@coordinates[, c("imagerow", "imagecol")][object$prediction == "aneuploid", ]
  object.tumor.coords$imagerow <- object.tumor.coords$imagerow * scf
  object.tumor.coords$imagecol <- object.tumor.coords$imagecol * scf
  object.tumor.coords$imagerow <- y.max - object.tumor.coords$imagerow + y.min
  
  pdf(paste0(output.name, "_", term, ".point.pdf"), width = 5, height = 5)
  p <- SpatialFeaturePlot(object, features = term, alpha = c(1)) + 
    # theme(legend.position = "top") + 
    geom_point(data = object.tumor.coords, 
               aes(x = imagecol, y = imagerow), 
               color = "red", shape = 2, size = 1, inherit.aes = F) 
  print(p)
  dev.off()
  
  pdf(paste0(output.name, "_", term, ".density.pdf"), width = 5, height = 5)
  p <- SpatialFeaturePlot(object, features = term, alpha = c(1)) + 
    # theme(legend.position = "top") + 
    geom_density_2d(data = object.tumor.coords, 
                    aes(x = imagecol, y = imagerow), 
                    color = "red",alpha = 0.8, bins = 5,
                    adjust = 0.2, size = 0.5,
                    inherit.aes = F)
  print(p)
  dev.off()
}

plot.overlay.density(object = pc, term = "GOBP_MESENCHYME_MORPHOGENESIS", output.name = "PC")
plot.overlay.density(object = bc, term = "GOBP_MESENCHYME_MORPHOGENESIS", output.name = "BC")
plot.overlay.density(object = cc, term = "GOBP_MESENCHYME_MORPHOGENESIS", output.name = "CC")
plot.overlay.density(object = ic, term = "GOBP_MESENCHYME_MORPHOGENESIS", output.name = "IC")
plot.overlay.density(object = oc, term = "GOBP_MESENCHYME_MORPHOGENESIS", output.name = "OC")
plot.overlay.density(object = pacc, term = "GOBP_MESENCHYME_MORPHOGENESIS", output.name = "PACC")

plot.overlay.density(object = pc, term = "BIOCARTA_GRANULOCYTES_PATHWAY", output.name = "PC")
plot.overlay.density(object = bc, term = "BIOCARTA_GRANULOCYTES_PATHWAY", output.name = "BC")
plot.overlay.density(object = cc, term = "BIOCARTA_GRANULOCYTES_PATHWAY", output.name = "CC")
plot.overlay.density(object = ic, term = "BIOCARTA_GRANULOCYTES_PATHWAY", output.name = "IC")
plot.overlay.density(object = oc, term = "BIOCARTA_GRANULOCYTES_PATHWAY", output.name = "OC")
plot.overlay.density(object = pacc, term = "BIOCARTA_GRANULOCYTES_PATHWAY", output.name = "PACC")

