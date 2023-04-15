library(Seurat)
library(SeuratData)
library(dplyr)
library(qpcR)
library(scRNAseq)
library(stringr)
sparsity <- function(matrix){
  sparsity <- sum(matrix == 0)/(dim(matrix)[1] * dim(matrix)[2])
  return(sparsity)
}
pang.markers <- as.data.frame(read.delim("PanglaoDB_markers_27_Mar_2020.tsv"))
cm2 <- read.csv("Cell_marker_Human.csv")
cm2m <- read.csv("Cell_marker_Mouse.csv")

data("pbmc3k")
pbmc3k.meta <- pbmc3k@meta.data

pbmc3k <- CreateSeuratObject(pbmc3k@assays$RNA@counts, min.cells = 10)
pbmc3k@meta.data <- pbmc3k.meta
Idents(pbmc3k) <- "seurat_annotations"
# remove cells annotated with NA; 
pbmc3k <- subset(pbmc3k, idents = c("Memory CD4 T", 
                                    "B",  
                                    "CD14+ Mono", 
                                    "NK",  "CD8 T", 
                                    "Naive CD4 T", 
                                    "FCGR3A+ Mono", 
                                    "DC", 
                                    "Platelet"))
# merge CD4 T cells to match marker list
pbmc3k$seurat_annotations <- recode(pbmc3k$seurat_annotations, "Naive CD4 T" = "CD4 T", "Memory CD4 T" = "CD4 T")
pbmc3k$cell.type <- pbmc3k$seurat_annotations
pbmc.markers <- read.csv("st4.csv")

# create marker list
pbmc.markers.list <- list("CD4 T" = intersect(pbmc.markers$CD4.T.cells, rownames(pbmc3k)), 
                          "B" = intersect(pbmc.markers$B.cells, rownames(pbmc3k)), 
                          "CD14+ Mono" = intersect(pbmc.markers$CD14.Monocytes, rownames(pbmc3k)), 
                          "NK" = intersect(pbmc.markers$NK.cells, rownames(pbmc3k)), 
                          "CD8 T" = intersect(pbmc.markers$CD8.T.cells, rownames(pbmc3k)), 
                          "FCGR3A+ Mono" = intersect(pbmc.markers$CD16.Monocytes, rownames(pbmc3k)), 
                          "DC" = intersect(pbmc.markers$DC, rownames(pbmc3k)), 
                          "Platelet" = intersect(pbmc.markers$Platelets, rownames(pbmc3k)))

# additional DC markers

DC.markers <- pang.markers[pang.markers$species == "Hs" & pang.markers$cell.type == "Dendritic cells", ]

pbmc.markers.list[["DC"]] <- intersect(unique(c(pbmc.markers.list[["DC"]], DC.markers$official.gene.symbol)), rownames(pbmc3k))

saveRDS(pbmc3k, "pbmc3k.rds")
saveRDS(pbmc.markers.list, "pbmc3k.markers.rds")
write.table(do.call(qpcR:::cbind.na, pbmc.markers.list), "pbmc3k.markers.csv", sep = ",")
print(length(unique(pbmc3k$cell.type)))
print(dim(pbmc3k))
sparsity(pbmc3k@assays$RNA@counts)
####
# InstallData("bmcite")
data("bmcite")

bmcite.meta <- bmcite@meta.data

bmcite <- CreateSeuratObject(bmcite@assays$RNA@counts, min.cells = 10)
bmcite@meta.data <- bmcite.meta
bmcite.small.markers <- read.csv("st4.csv")
# remove cells annotated with NA; 
# merge CD4 T cells to match marker list
bmcite$celltype.l2 <- recode(bmcite$celltype.l2, "CD4 Memory" = "CD4 T", "CD4 Naive" = "CD4 T", 
                             "CD8 Effector_1" = "CD8 T", "CD8 Effector_2" = "CD8 T", "CD8 Memory_1" = "CD8 T", 
                             "CD8 Memory_2" = "CD8 T", "CD8 Naive" = "CD8 T", "cDC2" = "cDC", "Memory B" = "B", 
                             "Naive B" = "B")
Idents(bmcite) <- "celltype.l2"
bmcite <- subset(bmcite, idents = c("CD4 T", 
                                    "B",  
                                    "CD14 Mono", 
                                    "NK",  "CD8 T", 
                                    "cDC", "GMP", "HSC", "pDC",
                                    "CD16 Mono"))
table(bmcite$celltype.l2)
bmcite.small <- subset(bmcite, downsample = 1500)
#table(bmcite.small@active.ident)

bmcite.small$cell.type <- as.character(bmcite.small@active.ident)

bmcite.small.markers.list <- list("CD4 T" = intersect(bmcite.small.markers$CD4.T.cells, rownames(bmcite.small)), 
                            "B" = intersect(bmcite.small.markers$B.cells, rownames(bmcite.small)), 
                            "CD14 Mono" = intersect(bmcite.small.markers$CD14.Monocytes, rownames(bmcite.small)), 
                            "NK" = intersect(bmcite.small.markers$NK.cells, rownames(bmcite.small)), 
                            "CD8 T" = intersect(bmcite.small.markers$CD8.T.cells, rownames(bmcite.small)), 
                            "CD16 Mono" = intersect(bmcite.small.markers$CD16.Monocytes, rownames(bmcite.small)), 
                            "cDC" = intersect(bmcite.small.markers$cDC, rownames(bmcite.small)), 
                            "GMP" = intersect(bmcite.small.markers$GMP, rownames(bmcite.small)),
                            "HSC" = intersect(bmcite.small.markers$HSC, rownames(bmcite.small)), 
                            "pDC" = intersect(bmcite.small.markers$pDC, rownames(bmcite.small)))

pDC.markers <- pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Plasmacytoid dendritic cells", ]

bmcite.small.markers.list[["pDC"]] <- intersect(unique(c(bmcite.small.markers.list[["pDC"]], pDC.markers$official.gene.symbol)), rownames(bmcite.small))

bmcite.small@assays$ADT <- NULL

saveRDS(bmcite.small, "bmcite.small.rds")
saveRDS(bmcite.small.markers.list, "bmcite.small.markers.rds")
write.table(do.call(qpcR:::cbind.na, bmcite.small.markers.list), "bmcite.small.markers.csv", sep = ",")

print(length(unique(bmcite.small$cell.type)))
print(dim(bmcite.small))
sparsity(bmcite.small@assays$RNA@counts)


##
#####
load("GSE183852_DCM_Cells")
heart <- subset(HDCM, subset = Condition == "Donor")
# DimPlot(heart, group.by = "Names")
Idents(heart) <- "Names"
heart <- subset(heart, downsample = 1000)
heart.meta <- heart@meta.data

heart <- CreateSeuratObject(heart@assays$RNA@counts, min.cells = 10)
heart@meta.data <- heart.meta

heart <- heart %>%
  NormalizeData(verbose = F) %>%
  FindVariableFeatures(verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = F) %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:10, verbose = F)

DimPlot(heart, group.by = "Names", label = T)

# str(bmcite)

heart$cell.type <- as.character(heart$Names)
table(heart$cell.type)

heart$cell.type <- recode(heart$cell.type, "NKT cell" = "NK cell", 
                            "Thick ascending limb of Loop of Henle" = "Loop of Henle") 
Idents(heart) <- "cell.type"


heart.markers.list <- list("Endothelium" = intersect(cm2[cm2$cell_name == "Endothelial cell" & cm2$tissue_class == "Heart", "Symbol"], 
                                                     rownames(heart)),
                           "Fibroblasts" = intersect(cm2[cm2$cell_name == "Fibroblast" & cm2$tissue_class == "Heart", "Symbol"], 
                                                     rownames(heart)),
                           "Macrophages" = intersect(cm2[cm2$cell_name == "Macrophage" & cm2$tissue_class == "Heart", "Symbol"], 
                                                     rownames(heart)),
                           "Smooth_Muscle" = intersect(cm2[cm2$cell_name == "Smooth muscle cell" & cm2$tissue_class == "Heart", "Symbol"], 
                                                     rownames(heart)),
                           "T_Cells" = intersect(cm2[cm2$cell_name == "T cell" & cm2$tissue_class == "Heart", "Symbol"], 
                                                     rownames(heart)),
                           "NK_Cells" = intersect(bmcite.markers$NK.cells, #here the bmcite markers is just XCell, not dedicated to bmcite
                                                      rownames(heart)), 
                             "Monocytes" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Monocytes", "official.gene.symbol"], 
                                                             rownames(heart)) 
                             )
# supply with
heart.markers.list$Macrophages <- unique(c(heart.markers.list$Macrophages, intersect(rownames(heart), 
                                                                     c("F13A1", "RBPJ", "MRC1", "CD163", "COLEC12", "RBM47", "FMN1", "MS4A6A", "STAB1", "FRMD4B"))))
heart.markers.list$Smooth_Muscle <- unique(c(heart.markers.list$Smooth_Muscle, intersect(rownames(heart), 
                                                                     c("MYH11", "ITGA8", "ACTA2", "CARMN", "KCNAB1", "TAGLN", "ZFHX3", "PRKG1", "NTRK3", "RCAN2"))))

heart.markers.list$T_Cells <- unique(c(heart.markers.list$T_Cells, intersect(rownames(heart), 
                                                                c("IL7R", "THEMIS", "ITK", "PARP8", "CAMK4", "SKAP1", "TC2N", "BCL11B", "PTPRC", "PRKCQ"))))



heart <- subset(heart, idents = names(heart.markers.list))


#table(heart@active.ident)
print(length(unique(heart$cell.type)))
print(dim(heart))
sparsity(heart@assays$RNA@counts)
DimPlot(heart, group.by = "cell.type", label = T)


saveRDS(heart, "heart.rds")
saveRDS(heart.markers.list, "heart.markers.rds")
write.table(do.call(qpcR:::cbind.na, heart.markers.list), "heart.markers.csv", sep = ",")



######
out <- listDatasets()
out

pancrea <- BaronPancreasData('human')
pancrea <- as.Seurat(pancrea, data = NULL)

pancrea.meta <- pancrea@meta.data

pancrea <- CreateSeuratObject(pancrea@assays$originalexp@counts, min.cells = 10)
pancrea@meta.data <- pancrea.meta

table(pancrea$donor)
Idents(pancrea) <- "donor"
pancrea <- subset(pancrea, idents = "GSM2230759")

pancrea$cell.type <- as.character(pancrea$label)
Idents(pancrea) <- "cell.type"

pancrea.markers.list <- list("acinar" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Acinar cells", "official.gene.symbol"], 
                                                 rownames(pancrea)), 
                            "alpha" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Alpha cells", "official.gene.symbol"], 
                                            rownames(pancrea)), 
                            "beta" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Beta cells", "official.gene.symbol"], 
                                               rownames(pancrea)), 
                            "delta" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Delta cells", "official.gene.symbol"], 
                                             rownames(pancrea)), 
                            "ductal" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Ductal cells", "official.gene.symbol"], 
                                                rownames(pancrea)), 
                            "epsilon" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Epsilon cells", "official.gene.symbol"], 
                                                    rownames(pancrea)), 
                            "gamma" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Gamma (PP) cells", "official.gene.symbol"], 
                                              rownames(pancrea)), 
                            "quiescent_stellate" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Pancreatic stellate cells", "official.gene.symbol"], 
                                              rownames(pancrea)),
                            "schwann" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Peri-islet Schwann cells", "official.gene.symbol"], 
                                                  rownames(pancrea)))
pancrea <- subset(pancrea, idents = names(pancrea.markers.list))

print(length(unique(pancrea$cell.type)))
print(dim(pancrea))
sparsity(pancrea@assays$RNA@counts)


pancrea <- readRDS("pancrea.rds")
pancrea <- pancrea %>%
  NormalizeData(verbose = F) %>%
  FindVariableFeatures(verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = F) %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:10, verbose = F)
DimPlot(pancrea, group.by = "cell.type", label = T)
# DimPlot(pancrea, group.by = "donor", label = T)

saveRDS(pancrea, "pancrea.rds")
saveRDS(pancrea.markers.list, "pancrea.markers.rds")
write.table(do.call(qpcR:::cbind.na, pancrea.markers.list), "pancrea.markers.csv", sep = ",")

##
liver.imm <- ZhaoImmuneLiverData()
liver.imm <- as.Seurat(liver.imm, data = NULL)

Idents(liver.imm) <- "sample"

liver.immune <- subset(liver.imm, idents = "donor 1 liver")

liver.immune.meta <- liver.immune@meta.data
matx <- liver.immune@assays$originalexp@counts
rownames(matx) <- liver.imm@assays$originalexp@meta.features$Symbol
liver.immune <- CreateSeuratObject(matx, min.cells = 10)
liver.immune@meta.data <- liver.immune.meta

liver.immune$cell.type <- as.character(liver.immune$fine)
Idents(liver.immune) <- "cell.type"
liver.immune <- liver.immune %>%
  NormalizeData(verbose = F) %>%
  FindVariableFeatures(verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = F) %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:10, verbose = F)
DimPlot(liver.immune, group.by = "cell.type", label = T)
table(liver.immune$cell.type)

liver.immune$cell.type <- recode(liver.immune$cell.type, "CCR7+ CD4" = "CD4", "CCR7+ CD8" = "CD8","CRCR3+ CD8" = "CD8","CRCR3+ CD4" = "CD4",
                                 "CX3CR1+ CD8" = "CD8", "CXCR6+ CD8" = "CD8","Naive B" = "Naive B", "cycling ASC" = "Plasma", "noncycline ASC" = "Plasma",
                                 "CD11c+ mem B" = "mem B","Classical mem B" = "mem B","CX3CR1+ NK" = "NK","CXCR6+ NK" = "NK", "Treg" = "Treg",
                                 "CD14+ Mo" = "CD14+ Mo", "TCL1A+ naive B" = "Naive B") 
Idents(liver.immune) <- "cell.type"

liver.immune.markers.list <- list("CD4" = intersect(pbmc.markers$CD4.T.cells, rownames(liver.immune)), 
                                  "Treg" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "T regulatory cells", "official.gene.symbol"], 
                                                          rownames(liver.immune)),
                                  "mem B" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "B cells memory", "official.gene.symbol"], 
                                                     rownames(liver.immune)),
                                  "Naive B" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "B cells naive", "official.gene.symbol"], 
                                                      rownames(liver.immune)), 
                                  "Plasma" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Plasma cells", "official.gene.symbol"], 
                                                        rownames(liver.immune)), 
                                  "CD14+ Mo" = intersect(pbmc.markers$CD14.Monocytes, rownames(liver.immune)), 
                                  "NK" = intersect(pbmc.markers$NK.cells, rownames(liver.immune)), 
                                  "CD8" = intersect(pbmc.markers$CD8.T.cells, rownames(liver.immune)))
liver.immune <- subset(liver.immune, idents = names(liver.immune.markers.list))

print(length(unique(liver.immune$cell.type)))
print(dim(liver.immune))
sparsity(liver.immune@assays$RNA@counts)


# liver.immune <- readRDS("liver.immune.rds")

# DimPlot(liver.immune, group.by = "donor", label = T)
liver.immune <- liver.immune %>%
  NormalizeData(verbose = F) %>%
  FindVariableFeatures(verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = F) %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:10, verbose = F)
DimPlot(liver.immune, group.by = "cell.type", label = T)
saveRDS(liver.immune, "liver.immune.rds")
saveRDS(liver.immune.markers.list, "liver.immune.markers.rds")
write.table(do.call(qpcR:::cbind.na, liver.immune.markers.list), "liver.immune.markers.csv", sep = ",")

##
spleen.immune <- subset(liver.imm, idents = "donor 1 spleen")

spleen.immune.meta <- spleen.immune@meta.data
matx <- spleen.immune@assays$originalexp@counts
rownames(matx) <- liver.imm@assays$originalexp@meta.features$Symbol
spleen.immune <- CreateSeuratObject(matx, min.cells = 10)
spleen.immune@meta.data <- spleen.immune.meta

spleen.immune$cell.type <- as.character(spleen.immune$fine)
Idents(spleen.immune) <- "cell.type"
spleen.immune <- spleen.immune %>%
  NormalizeData(verbose = F) %>%
  FindVariableFeatures(verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = F) %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:10, verbose = F)
DimPlot(spleen.immune, group.by = "cell.type", label = T)
table(spleen.immune$cell.type)

spleen.immune$cell.type <- recode(spleen.immune$cell.type, "CCR7+ CD4" = "CD4", "CCR7+ CD8" = "CD8","CRCR3+ CD8" = "CD8","CRCR3+ CD4" = "CD4",
                                 "CX3CR1+ CD8" = "CD8", "CXCR6+ CD8" = "CD8","Naive B" = "Naive B", "cycling ASC" = "Plasma", "noncycline ASC" = "Plasma",
                                 "CD11c+ mem B" = "mem B","Classical mem B" = "mem B","CX3CR1+ NK" = "NK","CXCR6+ NK" = "NK", "Treg" = "Treg",
                                 "CD14+ Mo" = "CD14+ Mo", "CD16+ Mo" = "CD16+ Mo", "TCL1A+ naive B" = "Naive B", "CD1c+ DCs" = "DC") 
Idents(spleen.immune) <- "cell.type"

spleen.immune.markers.list <- list("CD4" = intersect(pbmc.markers$CD4.T.cells, rownames(spleen.immune)), 
                                   "Megakaryocyte" = intersect(pbmc.markers$Megakaryocytes, rownames(spleen.immune)), 
                                  "Treg" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "T regulatory cells", "official.gene.symbol"], 
                                                     rownames(spleen.immune)),
                                  "mem B" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "B cells memory", "official.gene.symbol"], 
                                                      rownames(spleen.immune)),
                                  "Naive B" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "B cells naive", "official.gene.symbol"], 
                                                        rownames(spleen.immune)), 
                                  "Plasma" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Plasma cells", "official.gene.symbol"], 
                                                       rownames(spleen.immune)), 
                                  "DC" = intersect(pbmc.markers$DC, rownames(spleen.immune)), 
                                  "CD14+ Mo" = intersect(pbmc.markers$CD14.Monocytes, rownames(spleen.immune)), 
                                  "CD16+ Mo" = intersect(pbmc.markers$CD16.Monocytes, rownames(spleen.immune)), 
                                  "NK" = intersect(pbmc.markers$NK.cells, rownames(spleen.immune)), 
                                  "CD8" = intersect(pbmc.markers$CD8.T.cells, rownames(spleen.immune)))
DC.markers <- pang.markers[pang.markers$species == "Hs" & pang.markers$cell.type == "Dendritic cells", ]
spleen.immune.markers.list[["DC"]] <- intersect(unique(c(spleen.immune.markers.list[["DC"]], DC.markers$official.gene.symbol)), rownames(spleen.immune))
Megakaryocyte.markers <- pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Megakaryocytes", ]
spleen.immune.markers.list[["Megakaryocyte"]] <- intersect(unique(c(spleen.immune.markers.list[["Megakaryocyte"]], Megakaryocyte.markers$official.gene.symbol)), rownames(spleen.immune))



spleen.immune <- subset(spleen.immune, idents = names(spleen.immune.markers.list))

print(length(unique(spleen.immune$cell.type)))
print(dim(spleen.immune))
sparsity(spleen.immune@assays$RNA@counts)


# spleen.immune <- readRDS("spleen.immune.rds")

# DimPlot(spleen.immune, group.by = "donor", label = T)
spleen.immune <- spleen.immune %>%
  NormalizeData(verbose = F) %>%
  FindVariableFeatures(verbose = F) %>%
  ScaleData(verbose = F) %>%
  RunPCA(verbose = F) %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:10, verbose = F)
DimPlot(spleen.immune, group.by = "cell.type", label = T)
saveRDS(spleen.immune, "spleen.immune.rds")
saveRDS(spleen.immune.markers.list, "spleen.immune.markers.rds")
write.table(do.call(qpcR:::cbind.na, spleen.immune.markers.list), "spleen.immune.markers.csv", sep = ",")


hcortex.sce <- NowakowskiCortexData()
hcortex <- CreateSeuratObject(hcortex.sce@assays@data$tpm, min.cells = 10)
hcortex$cell.type <- hcortex.sce$WGCNAcluster


hcortex$cell.type <- recode(hcortex$cell.type, "EN-PFC1" = "EN", "EN-PFC2" = "EN","EN-PFC3" = "EN","EN-V1-1" = "EN",
                            "EN-V1-2" = "EN","EN-V1-3" = "EN","IN-CTX-CGE1" = "IN", "IN-CTX-CGE2" = "IN", "IN-CTX-MGE1" = "IN", 
                            "IN-CTX-MGE2" = "IN", "IN-STR" = "IN", "nIN1" = "IN", "nIN2" = "IN", "nIN3" = "IN", "nIN4" = "IN", "nIN5" = "IN", 
                            "IPC-div1" = "IPC", "IPC-div2" = "IPC", "IPC-nEN1" = "IPC", 
                            "IPC-nEN2" = "IPC", "IPC-nEN3" = "IPC", "MGE-IPC1" = "IPC", "MGE-IPC2" = "IPC", "MGE-IPC3" = "IPC", 
                            "MGE-RG1" = "RG", "MGE-RG2" = "RG", "RG-div1" = "RG", "RG-div2" = "RG","RG-early" = "RG", "tRG" = "RG",
                            "vRG" = "RG","oRG" = "RG")
Idents(hcortex) <- "cell.type"
hcortex.markers.list <- list("Astrocyte" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Astrocytes", "official.gene.symbol"], 
                                                               rownames(hcortex)), 
                            "Choroid" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Choroid plexus cells", "official.gene.symbol"], 
                                                       rownames(hcortex)), 
                            "Microglia" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Microglia", "official.gene.symbol"], 
                                                    rownames(hcortex)), 
                            "EN" = intersect(cm2[cm2$cell_name == "Excitatory neuron" & cm2$tissue_class == "Brain", "Symbol"], 
                                                           rownames(hcortex)), 
                            "Endothelial" = intersect(cm2[cm2$cell_name == "Endothelial cell" & cm2$tissue_class == "Brain", "Symbol"], 
                                                      rownames(hcortex)),
                            "IN" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Interneurons", "official.gene.symbol"], 
                                                    rownames(hcortex)), 
                            "IPC" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Neural stem/precursor cells", "official.gene.symbol"], 
                                                    rownames(hcortex)), 
                            "RG" = intersect(pang.markers[pang.markers$species != "Mm" & pang.markers$cell.type == "Radial glia cells", "official.gene.symbol"], 
                                              rownames(hcortex))
                            )
hcortex <- subset(hcortex, idents = names(hcortex.markers.list))
print(length(unique(hcortex$cell.type)))
print(dim(hcortex))

sparsity(hcortex@assays$RNA@counts)


saveRDS(hcortex, "hcortex.rds")
saveRDS(hcortex.markers.list, "hcortex.markers.rds")
write.table(do.call(qpcR:::cbind.na, hcortex.markers.list), "hcortex.markers.csv", sep = ",")

#
lung.mtx <- read.delim("trachea_10x_log2TPM.txt", row.names = 1)
lung.meta <- read.delim("trachea_10x_metadata.txt")
rownames(lung.meta) <- lung.meta$NAME
lung.meta <- lung.meta[colnames(lung.mtx), ]
identical(colnames(lung.mtx), lung.meta$NAME)

lung <- CreateSeuratObject(lung.mtx, min.cells = 10)
lung$cell.type <- lung.meta$cluster
table(lung$cell.type)

# lung$cell.type <- recode(lung$cell.type, "Enterocyte.Mature.Distal" = "Enterocytes", "Enterocyte.Mature.Proximal" = "Enterocytes")
Idents(lung) <- "cell.type"
lung.markers.list <- list("Basal" = intersect(str_to_title(tolower(cm2[cm2$cell_name == "Basal cell" & cm2$tissue_class == "Lung", "Symbol"])), 
                                              rownames(lung)),  
                          "Ciliated" = intersect(cm2m[cm2m$cell_name == "Ciliated cell" & cm2m$tissue_class == "Lung", "Symbol"], 
                                                             rownames(lung)), 
                          "Club" = intersect(cm2m[cm2m$cell_name == "Clara cell" & cm2m$tissue_class == "Lung", "Symbol"], 
                                                 rownames(lung)), 
                          "Ionocyte" = intersect(str_to_title(tolower(cm2[cm2$cell_name == "Ionocyte cell" & cm2$tissue_class == "Lung", "Symbol"])), 
                                                 rownames(lung)), 
                          "Neuroendocrine" = intersect(str_to_title(tolower(cm2[cm2$cell_name == "Neuroendocrine cell" & cm2$tissue_class == "Lung", "Symbol"])), 
                                                 rownames(lung)), 
                          "Goblet" = intersect(str_to_title(tolower(pang.markers[pang.markers$species != "Hs" & pang.markers$cell.type == "Airway goblet cells", "official.gene.symbol"])), 
                                                    rownames(lung)), 
                          "Tuft" = intersect(str_to_title(tolower(cm2[cm2$cell_name == "Brush cell (Tuft cell)" & cm2$tissue_class == "Lung", "Symbol"])), 
                                                       rownames(lung)))
lung <- subset(lung, idents = names(lung.markers.list))

print(length(unique(lung$cell.type)))
print(dim(lung))

sparsity(lung@assays$RNA@counts)


saveRDS(lung, "lung.rds")
saveRDS(lung.markers.list, "lung.markers.rds")
write.table(do.call(qpcR:::cbind.na, lung.markers.list), "lung.markers.csv", sep = ",")



