# cell chat
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(dplyr)
library(ComplexHeatmap)
# construct data
tnbc1.tu.obj <- readRDS("tnbc1.tu.obj.rds")
tnbc1.normal <- readRDS("tnbc1.normal.rds")

ar.p <- subset(tnbc1.tu.obj, subset = HALLMARK_G2M_CHECKPOINT == "positive")
ar.p$cell.type <- "tumor"
Idents(ar.p) <- "cell.type"
ar.n <- subset(tnbc1.tu.obj, subset = HALLMARK_G2M_CHECKPOINT == "negative")
ar.n$cell.type <- "tumor"
Idents(ar.n) <- "cell.type"
ar.n <- subset(ar.n, downsample = ncol(ar.p))

tnbc1.normal$cell.type <- Idents(tnbc1.normal)
# immune <- subset(tnbc1.normal, idents = "Immune")
# immune$cell.type <- "immune"

ar.p.cc <- merge(ar.p, tnbc1.normal)
ar.n.cc <- merge(ar.n, tnbc1.normal)

#
cellchat1 <- createCellChat(object = as.matrix(ar.p.cc@assays$RNA@data), meta = ar.p.cc@meta.data, group.by = "cell.type")
cellchat2 <- createCellChat(object = as.matrix(ar.n.cc@assays$RNA@data), meta = ar.n.cc@meta.data, group.by = "cell.type")

runcc <- function(ccobj){
  CellChatDB <- CellChatDB.human 
  ccobj@DB <- CellChatDB
  
  # process
  ccobj <- subsetData(ccobj)
  ccobj <- identifyOverExpressedGenes(ccobj)
  ccobj <- identifyOverExpressedInteractions(ccobj)
  # project gene expression data onto PPI network (optional)
  ccobj <- projectData(ccobj, PPI.human)
  # compute interaction
  ccobj <- computeCommunProb(ccobj)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  ccobj <- filterCommunication(ccobj, min.cells = 5)
  
  # pathway level
  ccobj <- computeCommunProbPathway(ccobj)
  
  ccobj <- aggregateNet(ccobj)
  ccobj <- netAnalysis_computeCentrality(ccobj, slot.name = "netP")
  
}

cellchat1 <- runcc(ccobj = cellchat1)
cellchat2 <- runcc(ccobj = cellchat2)

object.list <- list(neg = cellchat2, pos = cellchat1)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)
saveRDS(cellchat, "cellchat.rds")
