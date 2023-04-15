# quantile normalization for dyn data
library(preprocessCore)
library(Seurat)
library(dplyr)

design <- "5000.gene.4ct"
n.cells <- 4000
n.cells.per.type <- 1000
n.genes <- 5000
version.n <- 1


output_dyn_data <- function(design, version.n){
  gen <- read.csv(paste0(design, ".dyn.v", version.n, ".S.clean.csv"), header = F)
  rand <- read.csv(paste0(design, ".random.dyn.", version.n, ".txt"), header = F, row.names = 1)
  
  a <- normalize.quantiles.use.target(as.matrix(gen),target = as.vector(as.matrix(rand)),copy=TRUE,subset=NULL)
  rownames(a) <- rownames(a) - 1
  colnames(a) <- colnames(rand)
  
  gex <- rbind(a, rand)
  write.table(gex, paste0(design, ".v", version.n, ".dyn.clean.csv"), quote = F, row.names = F, col.names = F, sep = ",")
}

output_dyn_data(design = design, version.n = 3)
output_dyn_data(design = "5000.gene.6ct", version.n = 1)
output_dyn_data(design = "5000.gene.6ct", version.n = 2)
output_dyn_data(design = "5000.gene.6ct", version.n = 3)


# ``````````````````
# below are codes for code testing
gen <- read.csv(paste0(design, ".dyn.v", version.n, ".S.clean.csv"), header = F)
rand <- read.csv(paste0(design, ".random.dyn.", version.n, ".txt"), header = F, row.names = 1)

a <- normalize.quantiles.use.target(as.matrix(gen),target = as.vector(as.matrix(rand)),copy=TRUE,subset=NULL)
rownames(a) <- rownames(a) - 1
colnames(a) <- colnames(rand)

gex <- rbind(a, rand)
# colnames(gex) <- paste0("V", 1:n.cells)
write.table(gex, paste0(design, ".v", version.n, ".dyn.clean.csv"), quote = F, row.names = F, col.names = F, sep = ",")

# 
rownames(gex) <- paste0("gene", rownames(gex))
colnames(gex) <- paste0("cell_", colnames(gex))
sr <- CreateSeuratObject(gex)
sr$cell.type <- c(rep("cell.type.1", n.cells.per.type), rep("cell.type.2", n.cells.per.type), 
                  rep("cell.type.3", n.cells.per.type), rep("cell.type.4", n.cells.per.type))

sr <- sr %>%
  NormalizeData() %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:10)

DimPlot(sr, group.by = "cell.type")
