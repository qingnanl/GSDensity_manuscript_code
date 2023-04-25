library(Seurat)
library(matrixStats)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# n.master.regulators.r <- 0.005
# n.tight.strong.targets.r <- 0.1
# n.tight.medium.targets.r <- 0.1
# n.tight.weak.targets.r <- 0.1
# n.multi.targets.r <- 0.3
# n.rest.r <- 0.395
#

# design <- "5000.gene.3ct"
# n.genes <- 5000
# n.cell.types <- 3
# version.n <- 1

preprocess_regulation <- function(design, version.n, n.tight.strong.targets.r, 
                                  n.tight.medium.targets.r, n.tight.weak.targets.r, 
                                  n.multi.targets.r, n.rest.r, n.cell.types, n.genes){
  summary_list <- list()
  master.reg.file <- ifelse(file.exists(paste0(design, ".m.reg.", version.n, ".txt")), 
                            paste0(design, ".m.reg.", version.n, ".txt"), 
                            paste0(design, ".m.reg.dyn.", version.n, ".txt"))
  regulation.file <- ifelse(file.exists(paste0(design, ".reg.", version.n, ".txt")), 
                            paste0(design, ".reg.", version.n, ".txt"), 
                            paste0(design, ".reg.dyn.", version.n, ".txt"))

  master.reg <- read.delim(master.reg.file, sep = ",", row.names = 1, header = F)
  colnames(master.reg) <- paste0("cell", 1:n.cell.types)
  rownames(master.reg) <- paste0("gene", rownames(master.reg))
  
  regulation <- read.delim(regulation.file, sep = ",", row.names = 1, header = F)
  
  rownames(regulation) <- paste0("gene", rownames(regulation))
  
  tight.regulation <- regulation[1:(n.genes * (n.tight.strong.targets.r + n.tight.medium.targets.r + n.tight.weak.targets.r)), ]
  colnames(tight.regulation) <- c("n.regulators", "regulator", "strength", "hill.coeff")
  tight.regulation[, "regulator"] <- paste0("gene", tight.regulation[, "regulator"])
  tight.regulation <- as.data.frame(tight.regulation)
  # find cell specific regulators
  regulator.list <- list()
  for (cell.type in paste0("cell", 1:n.cell.types)){
    master.reg.exp1 <- master.reg[, cell.type]
    master.reg.exp2 <- rowMaxs(as.matrix(master.reg[, !names(master.reg) %in% cell.type]))
    specificity.ratio <- master.reg.exp1/master.reg.exp2
    names(specificity.ratio) <- rownames(master.reg)
    regulator.list[[cell.type]] <- specificity.ratio
  }

  target.list <- list()
  
  for (cell.type in paste0("cell", 1:n.cell.types)){
    master.reg.gene <- names(sort(regulator.list[[cell.type]], decreasing = T)[1:3]) # top 3 specific master regulators
    target.list[[paste0(cell.type, ".strong")]] <- rownames(tight.regulation[tight.regulation$regulator %in% master.reg.gene & 
                                                                               tight.regulation$strength >= 2, ])
    target.list[[paste0(cell.type, ".medium")]] <- rownames(tight.regulation[tight.regulation$regulator %in% master.reg.gene & 
                                                                               tight.regulation$strength <= 2 &
                                                                               tight.regulation$strength >= 1, ])
    target.list[[paste0(cell.type, ".weak")]] <- rownames(tight.regulation[tight.regulation$regulator %in% master.reg.gene & 
                                                                             tight.regulation$strength <= 1 &
                                                                             tight.regulation$strength >= 0.5, ])
  }
  
  summary_list[["target.list"]] <- target.list
  summary_list[["regulator.list"]] <- regulator.list
  summary_list[["regulation"]] <- regulation
  summary_list[["master.reg"]] <- master.reg
  return(summary_list)
}

# sum <- preprocess_regulation(design = "5000.gene.3ct", version.n = 1, n.tight.strong.targets.r = 0.1, 
#                             n.tight.medium.targets.r = 0.1, n.tight.weak.targets.r = 0.1, 
#                             n.multi.targets.r = 0.3, n.rest.r = 0.395, n.cell.types = 3, n.genes = 5000)



# heatmap for master regulator expression
heatmap.mr <- function(summary_list, title){
  hm <- pheatmap(sum$master.reg, main = title)
  return(hm)
}
# heatmap.mr(summary_list = sum, title = "5000 genes 3 ct version 1")
#

# 1d ordering plot for specificity ratio 
spec.ratio.plot <- function(summary_list, cell.type, title){
  plot.df <- data.frame(specificity_score = summary_list$regulator.list[[cell.type]], 
                        rank = rank(summary_list$regulator.list[[cell.type]]), 
                        gene_id = names(summary_list$regulator.list[[cell.type]]))
  p <- ggplot(data = plot.df, aes(x = rank, y = specificity_score, label = gene_id)) +        
    geom_text_repel() + 
    geom_point(alpha = 0.5) + 
    theme_classic() + ggtitle(title) +
    xlab("Rank of master regulators") + ylab("Specificity score") + 
    theme(
      plot.title = element_text(color="black", size=14, face="bold", hjust = 0.5),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold")
    )
  return(p)
}

# spec.ratio.plot(summary_list = sum, cell.type = "cell1", title = "5000 genes 3 ct version 1")

# sort(sum$regulator.list[["cell1"]])

# Create and preprocess Seurat objects

gen.cell.vector <- function(n.cell.types, n.cells.per.type){
  vec <- c()
  for (i in 1:n.cell.types){
    vec.i <- rep(paste0("cell", i), n.cells.per.type)
    vec <- c(vec, vec.i)
  }
  return(vec)
}

preprocess_sr <- function(design, version.n, drop.out, n.cell.types, n.cells.per.type){
  gex <- read.csv(paste0(design, ".v", version.n, ".", drop.out, ".csv"), header = F)
  n.cells <- ncol(gex)
  n.cells.per.type <- n.cells/n.cell.types
  rownames(gex) <- paste0("gene", as.numeric(rownames(gex)) - 1)
  colnames(gex) <- paste0("cell_", colnames(gex))
  sr <- CreateSeuratObject(gex)
  sr$cell.type <- gen.cell.vector(n.cell.types = n.cell.types, n.cells.per.type = n.cells.per.type)
  sr <- sr %>%
    NormalizeData() %>% 
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors() %>%
    RunUMAP(dims = 1:10)
  return(sr)
}

# sr <- preprocess_sr(design = "5000.gene.3ct", version.n = 1, drop.out = "dp1", 
#                    n.cell.types = 3, n.cells.per.type = 1500)

#DimPlot(sr, group.by = "cell.type")

#
#c(rep("cell.type.1", n.cells.per.type), rep("cell.type.2", n.cells.per.type), rep("cell.type.3", n.cells.per.type))
#drop.out <- "dp3"