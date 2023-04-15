# generate two files for SERGIO: one for master regulator expression for n.cell types; 
# one for regulator relationships
library(qpcR)
library(formattable)
library(pheatmap)
# set parameters
n.celltypes <- 4
n.cells.per.type <- 1000
n.genes <- 5000
n.master.regulators <- n.genes * 0.01
n.tight.strong.targets <- n.genes * 0.1
n.tight.medium.targets <- n.genes * 0.1
n.tight.weak.targets <- n.genes * 0.1
# n.multi.targets <- n.genes * 0.3
n.rest <- n.genes * 0.69

# n.tight.strong.targets
tight.strong.targets.info <- list()
for (i in 1:n.tight.strong.targets){
  cell <- i + n.master.regulators - 1
  # n.regs <- sample(1:3, 1, replace = FALSE)
  n.regs <- 1
  regs <- sample(0:(n.master.regulators-1), n.regs, replace = FALSE)
  strengths <- sample(c(-1, 1, 1, 1), size=length(regs), replace=TRUE) * runif(n = length(regs), min = 2, max = 5)
  hills <- rep(2.0, length(regs))
  target.info <- format(c(cell, n.regs, regs, strengths, hills), digits = 4, format = "f")
  # print(target.info)
  tight.strong.targets.info[[i]] <- target.info
}

tight.strong.targets.info.df <- do.call(qpcR:::rbind.na, tight.strong.targets.info)
#
# n.tight.medium.targets
tight.medium.targets.info <- list()
for (i in 1:n.tight.medium.targets){
  cell <- i + n.master.regulators + n.tight.strong.targets - 1
  # n.regs <- sample(1:3, 1, replace = FALSE)
  n.regs <- 1
  regs <- sample(0:(n.master.regulators-1), n.regs, replace = FALSE)
  strengths <- sample(c(-1,1,1,1), size=length(regs), replace=TRUE) * runif(n = length(regs), min = 1, max = 2)
  hills <- rep(2.0, length(regs))
  target.info <- format(c(cell, n.regs, regs, strengths, hills), digits = 4, format = "f")
  # print(target.info)
  tight.medium.targets.info[[i]] <- target.info
}

tight.medium.targets.info.df <- do.call(qpcR:::rbind.na, tight.medium.targets.info)

# n.tight.weak.targets
tight.weak.targets.info <- list()
for (i in 1:n.tight.weak.targets){
  cell <- i + n.master.regulators + n.tight.strong.targets + n.tight.medium.targets - 1
  # n.regs <- sample(1:3, 1, replace = FALSE)
  n.regs <- 1
  regs <- sample(0:(n.master.regulators-1), n.regs, replace = FALSE)
  strengths <- sample(c(-1,1,1,1), size=length(regs), replace=TRUE) * runif(n = length(regs), min = 0.5, max = 1)
  hills <- rep(2.0, length(regs))
  target.info <- format(c(cell, n.regs, regs, strengths, hills), digits = 4, format = "f")
  #  print(target.info)
  tight.weak.targets.info[[i]] <- target.info
}

tight.weak.targets.info.df <- do.call(qpcR:::rbind.na, tight.weak.targets.info)

# n.multi.targets
#multi.targets.info <- list()
#for (i in 1:n.multi.targets){
#  cell <- i + n.master.regulators + n.tight.strong.targets + n.tight.medium.targets + n.tight.weak.targets - 1
#  n.regs <- sample(2:4, 1, replace = FALSE)
#  # n.regs <- 1
#  regs <- sample(0:(n.master.regulators-1), n.regs, replace = FALSE)
#  strengths <- sample(c(-1,1,1,1), size=length(regs), replace=TRUE) * runif(n = length(regs), min = 0.5, max = 5)
#  hills <- rep(2.0, length(regs))
#  target.info <- format(c(cell, n.regs, regs, strengths, hills), digits = 4, format = "f")
#  #  print(target.info)
#  multi.targets.info[[i]] <- target.info
#}

#multi.targets.info.df <- do.call(qpcR:::rbind.na, multi.targets.info)
## n.rest
#rest.info <- list()
#for (i in 1:n.rest){
#  cell <- i + n.master.regulators + n.tight.strong.targets + n.tight.medium.targets + n.tight.weak.targets + n.multi.targets- 1
#  n.regs <- sample(5:10, 1, replace = FALSE)
#  regs <- sample(0:(n.master.regulators-1), n.regs, replace = FALSE)
#  strengths <- sample(c(-1, 1, 1, 1), size=length(regs), replace=TRUE) * runif(n = length(regs), min = 0.001, max = 1)
#  hills <- rep(2.0, length(regs))
#  target.info <- format(c(cell, n.regs, regs, strengths, hills), digits = 4, format = "f")
#  # print(target.info)
#  rest.info[[i]] <- target.info
#}

#rest.info.df <- do.call(qpcR:::rbind.na, rest.info)

#
reg.output <- qpcR:::rbind.na(tight.strong.targets.info.df, 
                              tight.medium.targets.info.df, 
                              tight.weak.targets.info.df)
#                              multi.targets.info.df,
#                              rest.info.df)
write.table(reg.output, "5000.gene.4ct.reg.dyn.3.txt", quote = F, col.names = F, row.names = F, sep = ",", na = "")

# master regulator expressions
M <- matrix(runif(n.celltypes * n.master.regulators, min = 0.1, max = 5), nrow = n.master.regulators)
m.reg <- as.data.frame(M)
pheatmap(m.reg)

rownames(m.reg) <- formattable(0:(n.master.regulators-1), digits = 4, format = "f")
write.table(m.reg, "5000.gene.4ct.m.reg.dyn.3.txt", quote = F, col.names = F, sep = ",")

# random gene expressions
total.items <- n.rest * n.celltypes * n.cells.per.type
random.numbers <- runif(total.items + 1, min = 0.1, max = 8)
Mr <- matrix(random.numbers, nrow = n.rest + 1)
m.r <- as.data.frame(Mr)
# pheatmap(m.reg)

rownames(m.r) <- formattable((n.master.regulators + 
                                n.tight.medium.targets + 
                                n.tight.strong.targets + 
                                n.tight.weak.targets):(n.genes-1), digits = 4, format = "f")
write.table(m.r, "5000.gene.4ct.random.dyn.3.txt", quote = F, col.names = F, sep = ",")



