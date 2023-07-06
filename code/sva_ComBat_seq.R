setwd("/data2/usr/LiuWW/project_3/de_batch_effect")

# package
library(sva)
library(tidyverse)

# data
load("/data2/usr/LiuWW/project_3/part5.global_analysis/genes_counts_matrix.RData")
load("/data2/usr/LiuWW/project_3/part5.global_analysis/sample_info.RData")

# batch effect adjustment
## 准备counts矩阵
counts_matrix <- column_to_rownames(.data = counts_new, var = "Gene_id")
## 准备batch和group信息
info <- sample_info %>% 
  filter(Sample_name %in% colnames(counts_matrix)) %>%
  mutate(Group = paste(Organ, Genotype, sep = "_"))

info <- column_to_rownames(.data = info, var = "Sample_name")

info <- info[colnames(counts_matrix),]

#########################
## ShootApex（SA）和TrueLeaf（TL）分开进行批次效应矫正
#########################


### SA和TL分别准备counts矩阵和样品info
#### counts矩阵
counts_matrix_SA <- counts_matrix[,c(1:203)]
counts_matrix_TL <- counts_matrix[,c(204:406)]
#### 样品info
info_SA <- info[c(1:203),]
info_TL <- info[c(204:406),]


### ComBat_seq for SA
batch <- as.factor(info_SA$Batch)

cov1 <- as.factor(info_SA$Organ)
cov2 <- as.factor(info_SA$Group)

covar_mat <- cbind(cov1, cov2)

adjusted_counts_SA <- sva::ComBat_seq(as.matrix(counts_matrix_SA), batch = batch, group = NULL, covar_mod = covar_mat)

### log for SA
# Found 11 batches
# Using null model in ComBat-seq.
# Adjusting for 1 covariate(s) or covariate level(s)
# Estimating dispersions
# Fitting the GLM model
# Shrinkage off - using GLM estimates for parameters
# Adjusting the data

# output for SA
## 不存在则新建目录adjusted_counts
if (! file.exists("./adjusted_counts")){
  dir.create("./adjusted_counts")
} 
## 结果分单个样本保存至adjusted_counts
for (i in colnames(adjusted_counts_SA)) {
  counts <- as.matrix(adjusted_counts_SA[, i])
  colnames(counts) <- i
  save(counts, file = paste0("./adjusted_counts/", i, ".RData"))
}


### ComBat_seq for TL
batch <- as.factor(info_TL$Batch)

cov1 <- as.factor(info_TL$Organ)
cov2 <- as.factor(info_TL$Group)

covar_mat <- cbind(cov1, cov2)

adjusted_counts_TL <- sva::ComBat_seq(as.matrix(counts_matrix_TL), batch = batch, group = NULL, covar_mod = covar_mat)

### log for TL
# Found 11 batches
# Using null model in ComBat-seq.
# Adjusting for 1 covariate(s) or covariate level(s)
# Estimating dispersions
# Fitting the GLM model
# Shrinkage off - using GLM estimates for parameters
# Adjusting the data

# output for TL
## 不存在则新建目录adjusted_counts
if (! file.exists("./adjusted_counts")){
  dir.create("./adjusted_counts")
} 
## 结果分单个样本保存至adjusted_counts
for (i in colnames(adjusted_counts_TL)) {
  counts <- as.matrix(adjusted_counts_TL[, i])
  colnames(counts) <- i
  save(counts, file = paste0("./adjusted_counts/", i, ".RData"))
}