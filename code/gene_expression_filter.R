setwd("/data2/usr/LiuWW/project_3/part6.global_analysis_debatch/gene_expression_filter")

# 加载包
library(tidyverse)

# 导入数据
data <- read.table(file="../../de_batch_effect/adjusted_merge_to_matrix/genes.sva.TMM.EXPR.matrix", header=T, sep="\t", row.names = 1, check.names = FALSE)

# 数据预处理
colnames(data) <- str_sub(colnames(data), 1, -9)
data <- round(data, digits = 3)

# 基因表达过滤
# discard genes with expression levels less than 1 (TPM < 1) in more than 98% of samples.(406*0.98=397.88)
nokeep <- rowSums(data < 1) >= 398

data_after_filter <- data[!nokeep,]

genes_high <- rownames(data_after_filter) 

genes_high <- sort(genes_high)

# 保存过滤后的表达矩阵为RData
TPM_genes_high <- data_after_filter[genes_high,]
save(TPM_genes_high, file = "TPM_genes_high.RData")

# 保存排序后的genes_high为RData
save(genes_high, file = "genes_high.RData")