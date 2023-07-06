setwd("/data2/usr/LiuWW/project_3/part5.global_analysis/sample_cor")

# 清除当前环境中的变量
rm(list=ls())

# 加载包
library(ComplexHeatmap)
require(circlize)
library(RColorBrewer)
library(ggplot2)

# 导入数据
data <- read.table(file="../genes_TPM_matrix.txt", header=T, sep="\t", row.names = 1, check.names = FALSE)

load("/data2/usr/LiuWW/project_3/part5.global_analysis/sample_info.RData")


# 数据预处理
###################
# data
###################

###################
# sample_info
###################
##添加一列Famliy
samp.inf <- sample_info %>% mutate(Family = origin_sample_name)

samp.inf$Family <- gsub(".*F", "F1", samp.inf$Family)
samp.inf$Family <- gsub(".*C.*", "Col-0", samp.inf$Family)
samp.inf$Family <- gsub("A.*|B.*", "Ecotype", samp.inf$Family)

##修改行名
row.names(samp.inf) <- samp.inf$Sample_name

##取406行
samp.inf <- samp.inf[colnames(data), ]

##取所需列
samp.inf <- samp.inf[,c("Batch", "Organ", "Family")]

##变量由字符型变因子型
samp.inf$Batch <- as.factor(samp.inf$Batch)
samp.inf$Organ <- as.factor(samp.inf$Organ)
samp.inf$Family <- as.factor(samp.inf$Family)


#############################
# SA和TL分开
#############################

data_SA <- data[, c(1:203)]
data_TL <- data[, c(204:406)]

samp.inf_SA <- samp.inf[c(1:203),]
samp.inf_TL <- samp.inf[c(204:406),]

############################
## plot for SA
############################

# 样本相关性矩阵计算
data1 <- data.matrix(data_SA)
cor <- round(cor(data1, method="pearson"), 4)
write.table(cor, file="SA_genes.TMM.EXPR.matrix_sample_cor_matrix_gene_32833.txt", sep="\t")

# 自定义颜色，colorRamp2函数来自于circlize包
mycol <- colorRamp2(c(min(cor), mean(cor), max(cor)), c("#6BAED6", "black", "#FEC44F"))

# HeatmapAnnotation函数对列进行注释，rowAnnotation函数对行进行注释
top_anno <- HeatmapAnnotation(df = samp.inf_SA, col = list(Batch = c("batch_01" = "#E64B35B2", "batch_02" = "#4DBBD5B2", "batch_03" = "#00A087B2", "batch_04" = "#3C5488B2", "batch_05" = "#F39B7FB2", "batch_06" = "#8491B4B2", "batch_07" = "#91D1C2B2", "batch_08" = "#DC0000B2", "batch_09" = "#7E6148B2", "batch_10" = "#B09C85B2", "batch_11" = "#FEB24C"), Organ = c("ShootApex" = "Gold", "TrueLeaf" = "LightSkyBlue"), Family = c("Col-0" = "#00A087B2", "F1" = "#3C5488B2", "Ecotype" = "#E64B35B2")))

left_anno <- rowAnnotation(df = samp.inf_SA, col = list(Batch = c("batch_01" = "#E64B35B2", "batch_02" = "#4DBBD5B2", "batch_03" = "#00A087B2", "batch_04" = "#3C5488B2", "batch_05" = "#F39B7FB2", "batch_06" = "#8491B4B2", "batch_07" = "#91D1C2B2", "batch_08" = "#DC0000B2", "batch_09" = "#7E6148B2", "batch_10" = "#B09C85B2", "batch_11" = "#FEB24C"), Organ = c("ShootApex" = "Gold", "TrueLeaf" = "LightSkyBlue"), Family = c("Col-0" = "#00A087B2", "F1" = "#3C5488B2", "Ecotype" = "#E64B35B2")))


##################做cluster##############################
h1 <- Heatmap(cor, 
        name = "pearson", # name参数设定图例标题
        col = mycol,
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = FALSE,
        # heatmap_width 和 heatmap_height 控制整个热图的宽度/高度，包括所有热图组件（不包括图例），而 width 和 height 仅控制 heamtap 主体的宽度/高度.
        heatmap_height = unit(1, "mm")*nrow(cor),
        heatmap_width = unit(1, "mm")*ncol(cor),
        height = 0.1,
        width = 0.1,
        top_annotation = top_anno,
        left_annotation = left_anno)

pdf(file="SA_genes.TMM.EXPR.matrix_sample_cor_gene_32833_cluster.pdf", height=20, width=20)
h1
dev.off()


##################不做cluster##############################
h2 <- Heatmap(cor, 
              name = "pearson", # name参数设定图例标题
              col = mycol,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = FALSE, show_column_names = FALSE,
              # heatmap_width 和 heatmap_height 控制整个热图的宽度/高度，包括所有热图组件（不包括图例），而 width 和 height 仅控制 heamtap 主体的宽度/高度.
              heatmap_height = unit(1, "mm")*nrow(cor),
              heatmap_width = unit(1, "mm")*ncol(cor),
              height = 0.1,
              width = 0.1,
              top_annotation = top_anno,
              left_annotation = left_anno)

pdf(file="SA_genes.TMM.EXPR.matrix_sample_cor_gene_32833_NOcluster.pdf", height=20, width=20)
h2
dev.off()


############################
## plot for TL
############################

# 样本相关性矩阵计算
data1 <- data.matrix(data_TL)
cor <- round(cor(data1, method="pearson"), 4)
write.table(cor, file="TL_genes.TMM.EXPR.matrix_sample_cor_matrix_gene_32833.txt", sep="\t")

# 自定义颜色，colorRamp2函数来自于circlize包
mycol <- colorRamp2(c(min(cor), mean(cor), max(cor)), c("#6BAED6", "black", "#FEC44F"))

# HeatmapAnnotation函数对列进行注释，rowAnnotation函数对行进行注释
top_anno <- HeatmapAnnotation(df = samp.inf_TL, col = list(Batch = c("batch_01" = "#E64B35B2", "batch_02" = "#4DBBD5B2", "batch_03" = "#00A087B2", "batch_04" = "#3C5488B2", "batch_05" = "#F39B7FB2", "batch_06" = "#8491B4B2", "batch_07" = "#91D1C2B2", "batch_08" = "#DC0000B2", "batch_09" = "#7E6148B2", "batch_10" = "#B09C85B2", "batch_11" = "#FEB24C"), Organ = c("ShootApex" = "Gold", "TrueLeaf" = "LightSkyBlue"), Family = c("Col-0" = "#00A087B2", "F1" = "#3C5488B2", "Ecotype" = "#E64B35B2")))

left_anno <- rowAnnotation(df = samp.inf_TL, col = list(Batch = c("batch_01" = "#E64B35B2", "batch_02" = "#4DBBD5B2", "batch_03" = "#00A087B2", "batch_04" = "#3C5488B2", "batch_05" = "#F39B7FB2", "batch_06" = "#8491B4B2", "batch_07" = "#91D1C2B2", "batch_08" = "#DC0000B2", "batch_09" = "#7E6148B2", "batch_10" = "#B09C85B2", "batch_11" = "#FEB24C"), Organ = c("ShootApex" = "Gold", "TrueLeaf" = "LightSkyBlue"), Family = c("Col-0" = "#00A087B2", "F1" = "#3C5488B2", "Ecotype" = "#E64B35B2")))


##################做cluster##############################
h1 <- Heatmap(cor, 
              name = "pearson", # name参数设定图例标题
              col = mycol,
              cluster_rows = TRUE, cluster_columns = TRUE,
              show_row_names = FALSE, show_column_names = FALSE,
              # heatmap_width 和 heatmap_height 控制整个热图的宽度/高度，包括所有热图组件（不包括图例），而 width 和 height 仅控制 heamtap 主体的宽度/高度.
              heatmap_height = unit(1, "mm")*nrow(cor),
              heatmap_width = unit(1, "mm")*ncol(cor),
              height = 0.1,
              width = 0.1,
              top_annotation = top_anno,
              left_annotation = left_anno)

pdf(file="TL_genes.TMM.EXPR.matrix_sample_cor_gene_32833_cluster.pdf", height=20, width=20)
h1
dev.off()


##################不做cluster##############################
h2 <- Heatmap(cor, 
              name = "pearson", # name参数设定图例标题
              col = mycol,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = FALSE, show_column_names = FALSE,
              # heatmap_width 和 heatmap_height 控制整个热图的宽度/高度，包括所有热图组件（不包括图例），而 width 和 height 仅控制 heamtap 主体的宽度/高度.
              heatmap_height = unit(1, "mm")*nrow(cor),
              heatmap_width = unit(1, "mm")*ncol(cor),
              height = 0.1,
              width = 0.1,
              top_annotation = top_anno,
              left_annotation = left_anno)

pdf(file="TL_genes.TMM.EXPR.matrix_sample_cor_gene_32833_NOcluster.pdf", height=20, width=20)
h2
dev.off()
