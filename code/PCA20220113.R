setwd("/data2/usr/LiuWW/project_3/part5.global_analysis/PCA")

# 加载包
library(ggplot2)
library(ggrepel)
library(ggsci)

# 导入数据
data <- read.table(file="../genes_TPM_matrix.txt", header=T, sep="\t", row.names = 1, check.names = FALSE)
samp.inf <- read.table(file="../sample_info.txt", header=T, sep="\t", row.names = NULL, check.names = FALSE)


# 数据处理

##添加一列Famliy
samp.inf <- samp.inf %>% mutate(Family = origin_sample_name)

samp.inf$Family <- gsub(".*F", "F1", samp.inf$Family)
samp.inf$Family <- gsub(".*C.*", "Col-0", samp.inf$Family)
samp.inf$Family <- gsub("A.*|B.*", "Ecotype", samp.inf$Family)

##修改行名
row.names(samp.inf) <- samp.inf$Sample_name


############################
# 筛选基因1000个
############################
sd <- apply(data, 1, sd)
sd <- as.matrix(sd)
sd <- as.data.frame(sd[order(sd[,1], decreasing = TRUE),])
gene <- rownames(sd)[1:1000]


# PCA
## 计算主成分
df <- as.data.frame(t(data))
df_pca <- prcomp(df[,gene], scale=T) 
summary(df_pca)

## 主成分和分组
df_pcs <- data.frame(df_pca$x, Genotype = samp.inf[row.names(df_pca$x),]$Genotype, Batch = samp.inf[row.names(df_pca$x),]$Batch, Sample = samp.inf[row.names(df_pca$x),]$Sample_name, Organ = samp.inf[row.names(df_pca$x),]$Organ, Family = samp.inf[row.names(df_pca$x),]$Family)

## 查看主成分结果
head(df_pcs,3)

# 绘图
## 添加PC1、PC2解释率
xlab <- paste("PC1"," (",round((summary(df_pca))$importance[2,1]*100,2),"%)",sep="")
ylab <- paste("PC2"," (",round((summary(df_pca))$importance[2,2]*100,2),"%)",sep="")

mytheme <- theme_bw() + theme(plot.background=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line= element_line(colour ="black"), plot.title = element_text(size=9, face = "bold", hjust = 0.5))

p <- ggplot(df_pcs, aes(x = PC1, y = PC2, color = Batch, shape = Organ, fill = Family, label = Genotype)) +
  geom_point(show.legend = TRUE, size = 3, stroke = 1.5) + 
  #geom_text_repel(size=2, col="black", segment.color="black", segment.size = .1, segment.alpha = I(1/3)) +
  xlab(xlab) + ylab(ylab) + 
  ggtitle("TPM_gene_sd_top1000") +
  coord_equal(ratio=1) + 
  mytheme + scale_color_manual(values = c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2", "#B09C85B2", "#FEB24C")) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"), breaks = c("Col-0", "Ecotype", "F1")) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

pdf(file="20220113_genes.TMM.EXPR.matrix_PCA_TPM_sd_top1000.pdf", height=10, width=8)
p
dev.off()
