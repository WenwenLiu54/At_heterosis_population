setwd("/home/ug0096/expr_DE")
library(edgeR)

# phenotype
load("/home/ug0096/barplot/phenotype.RData")

# 使用EdgeR设置bcv参数做没有样本重复的差异分析
# 获得表达矩阵
expr_count <- read.table(file = "./genes.sva.counts.matrix", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

C <- grep('Col', colnames(expr_count), value = T)
F1 <- grep('F1', colnames(expr_count), value = T)
need <- c(C, F1)
expr_count <- expr_count[, need]


# 获得分组信息
group <- need


# 没有重复的差异分析
# 创建DGEList
y <- DGEList(counts=expr_count, group=group)
# 差异分析, 注意pair参数，设置为 F1 vs Col
bcv = 0.1 #设置bcv为0.1

##---------------------- SA ----------------------------------
dir.create("./DE_SA")
for (f in F1[1:96]) {
  c <- paste0(str_sub(f, 1, 19), "Col-0")
  
  et <- exactTest(y, dispersion=bcv^2, pair = c(c, f))
  
  DEG_edgeR=as.data.frame(topTags(et, n = nrow(expr_count)))
  
  # 阈值|logFC|>=1，pvalue<0.05
  DEG_edgeR$change=ifelse(DEG_edgeR$logFC>=1&DEG_edgeR$FDR<0.05,"up",
                        ifelse(DEG_edgeR$logFC<=-1&DEG_edgeR$FDR<0.05,"down","not"))
  
  DEG_edgeR <- rownames_to_column(DEG_edgeR, var = "geneid")
  
  save(DEG_edgeR, file = paste0("./DE_SA/", f, ".RData"))
}

### SA MERGE
load(paste0("./DE_SA/", F1[1], ".RData"))
SA <- DEG_edgeR[, c("geneid", "change")]
colnames(SA) <- c("geneid", F1[1])
for (f in F1[2:96]) {
  load(paste0("./DE_SA/", f, ".RData"))
  tmp <- DEG_edgeR[, c("geneid", "change")]
  colnames(tmp) <- c("geneid", f)
  
  SA <- left_join(SA, tmp, by = "geneid")
}

save(SA, file = "SA_DE.RData")
load("./SA_DE.RData")

##---------------------- TL ----------------------------------
dir.create("./DE_TL")
for (f in F1[97:192]) {
  c <- paste0(str_sub(f, 1, 18), "Col-0")
  
  et <- exactTest(y, dispersion=bcv^2, pair = c(c, f))
  
  DEG_edgeR=as.data.frame(topTags(et, n = nrow(expr_count)))
  
  # 阈值|logFC|>=1，pvalue<0.05
  DEG_edgeR$change=ifelse(DEG_edgeR$logFC>=1&DEG_edgeR$FDR<0.05,"up",
                          ifelse(DEG_edgeR$logFC<=-1&DEG_edgeR$FDR<0.05,"down","not"))
  
  DEG_edgeR <- rownames_to_column(DEG_edgeR, var = "geneid")
  
  save(DEG_edgeR, file = paste0("./DE_TL/", f, ".RData"))
}

### TL MERGE
load(paste0("./DE_TL/", F1[97], ".RData"))
TL <- DEG_edgeR[, c("geneid", "change")]
colnames(TL) <- c("geneid", F1[97])
for (f in F1[98:192]) {
  load(paste0("./DE_TL/", f, ".RData"))
  tmp <- DEG_edgeR[, c("geneid", "change")]
  colnames(tmp) <- c("geneid", f)
  
  TL <- left_join(TL, tmp, by = "geneid")
}

save(TL, file = "TL_DE.RData")
load("./TL_DE.RData")


##----------------------------TWAS GENE PLOT-----------------------------
gene <- read.table(file = "/home/ug0096/part9.twas/paper_table_make/geneid_cor.txt", header = T, sep = "\t", check.names = F)

# PICK GENE
SA_POS <- gene %>% filter(TWAS编号 == "TWAS 1" & 相关系数 > 0) %>% select(基因ID)
SA_POS <- as.vector(SA_POS$基因ID)
SA_NEG <- gene %>% filter(TWAS编号 == "TWAS 1" & 相关系数 < 0) %>% select(基因ID)
SA_NEG <- as.vector(SA_NEG$基因ID) 

TL_POS <- gene %>% filter(TWAS编号 == "TWAS 7" & 相关系数 > 0) %>% select(基因ID)
TL_POS <- as.vector(TL_POS$基因ID)
TL_NEG <- gene %>% filter(TWAS编号 == "TWAS 7" & 相关系数 < 0) %>% select(基因ID)
TL_NEG <- as.vector(TL_NEG$基因ID)


# PICK DE OF GENE
colnames(SA) <- c("geneid", str_sub(colnames(SA)[2:97], 20, -4))
SA_POS_PLOT <- SA %>% filter(geneid %in% SA_POS) %>% select(str_sub(plot_data$Genotype, 4, -1))
SA_NEG_PLOT <- SA %>% filter(geneid %in% SA_NEG) %>% select(str_sub(plot_data$Genotype, 4, -1))

colnames(TL) <- c("geneid", str_sub(colnames(TL)[2:97], 19, -4)) 
TL_POS_PLOT <- TL %>% filter(geneid %in% TL_POS) %>% select(str_sub(plot_data$Genotype, 4, -1))
TL_NEG_PLOT <- TL %>% filter(geneid %in% TL_NEG) %>% select(str_sub(plot_data$Genotype, 4, -1))


# ADD FLAG FOR HEATMAP
SA_POS_PLOT$flag <- "SA_POS_PLOT"
SA_NEG_PLOT$flag <- "SA_NEG_PLOT"

TL_POS_PLOT$flag <- "TL_POS_PLOT"
TL_NEG_PLOT$flag <- "TL_NEG_PLOT"

# GET PLOT MATRIX
df <- rbind(SA_POS_PLOT, SA_NEG_PLOT, TL_POS_PLOT, TL_NEG_PLOT)

# 热图风格
ht_global_opt(
  legend_border = "black",
  heatmap_border = TRUE,
  annotation_border = TRUE)


h <- Heatmap(df[,1:96],
             name = 'F1vsC DE',
             col = c("up"="yellow", "down"="blue", "not"="gray"),
             show_row_names = F,
             cluster_columns = F,
             cluster_rows = F,
             column_title = "hybrids",
) + 
  Heatmap(df$flag,
          name = 'flag', 
          col = c("SA_POS_PLOT" = "#B388EB", 
                  "TL_POS_PLOT" = "#F7AEF8",
                  "SA_NEG_PLOT" = "#6699CC",
                  "TL_NEG_PLOT" = "#81B29A")
  )


pdf(file = "DE_heatmap.pdf", width = 9, height = 7)
draw(h, heatmap_legend_side = 'bottom')
dev.off()


##----------------------------plot barplot-----------------------------------
bar <- list()
for (i in c("SA_POS_PLOT", "SA_NEG_PLOT", "TL_POS_PLOT", "TL_NEG_PLOT")) {
  high <- get(i)[, c(65:96)]
  low <- get(i)[, c(1:32)]
  high_up <- apply(high, 2, function(x){length(which(x == "up"))})
  high_down <- apply(high, 2, function(x){length(which(x == "down"))})
  low_up <- apply(low, 2, function(x){length(which(x == "up"))})
  low_down <- apply(low, 2, function(x){length(which(x == "down"))})
  
  ##plot
  data <- rbind(data.frame(value = high_up, flag="high_up"), data.frame(value = low_up, flag="low_up"), data.frame(value = high_down, flag="high_down"), data.frame(value = low_down, flag="low_down"))
  rownames(data) <- NULL
  
  dm <- data %>% group_by(flag) %>%
    summarise(mean_count = mean(value),
              mean_count = mean(value),
              sd_count = sd(value),
              sd_count = sd(value))
  dm$class <- str_sub(dm$flag, -2, -1)
  dm$flag = factor(dm$flag, levels = c("high_up", "low_up", "high_down", "low_down"))
  
  up <- t.test(high_up, low_up, paired = FALSE)
  down <- t.test(high_down, low_down, paired = FALSE)
  
  bar[[i]] <- ggplot(dm, aes(x = flag, y = mean_count, fill = class)) + 
    geom_col(width = .7) + 
    geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count), position = position_dodge(0.9), width = .2) +
    labs(title = i,
         y = "count") +
    theme_classic() +
    scale_fill_npg() +
    geom_signif(annotations = c(signif(up$p.value, 3), signif(down$p.value, 3)), y_position = c(15, 15), xmin = c(1,3), xmax = c(2,4))
}

# 组图
p <- bar[["SA_POS_PLOT"]] + bar[["SA_NEG_PLOT"]] + bar[["TL_POS_PLOT"]] + bar[["TL_NEG_PLOT"]] + plot_layout(guides = "collect") & theme(legend.position = 'top')

ggsave(filename = "DE_barplot.pdf", plot = p, width = 20, height = 25, units = "cm")
