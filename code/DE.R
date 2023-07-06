setwd("/home/ug0096/DE_heatmap")
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

load("/home/ug0096/LA.RData")
load("/home/ug0096/CA.RData")
load("/home/ug0096/CN.RData")
load("/home/ug0096/96pop_Ecotype_name.RData")

# 删去重复数<3的组合
Ecotype_name <- Ecotype_name[-c(8, 94, 56)]

#--------LA---------
# 替换
c <- c("Col-1" ,"Col-2" ,"Col-3" ,"Col-4" ,"Col-5" ,"Col-6" ,"Col-7" ,"Col-8" ,"Col-9")
LA[LA$Genotype %in% c,]$Genotype <- str_replace(LA[LA$Genotype %in% c,]$Genotype, "Col-1" , "Col-01")
LA[LA$Genotype %in% c,]$Genotype <- str_replace(LA[LA$Genotype %in% c,]$Genotype, "Col-2" , "Col-02")
LA[LA$Genotype %in% c,]$Genotype <- str_replace(LA[LA$Genotype %in% c,]$Genotype, "Col-3" , "Col-03")
LA[LA$Genotype %in% c,]$Genotype <- str_replace(LA[LA$Genotype %in% c,]$Genotype, "Col-4" , "Col-04")
LA[LA$Genotype %in% c,]$Genotype <- str_replace(LA[LA$Genotype %in% c,]$Genotype, "Col-5" , "Col-05")
LA[LA$Genotype %in% c,]$Genotype <- str_replace(LA[LA$Genotype %in% c,]$Genotype, "Col-6" , "Col-06")
LA[LA$Genotype %in% c,]$Genotype <- str_replace(LA[LA$Genotype %in% c,]$Genotype, "Col-7" , "Col-07")
LA[LA$Genotype %in% c,]$Genotype <- str_replace(LA[LA$Genotype %in% c,]$Genotype, "Col-8" , "Col-08")
LA[LA$Genotype %in% c,]$Genotype <- str_replace(LA[LA$Genotype %in% c,]$Genotype, "Col-9" , "Col-09")

# 循环将每个组合数据存入列表，再收集DE p-value信息
la <- list()
sig_la <- data.frame(CF=numeric(), PF=numeric(), CP=numeric())
for (i in Ecotype_name) {
  j <- paste0("F1_", i)
  
  la[[i]] <- LA %>% filter(Genotype %in% c(i, j))
  
  Col <- LA %>% 
    filter(Genotype == paste0("Col-", str_sub(unique(la[[i]][,"V1"]), 7, 8)))
  
  la[[i]] <- rbind(la[[i]], Col)
  
  P <- as.vector(as.matrix(la[[i]] %>% filter(Genotype == i) %>% select(V4)))
  F1 <- as.vector(as.matrix(la[[i]] %>% filter(Genotype == j) %>% select(V4)))
  C <- as.vector(as.matrix(la[[i]] %>% filter(!Genotype %in% c(i,j)) %>% select(V4)))
  
  CF <- t.test(C, F1, paired = FALSE)
  PF <- t.test(P, F1, paired = FALSE)
  CP <- t.test(C, P, paired = FALSE)
  
  sig_la[i,] <- c(signif(CF$p.value, 3), signif(PF$p.value, 3), signif(CP$p.value, 3))
}


#--------CA---------
# 替换
c <- c("Col-1" ,"Col-2" ,"Col-3" ,"Col-4" ,"Col-5" ,"Col-6" ,"Col-7" ,"Col-8" ,"Col-9")
CA[CA$Genotype %in% c,]$Genotype <- str_replace(CA[CA$Genotype %in% c,]$Genotype, "Col-1" , "Col-01")
CA[CA$Genotype %in% c,]$Genotype <- str_replace(CA[CA$Genotype %in% c,]$Genotype, "Col-2" , "Col-02")
CA[CA$Genotype %in% c,]$Genotype <- str_replace(CA[CA$Genotype %in% c,]$Genotype, "Col-3" , "Col-03")
CA[CA$Genotype %in% c,]$Genotype <- str_replace(CA[CA$Genotype %in% c,]$Genotype, "Col-4" , "Col-04")
CA[CA$Genotype %in% c,]$Genotype <- str_replace(CA[CA$Genotype %in% c,]$Genotype, "Col-5" , "Col-05")
CA[CA$Genotype %in% c,]$Genotype <- str_replace(CA[CA$Genotype %in% c,]$Genotype, "Col-6" , "Col-06")
CA[CA$Genotype %in% c,]$Genotype <- str_replace(CA[CA$Genotype %in% c,]$Genotype, "Col-7" , "Col-07")
CA[CA$Genotype %in% c,]$Genotype <- str_replace(CA[CA$Genotype %in% c,]$Genotype, "Col-8" , "Col-08")
CA[CA$Genotype %in% c,]$Genotype <- str_replace(CA[CA$Genotype %in% c,]$Genotype, "Col-9" , "Col-09")

# 循环将每个组合数据存入列表，再收集DE p-value信息
ca <- list()
sig_ca <- data.frame(CF=numeric(), PF=numeric(), CP=numeric())
for (i in Ecotype_name) {
  j <- paste0("F1_", i)
  
  ca[[i]] <- CA %>% filter(Genotype %in% c(i, j))
  
  Col <- CA %>% 
    filter(Genotype == paste0("Col-", str_sub(unique(ca[[i]][,"V1"]), 7, 8)))
  
  ca[[i]] <- rbind(ca[[i]], Col)
  
  P <- as.vector(as.matrix(ca[[i]] %>% filter(Genotype == i) %>% select(V4)))
  F1 <- as.vector(as.matrix(ca[[i]] %>% filter(Genotype == j) %>% select(V4)))
  C <- as.vector(as.matrix(ca[[i]] %>% filter(!Genotype %in% c(i,j)) %>% select(V4)))
  
  CF <- t.test(C, F1, paired = FALSE)
  PF <- t.test(P, F1, paired = FALSE)
  CP <- t.test(C, P, paired = FALSE)
  
  sig_ca[i,] <- c(signif(CF$p.value, 3), signif(PF$p.value, 3), signif(CP$p.value, 3))
}


#--------CN---------
# 替换
c <- c("Col-1" ,"Col-2" ,"Col-3" ,"Col-4" ,"Col-5" ,"Col-6" ,"Col-7" ,"Col-8" ,"Col-9")
CN[CN$Genotype %in% c,]$Genotype <- str_replace(CN[CN$Genotype %in% c,]$Genotype, "Col-1" , "Col-01")
CN[CN$Genotype %in% c,]$Genotype <- str_replace(CN[CN$Genotype %in% c,]$Genotype, "Col-2" , "Col-02")
CN[CN$Genotype %in% c,]$Genotype <- str_replace(CN[CN$Genotype %in% c,]$Genotype, "Col-3" , "Col-03")
CN[CN$Genotype %in% c,]$Genotype <- str_replace(CN[CN$Genotype %in% c,]$Genotype, "Col-4" , "Col-04")
CN[CN$Genotype %in% c,]$Genotype <- str_replace(CN[CN$Genotype %in% c,]$Genotype, "Col-5" , "Col-05")
CN[CN$Genotype %in% c,]$Genotype <- str_replace(CN[CN$Genotype %in% c,]$Genotype, "Col-6" , "Col-06")
CN[CN$Genotype %in% c,]$Genotype <- str_replace(CN[CN$Genotype %in% c,]$Genotype, "Col-7" , "Col-07")
CN[CN$Genotype %in% c,]$Genotype <- str_replace(CN[CN$Genotype %in% c,]$Genotype, "Col-8" , "Col-08")
CN[CN$Genotype %in% c,]$Genotype <- str_replace(CN[CN$Genotype %in% c,]$Genotype, "Col-9" , "Col-09")

# 循环将每个组合数据存入列表，再收集DE p-value信息
cn <- list()
sig_cn <- data.frame(CF=numeric(), PF=numeric(), CP=numeric())
for (i in Ecotype_name) {
  j <- paste0("F1_", i)
  
  cn[[i]] <- CN %>% filter(Genotype %in% c(i, j))
  
  Col <- CN %>% 
    filter(Genotype == paste0("Col-", str_sub(unique(cn[[i]][,"V1"]), 7, 8)))
  
  cn[[i]] <- rbind(cn[[i]], Col)
  
  P <- as.vector(as.matrix(cn[[i]] %>% filter(Genotype == i) %>% select(CN)))
  F1 <- as.vector(as.matrix(cn[[i]] %>% filter(Genotype == j) %>% select(CN)))
  C <- as.vector(as.matrix(cn[[i]] %>% filter(!Genotype %in% c(i,j)) %>% select(CN)))
  
  CF <- t.test(C, F1, paired = FALSE)
  PF <- t.test(P, F1, paired = FALSE)
  CP <- t.test(C, P, paired = FALSE)
  
  sig_cn[i,] <- c(signif(CF$p.value, 3), signif(PF$p.value, 3), signif(CP$p.value, 3))
}

# save
save(sig_la, file = "sig_la.RData")
save(sig_ca, file = "sig_ca.RData")
save(sig_cn, file = "sig_cn.RData")

# DE relationship
load("/home/ug0096/DE_heatmap/sig_ca.RData")
load("/home/ug0096/DE_heatmap/sig_cn.RData")
load("/home/ug0096/DE_heatmap/sig_la.RData")
sig_la["Je-0",] <- c(1,1,1);sig_la["Pu2-23",] <- c(1,1,1);sig_la["Abd-0",] <- c(1,1,1)
sig_ca["Je-0",] <- c(1,1,1);sig_ca["Pu2-23",] <- c(1,1,1);sig_ca["Abd-0",] <- c(1,1,1)
sig_cn["Je-0",] <- c(1,1,1);sig_cn["Pu2-23",] <- c(1,1,1);sig_cn["Abd-0",] <- c(1,1,1)
load("/home/ug0096/96pop_Ecotype_name.RData")
sig_la <- sig_la[Ecotype_name,]
sig_ca <- sig_ca[Ecotype_name,]
sig_cn <- sig_cn[Ecotype_name,]
load("/home/ug0096/part0.phenotype/trait_32.RData")

data <- trait_32[Ecotype_name,c(5,6,9,10,1,2)]

data$LA_Col <- 1
data$CN_Col <- 1
data$CA_Col <- 1

# class = clear
data$LA_vs <- "NA"
data$CN_vs <- "NA"
data$CA_vs <- "NA"

# input class
for (i in 1:96) {
  ## LA
  if(data[i,"LA_F1"] > data[i,"LA_Col"] & data[i,"LA_Col"] > data[i,"LA_Ecotype"]){
    data[i,"LA_vs"] <- "F>C>E"
  }
  else if(data[i,"LA_F1"] > data[i,"LA_Ecotype"] & data[i,"LA_Ecotype"] > data[i,"LA_Col"]){
    data[i,"LA_vs"] <- "F>E>C"
  }
  else if(data[i,"LA_Col"] > data[i,"LA_F1"] & data[i,"LA_F1"] > data[i,"LA_Ecotype"]){
    data[i,"LA_vs"] <- "C>F>E"
  }
  else if(data[i,"LA_Col"] > data[i,"LA_Ecotype"] & data[i,"LA_Ecotype"] > data[i,"LA_F1"]){
    data[i,"LA_vs"] <- "C>E>F"
  }
  else if(data[i,"LA_Ecotype"] > data[i,"LA_F1"] & data[i,"LA_F1"] > data[i,"LA_Col"]){
    data[i,"LA_vs"] <- "E>F>C"
  }
  else if(data[i,"LA_Ecotype"] > data[i,"LA_Col"] & data[i,"LA_Col"] > data[i,"LA_F1"]){
    data[i,"LA_vs"] <- "E>C>F"
  }
  ## CN
  if(data[i,"CN_F1"] > data[i,"CN_Col"] & data[i,"CN_Col"] > data[i,"CN_Ecotype"]){
    data[i,"CN_vs"] <- "F>C>E"
  }
  else if(data[i,"CN_F1"] > data[i,"CN_Ecotype"] & data[i,"CN_Ecotype"] > data[i,"CN_Col"]){
    data[i,"CN_vs"] <- "F>E>C"
  }
  else if(data[i,"CN_Col"] > data[i,"CN_F1"] & data[i,"CN_F1"] > data[i,"CN_Ecotype"]){
    data[i,"CN_vs"] <- "C>F>E"
  }
  else if(data[i,"CN_Col"] > data[i,"CN_Ecotype"] & data[i,"CN_Ecotype"] > data[i,"CN_F1"]){
    data[i,"CN_vs"] <- "C>E>F"
  }
  else if(data[i,"CN_Ecotype"] > data[i,"CN_F1"] & data[i,"CN_F1"] > data[i,"CN_Col"]){
    data[i,"CN_vs"] <- "E>F>C"
  }
  else if(data[i,"CN_Ecotype"] > data[i,"CN_Col"] & data[i,"CN_Col"] > data[i,"CN_F1"]){
    data[i,"CN_vs"] <- "E>C>F"
  }
  ## CA
  if(data[i,"CA_F1"] > data[i,"CA_Col"] & data[i,"CA_Col"] > data[i,"CA_Ecotype"]){
    data[i,"CA_vs"] <- "F>C>E"
  }
  else if(data[i,"CA_F1"] > data[i,"CA_Ecotype"] & data[i,"CA_Ecotype"] > data[i,"CA_Col"]){
    data[i,"CA_vs"] <- "F>E>C"
  }
  else if(data[i,"CA_Col"] > data[i,"CA_F1"] & data[i,"CA_F1"] > data[i,"CA_Ecotype"]){
    data[i,"CA_vs"] <- "C>F>E"
  }
  else if(data[i,"CA_Col"] > data[i,"CA_Ecotype"] & data[i,"CA_Ecotype"] > data[i,"CA_F1"]){
    data[i,"CA_vs"] <- "C>E>F"
  }
  else if(data[i,"CA_Ecotype"] > data[i,"CA_F1"] & data[i,"CA_F1"] > data[i,"CA_Col"]){
    data[i,"CA_vs"] <- "E>F>C"
  }
  else if(data[i,"CA_Ecotype"] > data[i,"CA_Col"] & data[i,"CA_Col"] > data[i,"CA_F1"]){
    data[i,"CA_vs"] <- "E>C>F"
  }
}



# PLOT
plot_m1 <- data[,c(10,11,12)]

colnames(sig_la) <- paste0("LA_", colnames(sig_la))
colnames(sig_cn) <- paste0("CN_", colnames(sig_cn))
colnames(sig_ca) <- paste0("CA_", colnames(sig_ca))

de <- cbind(plot_m1, sig_la, sig_cn, sig_ca)
pe <- de

## de
for (i in c(1:96)) {
  # c=e=f
  if (de[i, "CA_CF"] > 0.05 & de[i, "CA_PF"] > 0.05 & de[i, "CA_CP"] > 0.05) {
    de[i, "CA_vs"] <- "F=C=E"
  }
  if (de[i, "CN_CF"] > 0.05 & de[i, "CN_PF"] > 0.05 & de[i, "CN_CP"] > 0.05) {
    de[i, "CN_vs"] <- "F=C=E"
  }
  if (de[i, "LA_CF"] > 0.05 & de[i, "LA_PF"] > 0.05 & de[i, "LA_CP"] > 0.05) {
    de[i, "LA_vs"] <- "F=C=E"
  }
  
  # c=e
  if (de[i, "CA_CF"] < 0.05 & de[i, "CA_PF"] < 0.05 & de[i, "CA_CP"] > 0.05) {
    de[i, "CA_vs"] <- str_replace(de[i, "CA_vs"], "C>E", "C=E")
    de[i, "CA_vs"] <- str_replace(de[i, "CA_vs"], "E>C", "C=E")
  }
  if (de[i, "CN_CF"] < 0.05 & de[i, "CN_PF"] < 0.05 & de[i, "CN_CP"] > 0.05) {
    de[i, "CN_vs"] <- str_replace(de[i, "CN_vs"], "C>E", "C=E")
    de[i, "CN_vs"] <- str_replace(de[i, "CN_vs"], "E>C", "C=E")
  }
  if (de[i, "LA_CF"] < 0.05 & de[i, "LA_PF"] < 0.05 & de[i, "LA_CP"] > 0.05) {
    de[i, "LA_vs"] <- str_replace(de[i, "LA_vs"], "C>E", "C=E")
    de[i, "LA_vs"] <- str_replace(de[i, "LA_vs"], "E>C", "C=E")
  }
  
  # c=f
  if (de[i, "CA_CF"] > 0.05 & de[i, "CA_PF"] < 0.05 & de[i, "CA_CP"] < 0.05) {
    de[i, "CA_vs"] <- str_replace(de[i, "CA_vs"], "C>F", "C=F")
    de[i, "CA_vs"] <- str_replace(de[i, "CA_vs"], "F>C", "C=F")
  }
  if (de[i, "CN_CF"] > 0.05 & de[i, "CN_PF"] < 0.05 & de[i, "CN_CP"] < 0.05) {
    de[i, "CN_vs"] <- str_replace(de[i, "CN_vs"], "C>F", "C=F")
    de[i, "CN_vs"] <- str_replace(de[i, "CN_vs"], "F>C", "C=F")
  }
  if (de[i, "LA_CF"] > 0.05 & de[i, "LA_PF"] < 0.05 & de[i, "LA_CP"] < 0.05) {
    de[i, "LA_vs"] <- str_replace(de[i, "LA_vs"], "C>F", "C=F")
    de[i, "LA_vs"] <- str_replace(de[i, "LA_vs"], "F>C", "C=F")
  }
  
  # e=f
  if (de[i, "CA_CF"] < 0.05 & de[i, "CA_PF"] > 0.05 & de[i, "CA_CP"] < 0.05) {
    de[i, "CA_vs"] <- str_replace(de[i, "CA_vs"], "E>F", "E=F")
    de[i, "CA_vs"] <- str_replace(de[i, "CA_vs"], "F>E", "E=F")
  }
  if (de[i, "CN_CF"] < 0.05 & de[i, "CN_PF"] > 0.05 & de[i, "CN_CP"] < 0.05) {
    de[i, "CN_vs"] <- str_replace(de[i, "CN_vs"], "E>F", "E=F")
    de[i, "CN_vs"] <- str_replace(de[i, "CN_vs"], "F>E", "E=F")
  }
  if (de[i, "LA_CF"] < 0.05 & de[i, "LA_PF"] > 0.05 & de[i, "LA_CP"] < 0.05) {
    de[i, "LA_vs"] <- str_replace(de[i, "LA_vs"], "E>F", "E=F")
    de[i, "LA_vs"] <- str_replace(de[i, "LA_vs"], "F>E", "E=F")
  }
  
  # ==
  if (isTRUE(de[i, "CA_CF"] < 0.05 & de[i, "CA_PF"] > 0.05 & de[i, "CA_CP"] > 0.05) | isTRUE(de[i, "CA_CF"] > 0.05 & de[i, "CA_PF"] < 0.05 & de[i, "CA_CP"] > 0.05) | isTRUE(de[i, "CA_CF"] > 0.05 & de[i, "CA_PF"] > 0.05 & de[i, "CA_CP"] < 0.05)) {
    de[i, "CA_vs"] <- "F=C=E"
  }
  if (isTRUE(de[i, "CN_CF"] < 0.05 & de[i, "CN_PF"] > 0.05 & de[i, "CN_CP"] > 0.05) | isTRUE(de[i, "CN_CF"] > 0.05 & de[i, "CN_PF"] < 0.05 & de[i, "CN_CP"] > 0.05) | isTRUE(de[i, "CN_CF"] > 0.05 & de[i, "CN_PF"] > 0.05 & de[i, "CN_CP"] < 0.05)) {
    de[i, "CN_vs"] <- "F=C=E"
  }
  if (isTRUE(de[i, "LA_CF"] < 0.05 & de[i, "LA_PF"] > 0.05 & de[i, "LA_CP"] > 0.05) | isTRUE(de[i, "LA_CF"] > 0.05 & de[i, "LA_PF"] < 0.05 & de[i, "LA_CP"] > 0.05) | isTRUE(de[i, "LA_CF"] > 0.05 & de[i, "LA_PF"] > 0.05 & de[i, "LA_CP"] < 0.05)) {
    de[i, "LA_vs"] <- "F=C=E"
  }
}

## pe
pe[, c(1:3)] <- "*"

for (i in c(1:96)) {
  # ==
  if (isTRUE(pe[i, "CA_CF"] < 0.05 & pe[i, "CA_PF"] > 0.05 & pe[i, "CA_CP"] > 0.05) | isTRUE(pe[i, "CA_CF"] > 0.05 & pe[i, "CA_PF"] < 0.05 & pe[i, "CA_CP"] > 0.05) | isTRUE(pe[i, "CA_CF"] > 0.05 & pe[i, "CA_PF"] > 0.05 & pe[i, "CA_CP"] < 0.05)) {
    pe[i, "CA_vs"] <- "#"
  }
  if (isTRUE(pe[i, "CN_CF"] < 0.05 & pe[i, "CN_PF"] > 0.05 & pe[i, "CN_CP"] > 0.05) | isTRUE(pe[i, "CN_CF"] > 0.05 & pe[i, "CN_PF"] < 0.05 & pe[i, "CN_CP"] > 0.05) | isTRUE(pe[i, "CN_CF"] > 0.05 & pe[i, "CN_PF"] > 0.05 & pe[i, "CN_CP"] < 0.05)) {
    pe[i, "CN_vs"] <- "#"
  }
  if (isTRUE(pe[i, "LA_CF"] < 0.05 & pe[i, "LA_PF"] > 0.05 & pe[i, "LA_CP"] > 0.05) | isTRUE(pe[i, "LA_CF"] > 0.05 & pe[i, "LA_PF"] < 0.05 & pe[i, "LA_CP"] > 0.05) | isTRUE(pe[i, "LA_CF"] > 0.05 & pe[i, "LA_PF"] > 0.05 & pe[i, "LA_CP"] < 0.05)) {
    pe[i, "LA_vs"] <- "#"
  }
}

# bind de and pe
dp <- cbind(de[,1:3], pe[,1:3])
colnames(dp) <- c("LA_vs", "CN_vs", "CA_vs", "LA_mark", "CN_mark", "CA_mark")
write.table(dp, file = "dp.txt", sep = "\t", quote = F, row.names = T)


## 确定行的顺序,并删除3个：Abd-0, Je-0, Pu2-23
order <- c("Ven-1","Vie-0","Chat-1","Er-0","Is-0","Kil-0","C24","Per-1","Bor-4","Ta-0","Tu-0","TueV-13", "HKT2.4","Ka-0","TueWa1-2","Ak-1","Westkar-4","Wl-0","Ws-2","Nok-3","Nw-0","Pna-17","ICE163","Sp-0","TDr-1","ICE150","ICE29","Dr-0", "Kro-0","Ler-1","Lm-2","Old-1","Pog-0","Sha","Stw-0","Ber","Kyoto","Le-0","Ob-0","Ts-1","Rubezhnoe-1","Zal-1","Koch-1","Jl-3","ICE119","Bch-1","Bla-1","Nemrut-1","Bs-1","Com-1","Ga-0","Vind-1","Ey15-2","Kar-1","Hovdala-2","Litva","Wei-0","Pu2-7","Knox-18","Fei-0","Yeg-1","Baa-1","Pla-0","Me-0","Rd-0","Pt-0","HR-5","Da1-12","El-0","Rou-0","Uk-1","Van-0","Ra-0", "Uod-1","Kin-0","Gel-1","Gu-0","RRS-7","Tuescha9","Bl-1","Tscha-1","Si-0","Pro-0","La-0","Cerv-1","Bay-0","ICE102","No-0","Ha-0", "Tottarp-2","Kl-5","Hn-0", "Mer-6")

## add heterosis
plot_m2 <- trait_32[order,c(8,7,12,11,4,3)]

## plot
### global plot settings by complexheatmap
ht_global_opt(
  legend_border = "black",
  heatmap_border = TRUE,
  annotation_border = TRUE)

h <- Heatmap(dp[order,1:3],
             show_column_names = T,
             show_row_names = T,
             # 聚类
             cluster_columns = F,
             cluster_rows = F,
             # 图例加名字
             name = "comparisons",
             # 加列标题
             column_title = "phenotype",
             column_title_side = "top",
             # 定义离散值颜色映射
             col = structure(c("#F6511D","#F6511D","#F6511D", "#DCDCDC", "#FFB400","#FFB400", "#9195D2","#9195D2", "#00A6ED","#00A6ED"), 
             names=c("F>C>E","F>E>C","F>C=E", "F=C=E", "C=F>E","E=F>C", "C>F>E","E>F>C", "E>C=F","C>E=F")),
             # 外边框
             border = 'black',
             # 网格线
             rect_gp = gpar(col = "black", lwd = 0.2),
             # *和#
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf(dp[order, 4:6][i, j]), x, y, gp = gpar(fontsize = 10))
             }
) +
  Heatmap(plot_m2,
          name = 'Heterosis level',
          col = colorRamp2(c(-0.4, 0, 1.2), c("#247BA0", "white", "#F25F5C")),
          show_column_names = T,
          show_row_names = T,
          cluster_columns = F,
          column_title = 'Heterosis',
          column_title_side = "top",
          # 外边框
          border = 'black',
          # 网格线
          rect_gp = gpar(col = "black", lwd = 0.2))


pdf(file = "Fig2b_phenotype_comparisons_heatmap.pdf", width = 5, height = 15)
draw(h, heatmap_legend_side = "right")
dev.off()
