setwd("/data2/usr/LiuWW/project_3/part9.twas/TWAS_TL/")

# Load nececarry packages
library(data.table)
library(tidyverse)
library(rtracklayer)
library(scales)

######################
# TWAS data prep #
######################

load("/data2/usr/LiuWW/project_3/part9.twas/TWAS_list.RData")

###########################################################################
# Prepare genotype data in rMVP format 主要是用K矩阵和Q矩阵
## The genotypic and phenotypic data can be formatted in one step using the code below. This chunk of code also calculates the kinship matrix, principle components, and saves genotypic data as a filebacked matrix

# If you have genotype data in PLINK Binary format:  
# fileBed, name of genotype data in PLINK Binary format
# fileKin, TRUE or FALSE, if TRUE, kinship matrix represents relationship among individuals will be calculated
# filePC, TRUE or FALSE, if TRUE, principal component analysis will be performed
# out, prefix of output file
# priority, "speed" or "memory", the "speed" mode is faster but uses more memory while "memory" is slower but uses less memory
# maxLine, number, if priority = "memory", it is the number of markers read into memory

# Full-featured function (Recommended)
MVP.Data(fileBed="/data2/usr/LiuWW/project_2/part4.GWAS/06.GWAS_GEMMA_96pop/snp_bfile",
         filePhe=NULL,
         fileKin=TRUE,
         filePC=TRUE,       
         #priority="speed",
         #maxLine=10000,
         out="mvp.plink"
)
###########################################################################

############################################################################
# Read in 96 expression matrix generated
## expression data
load("/data2/usr/LiuWW/project_3/part7.downstream_analysis/data/TPM_genes_high_TL_193.RData")

## expression data prepare
### 数据保留三位小数
TPM_genes_high_TL_193 <- round(as.matrix(TPM_genes_high_TL_193), digits = 3)

### split
#### Col-0
TPM_genes_high_TL_Col <- as.matrix(TPM_genes_high_TL_193[,193])

TPM_genes_high_TL_Col_96 <- TPM_genes_high_TL_Col

for (i in c(2:96)) {
  TPM_genes_high_TL_Col_96 <- cbind(TPM_genes_high_TL_Col_96, TPM_genes_high_TL_Col)
}

#### Ecotype
eco <- seq(from = 1, by = 2, length = 96)
TPM_genes_high_TL_eco_96 <- TPM_genes_high_TL_193[,eco]

################################
# 获得生态型名称顺序
################################
Ecotype_name <- str_sub(colnames(TPM_genes_high_TL_eco_96), start = 19, end = -1)

#### F1
F1 <- seq(from = 2, by = 2, length = 96)
TPM_genes_high_TL_F1_96 <- TPM_genes_high_TL_193[,F1]

#### F1_vs_MP
# TPM_genes_high_TL_F1_vs_MP_96 <- TPM_genes_high_TL_F1_96/((TPM_genes_high_TL_eco_96 + TPM_genes_high_TL_Col_96)/2)

#### F1_vs_HP
# TPM_genes_high_TL_HP_96 <- TPM_genes_high_TL_Col_96
# for (i in c(1:21014)) {
#  for (j in c(1:96)) {
#    TPM_genes_high_TL_HP_96[i, j] <- max(TPM_genes_high_TL_eco_96[i,j], TPM_genes_high_TL_Col_96[i,j])
#  }
# }
# TPM_genes_high_TL_F1_vs_HP_96 <- TPM_genes_high_TL_F1_96/TPM_genes_high_TL_HP_96

#### F1_vs_PP
TPM_genes_high_TL_F1_vs_PP_96 <- TPM_genes_high_TL_F1_96/TPM_genes_high_TL_eco_96

#### F1_vs_Col
TPM_genes_high_TL_F1_vs_Col_96 <- TPM_genes_high_TL_F1_96/TPM_genes_high_TL_Col_96

### 修改行名列名，方便绘图
colnames(TPM_genes_high_TL_Col_96) <- Ecotype_name
colnames(TPM_genes_high_TL_eco_96) <- Ecotype_name
colnames(TPM_genes_high_TL_F1_96) <- Ecotype_name
# colnames(TPM_genes_high_TL_F1_vs_MP_96) <- Ecotype_name
# colnames(TPM_genes_high_TL_F1_vs_HP_96) <- Ecotype_name
colnames(TPM_genes_high_TL_F1_vs_PP_96) <- Ecotype_name
colnames(TPM_genes_high_TL_F1_vs_Col_96) <- Ecotype_name

rownames(TPM_genes_high_TL_F1_96) <- paste0(rownames(TPM_genes_high_TL_eco_96), "_F1")
# rownames(TPM_genes_high_TL_F1_vs_MP_96) <- paste0(rownames(TPM_genes_high_TL_eco_96), "_F1_vs_MP")
# rownames(TPM_genes_high_TL_F1_vs_HP_96) <- paste0(rownames(TPM_genes_high_TL_eco_96), "_F1_vs_HP")
rownames(TPM_genes_high_TL_F1_vs_PP_96) <- paste0(rownames(TPM_genes_high_TL_eco_96), "_F1_vs_PP")
rownames(TPM_genes_high_TL_F1_vs_Col_96) <- paste0(rownames(TPM_genes_high_TL_eco_96), "_F1_vs_Col")
rownames(TPM_genes_high_TL_Col_96) <- paste0(rownames(TPM_genes_high_TL_Col_96), "_Col")
rownames(TPM_genes_high_TL_eco_96) <- paste0(rownames(TPM_genes_high_TL_eco_96), "_eco")
############################################################################

############################################################################
# gtf文件处理得到gene loc
# Keep only protein coding genes gff file
gtf <- rtracklayer::import('/data2/usr/LiuWW/project_3/part0.ref/genes.gtf')
gtf_df <- as.data.frame(gtf)
gff_gene <- gtf_df[gtf_df$type=="gene",] # Protein coding genes

# Keep needed columns from gff file
gff_keep_cols = c("seqnames","source","type","start","end","strand","gene_id","gene_biotype")
out = gff_gene[,colnames(gff_gene)%in%gff_keep_cols]
genes_keep = out$gene_id
length(unique(genes_keep)) # 32,833 protein coding genes

# Gene coordinates of expression matrix
gene_loc = out[,c("gene_id","seqnames","start","end")]
colnames(gene_loc) = c("geneid","chr","s1","s2")
############################################################################


############################################################################
# 循环内容：每一个表达表型对每一个测量表型做TWAS
# 表达表型有5种：
EE <- c("TPM_genes_high_TL_eco_96",
        "TPM_genes_high_TL_F1_96",
        "TPM_genes_high_TL_F1_vs_Col_96",
        "TPM_genes_high_TL_F1_vs_PP_96")

# 测量表型根据TWAS_list确定

for (E in EE) {
  # 每种表达表型建一个目录,用于保存输出结果
  ifelse(!dir.exists(E), dir.create(E), FALSE)
  
  # 初始表达数据
  #----------------------
  expr <- get(E)
  expr <- rownames_to_column(as.data.frame(expr), var = "id")
  expr$id <- str_sub(string = expr$id, start = 1, end = 9)
  rownames(expr) <- expr$id
  #----------------------
  
  
  # 对表达数据进行过滤等处理
  # Keep only protein coding genes from .gff file
  expr_sub = expr[expr$id%in%gene_loc$geneid,]
  
  # Check variance across all lines for expression in each data file
  variance = data.frame(apply(expr_sub[,2:ncol(expr_sub)], 1, var))
  variance$names = rownames(variance)
  colnames(variance) = c("Var", "Gene")
  
  # Remove genes with zero variance
  noVar = subset(variance, Var == 0)
  noVar.genes = noVar$Gene
  dim(expr_sub)
  expr_sub = expr_sub[!expr_sub$id%in%noVar.genes,]
  dim(expr_sub) # 20,989 genes
  
  # Remove genes not expressed in at least 5% of lines
  low_expr = data.frame(apply(expr_sub[,2:ncol(expr_sub)], 1, function(x){sum(x == 0)}))
  low_expr$names = rownames(low_expr)
  colnames(low_expr) = c("Low_expr", "Gene")
  Low = subset(low_expr, Low_expr > (ncol(expr_sub)-1)*0.95)
  Low.genes = Low$Gene
  expr_sub = expr_sub[!rownames(expr_sub)%in%Low.genes,]
  dim(expr_sub) # 20,883 genes
  
  ############################################################################
  # Read in phenotypic data
  load(file = "/data2/usr/LiuWW/project_2/part0.phenotype/trait_32.RData")
  phenos <- trait_32[Ecotype_name, c(1,2,5,6,9,10,22,24,26,28,30,32)]
  
  colnames(phenos) <- c("CA_Ecotype", "CA_F1", "LA_Ecotype", "LA_F1", "CN_Ecotype", "CN_F1", "CA_F1vsC", "CA_F1vsP", "LA_F1vsC", "LA_F1vsP", "CN_F1vsC", "CN_F1vsP")
  
  phenos$Ecotype <- rownames(phenos)
  
  phenos <- dplyr::select(phenos, Ecotype, everything())
  ############################################################################
  
  # Match name order to phenotype file
  phenos_order  = c("id",phenos$Ecotype)
  expr_sub = expr_sub[,order(match(colnames(expr_sub),phenos_order))]
  identical(as.character(colnames(expr_sub)[2:ncol(expr_sub)]), as.character(phenos$Ecotype))
  gene_loc = gene_loc[gene_loc$geneid%in%expr_sub$id,]
  
  # Make quantile rank normalized gene expression matrix
  mat = expr_sub[,2:ncol(expr_sub)]
  dim(mat)
  mat = t(apply(mat, 1, rank, ties.method = "average",na.last = "keep"))
  mat = qnorm(mat /ncol(expr_sub))
  
  # Scale TPM values to be between -1 and 1 to be compatible with rrBLUP
  mat = t(apply(mat, 1, rescale, to = c(-1,1)))
  mat = cbind(expr_sub$id, mat)
  colnames(mat)[1] = "id"
  
  # Make output file
  mat_out = left_join(gene_loc, mat, by = c("geneid"="id"), copy = TRUE)
  
  # Keep only genes on chromosomes, remove mitochondrial and chrolorplast genes
  chrm_keep = as.factor(seq(1,5,1))
  mat_out = mat_out[mat_out$chr%in%chrm_keep,]
  dim(mat_out) # 20,738 genes
  
  # Sort by chromosome and transcription start site from gff file
  mat_out = mat_out %>%
    arrange(chr,s1)
  
  # Write the file
  fwrite(mat_out, file = paste0(E, "/", E, "_scaled_quant_rank_TPM.txt"), col.names = TRUE, sep = "\t", quote = FALSE)
  
  # Look at histogram of expression of AT2G20370
  # hist(as.numeric(mat_out[mat_out$geneid == "AT2G20370", 5:ncol(mat_out)]))
  
  ######################
  # TWAS analysis #
  ######################
  # Load nececarry packages
  library(data.table)
  library(rrBLUP)
  library(tidyverse)
  library(rMVP)
  
  # Read in phenotypes
  colnames(phenos)[1] <- "Taxa"
  
  # Read in expression data and put in the format for rrBLUP
  expr = fread(paste0(E, "/", E, "_scaled_quant_rank_TPM.txt"), data.table = FALSE)
  expr = expr[,-4]
  colnames(expr)[1:3] = c("marker","chrom","pos")
  
  # Read A matrix
  A.widiv = attach.big.matrix("./mvp.plink.kin.desc")
  A.widiv = A.widiv[,]
  mean(diag(A.widiv))  # = 2 because inbred lines 自己的数据是1.891521
  
  # Check files are in the same order
  identical(as.character(phenos$Taxa), as.character(colnames(expr[,4:ncol(expr)])))
  
  # Run rrBLUP for each corresponding trait in TWAS_list and write out results
  pheno_names <- TWAS_list[TWAS_list$V1==E,2]
  
  ## K only model
  for(i in pheno_names){
    fitlist = list()
    fitlist[[i]] = GWAS(pheno = phenos[,c("Taxa",i)], geno = expr, K=A.widiv, min.MAF=0.0, n.core=10, P3D=TRUE, plot=FALSE)
    tmp = fitlist[[i]]
    colnames(tmp) = c("geneid","chrom","pos","p_value_log10")
    tmp$trait = i
    fwrite(tmp, file = paste0(E, "/", E, "_", i, "_rrBLUP_TWAS_K", ".txt"),
           col.names = TRUE, sep = "\t", quote = FALSE, na = "NA")
    rm(fitlist)
  }
  
  # ## K and SNP calculated PCs from K matrix in model
  # for(i in pheno_names){
  #   fitlist = list()
  #   fitlist[[i]] = GWAS(pheno = phenos[,c("Taxa",i)], geno = expr, n.PC = 3, K=A.widiv, min.MAF=0.0, n.core=10, P3D=TRUE, plot=FALSE)
  #   tmp = fitlist[[i]]
  #   colnames(tmp) = c("geneid","chrom","pos","p_value_log10")
  #   tmp$trait = i
  #   fwrite(tmp, file = paste0(E, "/", E, "_", i, "_rrBLUP_TWAS_KQ", ".txt"),
  #          col.names = TRUE, sep = "\t", quote = FALSE, na = "NA")
  #   rm(fitlist)
  # }
}



###################
# 画曼哈顿图
###################
## # 每种表达表型对TWAS_list里相应的测量表型的TWAS结果统计表合并为一个大表,并绘图
for (E in EE) {
  pheno_names <- TWAS_list[TWAS_list$V1==E,2]
  ## K only model
  ### get list
  tmp <- fread(paste0(E, "/", E, "_", pheno_names[1], "_rrBLUP_TWAS_K", ".txt"), data.table = FALSE)
  tmp <- tmp[,1:4]
  
  for(i in pheno_names[2:length(pheno_names)]){
    right <- fread(paste0(E, "/", E, "_", i, "_rrBLUP_TWAS_K", ".txt"), data.table = FALSE)
    right <- right[,1:4]
    tmp <- left_join(tmp, right, by = c("geneid", "chrom", "pos"))
  }
  
  colnames(tmp) = c("geneid","chrom","pos", pheno_names)
  
  # f函数用于将-log10(P)转换为P
  f <- function(x) {
    10^(-x)
  }
  # 执行f函数
  tmp <-as.matrix(tmp)
  for (i in 4:(3+length(pheno_names))) {
    for (j in 1:nrow(tmp)) {
      tmp[j,i] <- f(as.numeric(tmp[j,i]))
    }
  }
  
  fwrite(tmp, file = paste0(E, "/", E, "_", "all", "_rrBLUP_TWAS_K", ".txt"),
         col.names = TRUE, sep = "\t", quote = FALSE, na = "NA")
  
  ### get plot
  # MVP.Report(tmp,plot.type="c",r=5,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:5),sep=""),
  #            threshold=c(0.01/nrow(tmp),0.05/nrow(tmp)),cir.chr.h=15,amplify=TRUE,threshold.lty=c(1,2),
  #            threshold.col=c("red", "blue"),signal.line=1,signal.col=c("red","green"),
  #            chr.den.col=c("darkgreen","yellow","red"),
  #            bin.size=1e6,outward=FALSE,
  #            outpath = E,
  #            file.type="jpg",dpi = 300,
  #            memo="K",
  #            cir.band=4,
  #            H=15,
  #            cir.legend=TRUE,
  #            height = 14,width = 14)
  MVP.Report(tmp, plot.type = "m", 
             col = c("#4197d8", "#f8c120", "#413496", "#495226", "#d60b6f", 
                     "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d"),
             chr.labels=paste("Chr",c(1:5),sep=""), 
             threshold=c(0.01/nrow(tmp),0.05/nrow(tmp),0.01),
             threshold.col=c("red", "red", "red"),
             threshold.lty=c(2),
             signal.line=1,signal.col=c("red","red","red"),
             chr.den.col=c("darkgreen","yellow","red"),
             bin.size=1e6,outward=FALSE,
             outpath = E,
             file.type="pdf",
             dpi = 300,
             memo="K",
             multracks=TRUE,
             cex.axis = 1.5, 
             cex.lab = 1.5, 
             xlab = "Chromosome",
             ylab = expression(-log[10](italic(p)))
             )
  
  # ## K and SNP calculated PCs from K matrix in model
  # ### get list
  # tmp <- fread(paste0(E, "/", E, "_", pheno_names[1], "_rrBLUP_TWAS_KQ", ".txt"), data.table = FALSE)
  # tmp <- tmp[,1:4]
  # 
  # for(i in pheno_names[2:length(pheno_names)]){
  #   right <- fread(paste0(E, "/", E, "_", i, "_rrBLUP_TWAS_KQ", ".txt"), data.table = FALSE)
  #   right <- right[,1:4]
  #   tmp <- left_join(tmp, right, by = c("geneid", "chrom", "pos"))
  # }
  # 
  # colnames(tmp) = c("geneid","chrom","pos", pheno_names)
  # 
  # # f函数用于将-log10(P)转换为P
  # f <- function(x) {
  #   10^(-x)
  # }
  # # 执行f函数
  # tmp <-as.matrix(tmp)
  # for (i in 4:(3+length(pheno_names))) {
  #   for (j in 1:nrow(tmp)) {
  #     tmp[j,i] <- f(as.numeric(tmp[j,i]))
  #   }
  # }
  # 
  # fwrite(tmp, file = paste0(E, "/", E, "_", "all", "_rrBLUP_TWAS_KQ", ".txt"),
  #        col.names = TRUE, sep = "\t", quote = FALSE, na = "NA")
  # 
  # ### get plot
  # MVP.Report(tmp,plot.type="c",r=5,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:5),sep=""),
  #            threshold=c(0.01/nrow(tmp),0.05/nrow(tmp)),cir.chr.h=15,amplify=TRUE,threshold.lty=c(1,2),
  #            threshold.col=c("red", "blue"),signal.line=1,signal.col=c("red","green"),
  #            chr.den.col=c("darkgreen","yellow","red"),
  #            bin.size=1e6,outward=FALSE,
  #            outpath = E,
  #            file.type="jpg",dpi = 300,
  #            memo="KQ",
  #            cir.band=4,
  #            H=15,
  #            cir.legend=TRUE,
  #            height = 14,width = 14)
}



###################
# 画显著点的相关性散点图
###################
## # 每种表达表型对相应每种测量表型的TWAS结果统计表筛选显著位点
for (E in EE) {
  pheno_names <- TWAS_list[TWAS_list$V1==E,2]
  for (P in pheno_names) {
    genes_cor <- data.frame(geneid = character(0), expr = character(0), pheno = character(0), cor = numeric(0), pvalue = numeric(0))
    
    expr = fread(paste0(E, "/", E, "_scaled_quant_rank_TPM.txt"), data.table = FALSE)
    expr = expr[,-4]
    colnames(expr)[1:3] = c("marker","chrom","pos")
    
    ## Read in results for K model and make a scatterplot for each significant gene
    k_q_model = fread(paste0(E, "/", E, "_", P, "_rrBLUP_TWAS_K", ".txt"))
    
    sig_genes = filter(k_q_model, p_value_log10 > 2) #即p_value<0.01
    
    genes <- sig_genes$geneid #得到geneid
    
    for (gene in genes) {
      # make a scatterplot for each significant gene
      sig_gene_expr = expr[expr$marker==gene,]
      sig_gene_expr = data.frame(t(sig_gene_expr[,-c(2:3)]))
      colnames(sig_gene_expr) = sig_gene_expr[1,1]
      sig_gene_expr$Taxa = rownames(sig_gene_expr)
      sig_gene_expr = sig_gene_expr[-1,]
      plot_data = left_join(phenos, sig_gene_expr)
      plot_data[,2:14] = as.data.frame(lapply(plot_data[,2:14], as.numeric))
      colnames(plot_data)[14] <- "gene"
      
      ## Plot it
      p <- ggplot(plot_data, aes(gene, get(P))) +
        geom_point(alpha = 0.5, fill = "plum3")+
        geom_smooth(method = 'lm',col="darkred", aes(group=1), se=FALSE, formula = y ~ x)+
        theme_test(base_size = 15)+
        ylab(P)+
        xlab(paste0("Scaled expression phenotype of ", E)) +
        ggtitle(colnames(sig_gene_expr)[1]) +
        theme(plot.title = element_text(hjust = 0.5))
      
      model.lm <- lm(gene ~ get(P), data = plot_data)
      summary(model.lm)
      
      l <- list(a = format(coef(model.lm)[1], digits = 4),
                b = format(abs(coef(model.lm)[2]), digits = 4),
                r2 = format(summary(model.lm)$r.squared, digits = 4),
                p = format(summary(model.lm)$coefficients[2,4], digits = 4))
      
      eq <- substitute(~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
      
      p1 <- p + geom_text(aes(x = (max(gene)+min(gene))/2, y = max(get(P)+(max(get(P))-min(get(P)))/5), label = as.character(as.expression(eq))), parse = TRUE)
      
      ggsave(filename = paste0(E, "/", E, "_", P, "_", gene, "_TWAS_signif_scatterplot.pdf"), plot = p1, width = 8, height = 5)
      
      # 获得统计参数
      cor <- cor.test(plot_data$gene, plot_data[,P],alternative = "two.sided",method = "pearson",conf.level = 0.95)
      genes_cor[gene,] <- c(gene, E, P, cor$estimate["cor"], cor$p.value)
      
      write.table(x = genes_cor, file = paste0(E, "/", E, "_", P, "_TWAS_signif_list.txt"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
    }
  }
}



######################
# rbind 2E*3P个TWAS_signif_list.txt
######################
vigor <- data.frame(geneid = character(0), expr = character(0), pheno = character(0), cor = numeric(0), pvalue = numeric(0))

for (E in EE[3:4]) {
  pheno_names <- TWAS_list[TWAS_list$V1==E,2]
  for (P in pheno_names) {
    tmp <- read.table(file = paste0(E, "/", E, "_", P, "_TWAS_signif_list.txt"), header = TRUE, sep = "\t", row.names = NULL)
    vigor <- rbind(vigor, tmp)
  }
}

# save results
vigor_genes <- unique(vigor$geneid)
vigor_genes <- sort(vigor_genes)
write.table(x = vigor_genes, file = "vigor_genes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
save(vigor, file = "vigor.RData")
save(vigor_genes, file = "vigor_genes.RData")


#############################
# cor heatmap plot
#############################
## ----------------------------------------------------
library(tidyverse)
library(ComplexHeatmap) 
library(circlize)

# global plot settings by complexheatmap
ht_global_opt(
  legend_border = "black",
  heatmap_border = TRUE,
  annotation_border = TRUE)

load("/data2/usr/LiuWW/project_3/part9.twas/TWAS_TL/vigor.RData")
load("/data2/usr/LiuWW/project_3/part9.twas/TWAS_TL/vigor_genes.RData")

# prepare plot data
## cor data
long_data <- vigor %>% mutate(expr = str_sub(string = expr, start = 16, end = 25)) %>%
  mutate(expr_geneid = paste0(expr, "_", geneid)) %>%
  dplyr::select(expr_geneid, pheno, cor)
    
wide_data <- long_data %>% spread(pheno, cor)

head(wide_data)

plot_data <- wide_data %>% mutate(expr = str_sub(string = expr_geneid, start = 1, end = 10)) %>%
  mutate(geneid = str_sub(string = expr_geneid, start = -9, end = -1))


plot_data1 <- data.frame(geneid = vigor_genes,
                       expr = "TL_F1_vs_C") %>%
  mutate(expr_geneid = paste0(expr, "_", geneid)) %>% 
  left_join(y=plot_data, by = c("expr_geneid" = "expr_geneid")) %>%
  dplyr::select(-expr.y, -geneid.y)
colnames(plot_data1)[1:2] <- c("geneid", "expr")


plot_data2 <- data.frame(geneid = vigor_genes,
                         expr = "TL_F1_vs_P") %>%
  mutate(expr_geneid = paste0(expr, "_", geneid)) %>% 
  left_join(y=plot_data, by = c("expr_geneid" = "expr_geneid")) %>%
  dplyr::select(-expr.y, -geneid.y)
colnames(plot_data2)[1:2] <- c("geneid", "expr")


# df转为矩阵
rownames(plot_data1) <- plot_data1$geneid
rownames(plot_data2) <- plot_data2$geneid

m1 <- as.matrix(plot_data1[,c(4,6,8)])
m2 <- as.matrix(plot_data2[,c(5,7,9)])

# 填充缺失值
m1[is.na(m1)] <- 0
m2[is.na(m2)] <- 0


## column_info
### make
column_info <- data.frame(pheno = c("LA_F1vsC", "CN_F1vsC", "CA_F1vsC",
                                    "LA_F1vsP", "CN_F1vsP", "CA_F1vsP"),
                          heterosis = c("LA", "CN", "CA",
                                        "LA", "CN", "CA"),
                          expr = c("TL_F1_vs_C", "TL_F1_vs_C", "TL_F1_vs_C",  
                                   "TL_F1_vs_P", "TL_F1_vs_P", "TL_F1_vs_P"))

rownames(column_info) <- column_info$pheno

myorder <- c("LA_F1vsC", "LA_F1vsP",
             "CN_F1vsC", "CN_F1vsP",
             "CA_F1vsC", "CA_F1vsP")

column_info <- column_info[myorder,]

### 转为因子
column_info$pheno <- factor(column_info$pheno, levels = myorder)
column_info$heterosis <- factor(column_info$heterosis, levels = c("LA", "CN", "CA"))
column_info$expr <- factor(column_info$expr, levels = c("TL_F1_vs_C", "TL_F1_vs_P"))

m <- cbind(m1, m2)

m <- m[,myorder]

# 顶部注释
column_ha_t <- HeatmapAnnotation(trait = column_info$heterosis,
                               expr2pheno = column_info$expr,
                               col = list(trait = c("LA" = "#E64B35FF",
                                                        "CN" = "#00A087FF",
                                                        "CA" = "#3C5488FF"),
                                          expr2pheno = c("TL_F1_vs_C" = "#FDE74C",
                                                   "TL_F1_vs_P" = "#5BC0EB")
                               ),
                               gp = gpar(col = "black"),
                               annotation_name_side = "right")



## 统计基因数目
m_pos_info <- as.data.frame(colSums(m > 0))
m_neg_info <- as.data.frame(colSums(m < 0))

colnames(m_pos_info) <- "positive_count"
colnames(m_neg_info) <- "negative_count"

## 汇总
column_info <- cbind(column_info[,1:3], m_pos_info, m_neg_info)

# 底部注释
column_ha_b <- HeatmapAnnotation(positive_number = anno_barplot(column_info$positive_count,
                                                                height = unit(0.6, "cm"),
                                                                axis_param = list(
                                                                  at = c(0,80,160),
                                                                  labels = c("0","80","160"),
                                                                  gp = gpar(fontsize = 4)
                                                                )),
                                 negative_number = anno_barplot(column_info$negative_count,
                                                                height = unit(0.6, "cm"),
                                                                axis_param = list(
                                                                  at = c(0,80,160),
                                                                  labels = c("0","80","160"),
                                                                  gp = gpar(fontsize = 4)
                                                                )))
# general plot
# sort
m <- as.data.frame(m)
m_sort <- dplyr::arrange(m, m[,1], m[,2], m[,3], m[,4], m[,5], m[,6])

h <- Heatmap(m_sort,
             show_column_names = F,
             show_row_names = F,
             cluster_columns = F,
             cluster_rows = F,
             # 图例加名字
             name = "Correlation level",
             # 加列标题
             column_title = "Correlation level",
             column_title_side = "top",
             # 修改颜色
             col = colorRamp2(
               # 将 0，10，20
               breaks = c(-0.6, 0, 0.6),
               # 分别对应给三种颜色
               colors = c('blue','white', 'red')
             ),
             # 对列分割
             column_split = column_info$heterosis,
             #column_km = 1, 
             show_parent_dend_line = FALSE,
             # 外边框
             border = 'black',
             # 组间 gap
             column_gap = unit(0.01, 'npc'),
             # 在上方添加列注释
             top_annotation = column_ha_t,
             # 在下方添加列注释
             bottom_annotation = column_ha_b
             )
  
pdf(file = "Fig1_TWAS_TL_vigor_genes_cor_heatmap.pdf", width = 6, height = 8)  
draw(h, heatmap_legend_side = "right")
dev.off()

##################
# 候选基因韦恩图
##################
library(UpSetR)

## 读入每一个集合，并转化为list
vdata <- m

LA_F1vsC <- as.list(rownames(vdata[vdata$LA_F1vsC != 0,])) 
CN_F1vsC <- as.list(rownames(vdata[vdata$CN_F1vsC != 0,]))   
CA_F1vsC <- as.list(rownames(vdata[vdata$CA_F1vsC != 0,])) 

LA_F1vsP <- as.list(rownames(vdata[vdata$LA_F1vsP != 0,]))
CN_F1vsP <- as.list(rownames(vdata[vdata$CN_F1vsP != 0,]))
CA_F1vsP <- as.list(rownames(vdata[vdata$CA_F1vsP != 0,]))

## 用 Upset 自带的函数转换数据结构，符合画图要求
example <-list(LA_F1vsC = LA_F1vsC, 
               CN_F1vsC = CN_F1vsC, 
               CA_F1vsC = CA_F1vsC, 
               LA_F1vsP = LA_F1vsP,
               CN_F1vsP = CN_F1vsP,
               CA_F1vsP = CA_F1vsP) # 合并列表   
example <- fromList(example) # Upset 自带函数转化数据结构   

## 画图
setsBarColors <- c('#4DBBD5FF', '#91D1C2FF', '#F39B7FFF', '#3C5488FF', '#00A087FF', '#E64B35FF') # 设置集合颜色

pdf(file = "Fig3cd_TWAS_TL_vigor_genes_venn.pdf", width = 8, height = 6) 
upset(example,  
      nsets = length(example),  
      nintersects = 1000,  
      sets = c('CA_F1vsP','CN_F1vsP',"LA_F1vsP","CA_F1vsC","CN_F1vsC","LA_F1vsC"), # 这里的集合顺序和图上呈现的是相反的   
      keep.order = TRUE,  
      point.size = 4,  
      line.size = 1,  
      number.angles = 0,  
      text.scale = c(1, 1, 1, 1, 1, 1),  
      order.by = c("freq","degree"),
      decreasing = c(TRUE, FALSE),
      matrix.color ="#4285F4",  
      main.bar.color = 'black',  
      mainbar.y.label = 'candidate gene',
      mb.ratio = c(0.55, 0.45), # 控制上方条形图以及下方点图的比例
      sets.bar.color = setsBarColors)
dev.off()


#--------------------------------------------------------------------------------------
##########################
# Functional analysis
##########################
library(tidyverse)

## input tair10 gene function description
tair10 <- read.table(file = "../../TAIR10_functional_descriptions", header = TRUE, sep = "\t", quote = "", fill = TRUE)

tair10$Model_name <- str_sub(tair10$Model_name, start = 1, end = 9)

desc <- tair10[!duplicated(tair10$Model_name),]

rownames(desc) <- desc$Model_name

## get result
gene_func <- desc[vigor_genes,]

## output
write.table(x = gene_func, file = "Table1_TWAS_TL_vigor_genes_description.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#--------------------------------------------------------------------------------------
##########################
# GO analysis
##########################
## use clusterProfiler for GO analysis
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.At.tair.db))

##########################
# 全部vigor genes一起GO分析
##########################
ego <- enrichGO(gene = vigor_genes,
                #universe = vigor_genes,
                OrgDb = "org.At.tair.db",
                keyType="TAIR",
                ont="BP",
                pAdjustMethod = "fdr",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05
)
ego <- simplify(ego)
ego_result <- ego@result

## Save the result for future use
save(ego_result, ego_result, file = "ego_result.RData")
load("/data2/usr/LiuWW/project_3/part9.twas/TWAS_TL/ego_result.RData")

## plot
emapplot(ego, 
         # showCategory = 30,
         color = "p.adjust",
         layout = "nicely",
         by = count,
         includeAll = TRUE,
         line_scale = 0.6,
         pie_scale = 2)
## export the plots to "Fig2_TWAS_TL_vigor_genes_ego_emapplot"

# match GO term and genes
gene2term <- ego_result %>% dplyr::select(ID, Description, geneID, p.adjust)

## the first
for (term in gene2term$Description[1]) {
  gene <- strsplit(x = gene2term[gene2term$Description == term, "geneID"], split = "/")
  gene2term_long <- data.frame(Term = term, GeneID = unlist(gene), p.adjust = gene2term[gene2term$Description == term, "p.adjust"])
}

## the others
for (term in gene2term$Description[2:length(gene2term$Description)]) {
  gene <- strsplit(x = gene2term[gene2term$Description == term, "geneID"], split = "/")
  tmp <- data.frame(Term = term, GeneID = unlist(gene), p.adjust = gene2term[gene2term$Description == term, "p.adjust"])
  gene2term_long <- rbind(gene2term_long, tmp)
}

## Save the result for future use
save(gene2term_long, gene2term_long, file = "gene2term_long.RData")
load("/data2/usr/LiuWW/project_3/part9.twas/TWAS_TL/gene2term_long.RData")




##########################
# plot only known genes heatmap
##########################
#------------------ make row_anatation of GO_info
# long to wide
gene2term_wide <- spread(gene2term_long, key = "Term", value = "p.adjust")
GO_info <- column_to_rownames(gene2term_wide, var = "GeneID")

#------------------ get known_gene
known_gene <- unique(gene2term_wide$GeneID)

#------------------ get and sort mm data containing only known genes
mm <- m[rownames(m) %in% known_gene,]

mm <- dplyr::arrange(mm, mm[,1], mm[,2], mm[,3], mm[,4], mm[,5], mm[,6])

#------------------ plot all GO term
hh <- Heatmap(mm,
              show_column_names = F,
              show_row_names = F,
              cluster_columns = F,
              cluster_rows = F,
              # 图例加名字
              name = "Correlation level",
              # 加列标题
              column_title = "Correlation level",
              column_title_side = "top",
              # 修改颜色
              col = colorRamp2(
                # 将 0，10，20
                breaks = c(-0.6, 0, 0.6),
                # 分别对应给三种颜色
                colors = c('blue','white', 'red')
              ),
              # 对列分割
              column_split = column_info$heterosis,
              # 在上方添加列注释
              top_annotation = column_ha_t,
              # 在下方添加列注释
              #bottom_annotation = column_ha_b2,
              width = unit(8, "cm"),
              height = unit(12, "cm")
) +
  Heatmap(GO_info[rownames(mm),],
          show_row_names = F,
          cluster_columns = F,
          cluster_rows = F,
          # 图例加名字
          name = "GO term enrichment",
          # 加列标题
          column_title = "GO term for each gene",
          column_title_side = "top",
          # 列名角度
          column_names_rot = 45,
          # 设置颜色
          col = colorRamp2(c(0.05, 0.01),
                           c("#C5C3C6", "#0EAD69")
          ),
          na_col = "white",
          width = unit(16, "cm"),
          height = unit(12, "cm"))

pdf(file = "Fig3_TWAS_TL_known_vigor_genes_cor_heatmap.pdf", width = 15, height = 15)
draw(hh, heatmap_legend_side = "right")
dev.off()

#------------------ plot selected GO term
GO_info_select <- GO_info[rownames(mm), c(8,5,6,14,  4,11,12,13)]

# 分割列的信息
column_info2 <- data.frame(term = colnames(GO_info_select), class = c(rep("A",4), rep("B",4)))

column_info2 <- column_to_rownames(column_info2, var = "term")


#############定义函数删除全为NA的行##########
removeRowsAllNa <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
#############################################

#------------------ 删除GO_info_select中selected GO term全为NA的行
GO_info_select <- GO_info_select[rownames(removeRowsAllNa(GO_info_select)),]

#------------------ 删除mm中selected GO term全为NA的行
mmm <- mm[rownames(removeRowsAllNa(GO_info_select)),]

# 分割行的信息
row_info <- data.frame(gene = rownames(mmm), class = c(rep("A",12), rep("B",8), rep("C",25), rep("D",13), rep("E",3), rep("F",45)))

row_info <- column_to_rownames(row_info, var = "gene")



# 基因ID转变 TAIR to SYMBOL
library(clusterProfiler)
library(org.At.tair.db)

gene <- rownames(mmm)

gene.df <- bitr(gene, 
                fromType = "TAIR", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("SYMBOL"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.At.tair.db) #Orgdb是指对应的注释包是哪个
head(gene.df)

# 去除重复SYMBOL
gene.df <- gene.df %>% distinct(TAIR, .keep_all = TRUE) 

# get mark 作为行注释
x <- rownames_to_column(.data = mmm, var = "TAIR")

mark <- left_join(x, gene.df, by = "TAIR") %>%
  dplyr::select(TAIR, SYMBOL)
mark <- column_to_rownames(.data = mark, var = "TAIR")
mark$SYMBOL <- paste0(rownames(mark), " (", mark$SYMBOL, ")")

row_ha <- rowAnnotation(
  key_genes = anno_mark(which(rownames(mmm) %in% rownames(mmm)),
                        labels = mark$SYMBOL,
                        labels_gp = gpar(fontsize = 6))
)



hhs <- Heatmap(mmm,
              show_column_names = F,
              show_row_names = F,
              cluster_columns = F,
              cluster_rows = F,
              # 图例加名字
              name = "Correlation level",
              # 加列标题
              column_title = "Correlation level",
              column_title_side = "top",
              # 修改颜色
              col = colorRamp2(
                # 将 0，10，20
                breaks = c(-0.6, 0, 0.6),
                # 分别对应给三种颜色
                colors = c('blue','white', 'red')
              ),
              # 对列分割
              column_split = column_info$heterosis,
              # 对行分割
              row_split = row_info$class,
              # 在上方添加列注释
              top_annotation = column_ha_t,
              # 在右侧添加行注释
              right_annotation = row_ha,
              # cell边框颜色
              rect_gp = gpar(col = "white", lwd = 0.02),
              width = unit(10, "cm"),
              height = unit(40, "cm")
) +
  Heatmap(GO_info_select,
          show_row_names = F,
          cluster_columns = F,
          cluster_rows = F,
          # 图例加名字
          name = "GO term enrichment (p.adjust)",
          # 加列标题
          column_title = "GO term for each gene",
          column_title_side = "top",
          # 列名角度
          column_names_rot = 45,
          # 分割列
          column_split = column_info2$class,
          # 设置颜色
          col = colorRamp2(c(0.05, 0.01),
                           c("#C5C3C6", "#0EAD69")
          ),
          na_col = "white",
          # cell边框颜色
          rect_gp = gpar(col = "white", lwd = 0.02),
          width = unit(10, "cm"),
          height = unit(40, "cm"))

pdf(file = "Fig3_TWAS_TL_known_vigor_genes_cor_heatmap_selected_GOterm.pdf", width = 15, height = 26)
draw(hhs, heatmap_legend_side = "right")
dev.off()







######################################
# 以vigor_genes为背景进行GO分析
######################################

ID_list <- list()

classes <- c("LA_pos_F1vsC", "LA_pos_F1vsP", "CN_pos_F1vsC", "CN_pos_F1vsP", "CA_pos_F1vsC", "CA_pos_F1vsP",
             "LA_neg_F1vsC", "LA_neg_F1vsP", "CN_neg_F1vsC", "CN_neg_F1vsP", "CA_neg_F1vsC", "CA_neg_F1vsP")

for (class in classes){
  
  col <- paste0(str_sub(class,1,3), str_sub(class,-5,-1))
  
  if (grepl("pos",class) == TRUE) {
    ID_list[[class]] <- rownames(m[m[,col]>0,])
  } else {
    ID_list[[class]] <- rownames(m[m[,col]<0,])
  }
}



# use clusterProfiler for GO analysis
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.At.tair.db))

# enricher
GO_list <- list()
for (class in classes){
  GENE_ID <- ID_list[[class]]
  ego <- enrichGO(gene = GENE_ID, 
                  OrgDb = "org.At.tair.db",
                  keyType="TAIR",
                  ont="BP",
                  pAdjustMethod = "fdr",
                  universe = vigor_genes
  )
  ego <- simplify(ego)
  GO_list[[class]] <- ego@result
}

# Save the result for future use
save(ID_list, GO_list, file = "ID_GO_list.RData")
load(file = "ID_GO_list.RData")

# filter by p.adjust
for (class in classes) {
  GO_list[[class]] <- GO_list[[class]][which(GO_list[[class]][,6] < 0.05),]
}
# extract Term name and p.adjust
for (class in classes) {
  GO_list[[class]] <- GO_list[[class]][,c(2,6)]
  colnames(GO_list[[class]]) <- c("Term", class)
}
# merge
temp <- GO_list[[1]]
for (i in names(GO_list)[-1]) {
  df2 <- GO_list[[i]]
  temp <- merge(x = temp, y = df2, by = "Term", all.x = TRUE, all.y = TRUE)  
}

save(temp, file = "GO_merge.RData")
write.table(temp, file = "GO_merge.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# plot heat map
temp[,-1][!is.na(temp[,-1])] <- -log10(temp[,-1][!is.na(temp[,-1])])
##sort
library(tidyverse)
temp <- dplyr::arrange(temp, desc(LA_pos_F1vsC), desc(LA_pos_F1vsP), desc(CN_pos_F1vsC), desc(CN_pos_F1vsP), desc(CA_pos_F1vsC), desc(CA_pos_F1vsP), desc(LA_neg_F1vsC), desc(LA_neg_F1vsP), desc(CN_neg_F1vsC), desc(CN_neg_F1vsP), desc(CA_neg_F1vsC), desc(CA_neg_F1vsP))
##plot
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(require(circlize)) 
suppressPackageStartupMessages(library(dendextend))
mat <- as.matrix(temp[, -1])
rownames(mat) <-  temp[,"Term"]

##save plot
# 顶部注释
column_info_go <- rbind(column_info, column_info)[,1:3]
column_info_go$class <- str_sub(classes, 4,6)

column_ha_t <- HeatmapAnnotation(trait = column_info_go$heterosis,
                                 expr2pheno = column_info_go$expr,
                                 class = column_info_go$class,
                                 col = list(trait = c("LA" = "#E64B35FF",
                                                      "CN" = "#00A087FF",
                                                      "CA" = "#3C5488FF"),
                                            expr2pheno = c("TL_F1_vs_C" = "#FDE74C",
                                                           "TL_F1_vs_P" = "#5BC0EB"),
                                            class = c("pos" = "red",
                                                      "neg" = "blue")
                                 ),
                                 gp = gpar(col = "black"),
                                 annotation_name_side = "right")


pdf(file = "Fig2G_class_GO.pdf", width = 6, height = 12)
h <- Heatmap(mat, 
             col = colorRamp2(c(1.3, 5), c("yellow", "red")), 
             name = "-log10(p.adjust)", 
             na_col = "grey", 
             cluster_rows = FALSE, cluster_columns = FALSE, 
             row_title = "GO Terms",  
             column_names_side = "top",
             row_title_gp = gpar(fontsize = 10), 
             column_title_gp = gpar(fontsize = 3), 
             row_names_gp = gpar(fontsize =3), 
             column_names_gp = gpar(fontsize = 3), 
             column_title_side= "bottom", 
             heatmap_legend_param = list(title= "-log10(p.adjust)", 
                                         title_position ="leftcenter-rot", 
                                         legend_height=unit(2,"cm"), 
                                         legend_direction="vertical"), 
             row_names_max_width = unit(100, "mm"), 
             heatmap_width = unit(8, "cm"), 
             heatmap_height = unit(12, "cm"), 
             top_annotation = column_ha_t,
             # cell边框颜色
             rect_gp = gpar(col = "white", lwd = 0.01))
draw(h, heatmap_legend_side = "left")
dev.off()


