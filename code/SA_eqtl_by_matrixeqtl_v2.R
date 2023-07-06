setwd("/data2/usr/LiuWW/project_3/part8.eqtl")

##################
# reference: Getting started with Matrix eQTL  http://bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
# Matrix eQTL by Andrey A. Shabalin
# paper: Shabalin, A.A. Matrix eQTL: Ultra fast eQTL analysis via large matrix operations. Bioinformatics 28, no. 10 (2012): 1353-1358.
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Be sure to use an up to date version of R and Matrix eQTL.
##################

#--------------------------------------------------------------------------
# The toy dataset consists of five files: genotype SNP.txt, expression GE.txt, a file Covariates.txt with two covariates, gender and age, and files geneloc.txt and snpsloc.txt with gene and SNP location information. 

# First three files contain artificial information for 16 samples. The sample genotype file contains random values 0, 1, and 2 although the input genotype file is allowed to contain any real values.
#--------------------------------------------------------------------------

# 创建循环对象
SA_list <- dir(path="/data2/usr/LiuWW/project_3/part8.eqtl/data/GE_SA")

SA_list <- str_sub(SA_list[c(5)], start = 1, end = -5)

# 对循环对象进行循环
for (i in SA_list) {
  ############################################################
  # Test all gene-SNP pairs and plot a histogram of all p-values
  ############################################################
  
  #############
  # section 1
  #############
  
  # source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
  # load the package
  library(MatrixEQTL)
  
  ## Location of the package with the data files.
  # base.dir = find.package('MatrixEQTL')
  base.dir = '.'
  
  ## Settings [# set the parameters such as selected linear model and names of genotype and expression data files.]
  
  # Linear model to use, modelANOVA, or modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Genotype file name
  SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="")
  
  # Gene expression file name
  expression_file_name = paste(base.dir, "/data/GE_SA/", i, ".txt", sep="");
  
  # Covariates file name [A separate file may be provided with extra covariates.]
  # Set to character() for no covariates [In case of no covariates set the variable covariates_file_name to character()]
  covariates_file_name = character()
  # covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="")
  
  # Output file name
  output_file_name = tempfile()
  
  outdir = paste("./TWAS_candidate_eQTL/GE_SA/", i, sep = "")
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  # Only associations significant at this level will be saved [# The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset may cause excessively large output files.]
  # pvOutputThreshold = 1e-2;
  pvOutputThreshold = 1e-6;
  
  # Error covariance matrix [# Finally, define the covariance matrix for the error term. This parameter is rarely used. If the covariance matrix is a multiple of identity, set it to numeric().]
  # Set to numeric() for identity.
  # errorCovariance = numeric()
  # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
  
  
  
  #############
  # section 2
  #############
  # This section contains three very similar parts loading the files with genotype, gene expression, and covariates. In each part one can set the file delimiter (i.e. tabulation "\t", comma ",", or space " "), the string representation for missing values, the number of rows with column labels, and the number of columns with row labels. Finally, one can change the number of the variables in a slice for the file reading procedure (do not change if not sure).
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  
  # cvrt = SlicedData$new();
  # cvrt$fileDelimiter = "\t";      # the TAB character
  # cvrt$fileOmitCharacters = "NA"; # denote missing values;
  # cvrt$fileSkipRows = 1;          # one row of column labels
  # cvrt$fileSkipColumns = 1;       # one column of row labels
  # if(length(covariates_file_name)>0) {
  #   cvrt$LoadFile(covariates_file_name);
  # }
  
  
  #############
  # section 3: Run the analysis using the main Matrix eQTL function
  #############
  # Each significant gene-SNP association is recorded in a separate line in the output file and in the returned object me. In case of cis/trans eQTL analysis, two output files are produced, one with cis-eQTLs, another only with trans. Every record contains a SNP name, a transcript name, estimate of the effect size, t- or F-statistic, p-value, and FDR.
  
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    #cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    #errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  unlink(output_file_name);
  
  ## Results:
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
  cat('Detected eQTLs:', '\n');
  show(me$all$eqtls)
  ## save main results
  save(me, file = paste0(outdir, "/me.RData"))
  
  ###########################
  # Fig0 Plot the histogram of all p-values
  ###########################
  pdf(file = paste0(outdir, "/Fig0.Plot_the_histogram_of_all_p-values.pdf"), width = 12, height = 5)
  plot(me)
  dev.off()
  
  
  ############################################################
  # Test local and distant gene-SNP pairs separately and plot Q-Q plots of local and distant p-values [Cis- and trans- eQTL analysis]
  ############################################################
  # Matrix eQTL can distinguish local (cis-) and distant (trans-) eQTLs and perform separate correction for multiple comparisons for those groups. The main Matrix eQTL function Matrix_eQTL_main requires several extra parameters for cis/trans analysis:
  
  # pvOutputThreshold.cis – p-value threshold for cis-eQTLs.
  # output_file_name.cis – detected cis-eQTLs are saved in this file.
  # cisDist – maximum distance at which gene-SNP pair is considered local.
  # snpspos – data frame with information about SNP locations, must have 3 columns - SNP name, chromosome, and position. See sample SNP location file.
  # genepos – data frame with information about gene locations, must have 4 columns - the name, chromosome, and positions of the left and right ends. See sample gene location file.
  
  # source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
  # library(MatrixEQTL)
  
  ## Location of the package with the data files.
  # base.dir = find.package('MatrixEQTL');
  # base.dir = '.';
  
  ## Settings
  
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  # useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Genotype file name
  # SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
  snps_location_file_name = paste(base.dir, "/data/snpsloc.txt", sep="");
  
  # Gene expression file name
  # expression_file_name = paste(base.dir, "/data/GE.txt", sep="");
  gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="");
  
  # Covariates file name
  # Set to character() for no covariates
  # covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");
  
  # Output file name
  output_file_name_cis = tempfile();
  output_file_name_tra = tempfile();
  
  # Only associations significant at this level will be saved
  # pvOutputThreshold_cis = 2e-2;
  # pvOutputThreshold_tra = 1e-2;
  pvOutputThreshold_cis = 1e-6;
  pvOutputThreshold_tra = 1e-6;
  
  # Error covariance matrix
  # Set to numeric() for identity.
  # errorCovariance = numeric();
  # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
  
  # Distance for local gene-SNP pairs
  # cisDist = 1e6;
  cisDist = 1e6;
  
  ## Load genotype data
  
  # snps = SlicedData$new();
  # snps$fileDelimiter = "\t";      # the TAB character
  # snps$fileOmitCharacters = "NA"; # denote missing values;
  # snps$fileSkipRows = 1;          # one row of column labels
  # snps$fileSkipColumns = 1;       # one column of row labels
  # snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  # snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  # gene = SlicedData$new();
  # gene$fileDelimiter = "\t";      # the TAB character
  # gene$fileOmitCharacters = "NA"; # denote missing values;
  # gene$fileSkipRows = 1;          # one row of column labels
  # gene$fileSkipColumns = 1;       # one column of row labels
  # gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  # gene$LoadFile(expression_file_name);
  
  ## Load covariates
  
  # cvrt = SlicedData$new();
  # cvrt$fileDelimiter = "\t";      # the TAB character
  # cvrt$fileOmitCharacters = "NA"; # denote missing values;
  # cvrt$fileSkipRows = 1;          # one row of column labels
  # cvrt$fileSkipColumns = 1;       # one column of row labels
  # if(length(covariates_file_name)>0) {
  #   cvrt$LoadFile(covariates_file_name);
  # }
  
  
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  
  
  ###########################
  # get lead snp
  ###########################
  df <- me[["all"]][["eqtls"]]
  
  # The significant SNPs for each trait (gene) were grouped into clusters with a maximum distance of 10 kb between two consecutive SNPs, and only those clusters with more than three significant SNPs were retained as putative eQTL.
  
  df <- df %>% 
    separate(snps, into = c("snp_chr", "snp_loc"),sep = ":",remove = FALSE) %>% 
    dplyr::arrange(gene, snp_chr, snp_loc)
  
  # 将df转为list,按gene、snp_chr进行分组
  list <- split(df, df$gene)
  
  for (i in c(1:length(list))) {
    list[[i]] <- split(list[[i]], list[[i]]$snp_chr)
  }
  
  
  # 对每个分组进行分cluster
  for (j in c(1:length(list))) {
    for (k in c(1:length(list[[j]]))) {
      #模块1
      group <- list[[j]][[k]]
      group$cluster <- "no"
      group$cluster[1] <- "yes"
      
      #模块2
      if (nrow(group) > 2) {
        #计算距离
        for (m in c(2:nrow(group))) {
          n <- m-1
          if (as.numeric(group$snp_loc[m])-as.numeric(group$snp_loc[n])>=10001) {
            group$cluster[m] <- "yes"
          }
        }
        
        #重编行名
        rownames(group) <- c(1:nrow(group))
        #每个cluster的首行行号
        co <- rownames(group[group$cluster == "yes",])
        #cluster编号
        coc <- paste0("cluster", 1:length(co))
        #置换首行的yes标记为cluster行号
        for (p in c(1:length(co))) {
          group[co[p],"cluster"] <- coc[p]
        }
        
        #每个cluster首行的随后行也随着标为cluster号
        for (b in c(2:nrow(group))) {
          d <- b-1
          if (group$cluster[b] == "no") {
            group$cluster[b] <- group$cluster[d]
          }
        }
        
        list[[j]][[k]] <- group
      } else {
        group$cluster <- "less_snp"
        list[[j]][[k]] <- group
      }
    }
  }
  
  # 将list转为df
  ## 转第一层
  list2 <- list()
  for (s in c(1:length(list))) {
    list2[[s]] <- list[[s]][[1]]
    if (length(list[[s]]) >= 2) {
      for (x in c(2:length(list[[s]]))) {
        list2[[s]] <- rbind(list2[[s]], list[[s]][[x]])
      }
    }
  }
  
  ## 转第二层
  df2 <- list2[[1]]
  for (x in c(2:length(list2))) {
    df2 <- rbind(df2, list2[[x]])
  }
  
  # 删除cluster内snp数目<=3个的cluster
  df2$flag1 <- paste(df2$gene, df2$snp_chr, df2$cluster, sep=",")
  
  ## 计数
  flag1 <- as.data.frame(table(df2$flag1))
  
  filter <- flag1[flag1$Freq <= 3,]
  
  ## filter
  df3 <- df2 %>% filter(!(flag1 %in% filter$Var1))
  
  
  # 计算每个cluster中P值最小的snp为该eQTL的lead snp
  lead_p <- aggregate(df3$pvalue, by=list(type=df3$flag1), min)
  
  ## 找出所有可作为lead的snp，其中部分snp的p值相同小，可稍后再任选一个
  df4 <- NULL
  for (g in 1:nrow(lead_p)) {
    tmp <- df3[df3$pvalue == lead_p$x[g] & df3$flag1 == lead_p$type[g],]
    df4 <- rbind(df4, tmp)
  }
  
  ## snp的p值相同小的选第一个,得到最终的lead SNP代表eQTL
  lead_snp <- df4[!duplicated(df4$flag1),]
  
  ## save main results
  save(lead_snp, file = paste0(outdir, "/lead_snp.RData"))
  
  
  ###########################
  # Fig1 散点图
  ###########################
  # 从结果中提取所需要的列
  library(tidyverse)
  
  data <- lead_snp
  
  data <- data %>% dplyr::select(c("snps", "gene", "pvalue")) 
  
  # 获得绘散点图数据
  # part1
  data <- left_join(x = data, y = genepos, by = c("gene" = "geneid"))
  
  data <- data %>% dplyr::select(-c("s2"))
  
  colnames(data) <- c("snps", "gene", "pvalue", "chr_gene", "left_gene")
  
  # part2
  data <- left_join(x = data, y = snpspos, by = c("snps" = "snp"))
  
  colnames(data) <- c("snps", "gene", "pvalue", "chr_gene", "left_gene", "chr_snp", "pos_snp")
  
  # 得到绘制散点图原始数据
  pdata <- data %>% dplyr::select(c("chr_gene", "left_gene", "chr_snp", "pos_snp", "pvalue"))
  
  entire <- data.frame(chr_gene=c("chr1","chr2","chr3","chr4","chr5", "chr1","chr2","chr3","chr4","chr5"),
                       left_gene=c(0,0,0,0,0, 30427671, 19698289, 23459830, 18585056, 26975502),
                       chr_snp=c("chr1","chr2","chr3","chr4","chr5", "chr1","chr2","chr3","chr4","chr5"),
                       pos_snp=c(0,0,0,0,0, 30427671, 19698289, 23459830, 18585056, 26975502),
                       pvalue=c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
  pdata <- rbind(pdata, entire)
  
  pdata$chr_gene <- as.factor(pdata$chr_gene)
  pdata$chr_snp <- as.factor(pdata$chr_snp)
  
  # 绘制散点图
  library(ggplot2)
  library(ggsci)
  library(cowplot)
  
  ## 图片主体
  sp <- ggplot(pdata, aes(x=pos_snp, y=left_gene)) + 
    geom_point(shape=20, size=1, aes(color=pvalue)) +
    labs(x="snp", y="gene") +
    theme_test() +
    scale_color_gsea()
  
  ## 分面【对chr_gene进行垂直分割, 对chr_snp进行水平分割】
  ### 自定义分面标签顺序
  pdata <- within(pdata, chr_snp <- factor(chr_snp, levels = c("chr1", "chr2", "chr3", "chr4", "chr5")))
  with(pdata, levels(chr_snp))
  
  pdata <- within(pdata, chr_gene <- factor(chr_gene, levels = c("chr5", "chr4", "chr3", "chr2", "chr1")))
  with(pdata, levels(chr_gene))
  
  ### 分面
  sp <- sp + facet_grid(chr_gene ~ chr_snp, scales="free", space="free", switch = "both") +
    theme(#调整panel标签的字体
      strip.text.x = element_text(size=8),
      strip.text.y = element_text(size=8),
      #调整panel标签的颜色和底纹
      strip.background = element_rect(colour="#000000", fill="#E4FDE1"),
      #调整panel的间距
      panel.spacing = unit(0,"cm"),
      #调整panel和panel标签之间的间距
      strip.placement = "inside",
      strip.switch.pad.grid = unit(0, "cm"),
      #调整panel边框颜色
      panel.border = element_rect(color = "black", size=0.2, linetype="dashed")
    )
  
  ## 保存
  ggsave(sp, file=paste0(outdir, "/Fig1.sp.pdf"), width=7, height=6) # cm
  
  
  ## Run the analysis
  me2 = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    #cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    #errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  unlink(output_file_name_tra);
  unlink(output_file_name_cis);
  
  ## Results:
  cat('Analysis done in: ', me2$time.in.sec, ' seconds', '\n');
  cat('Detected local eQTLs:', '\n');
  show(me2$cis$eqtls)
  cat('Detected distant eQTLs:', '\n');
  show(me2$trans$eqtls)
  
  ## save main results
  save(me2, file = paste0(outdir, "/me2.RData"))
  
  ## Plot the Q-Q plot of local and distant p-values
  pdf(file = paste0(outdir, "/Fig2.Q-Q_plot_of_local_and_distant_p-values.pdf"), width = 12, height = 5)
  plot(me2)
  dev.off()
}

