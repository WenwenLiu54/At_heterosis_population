#!/bin/bash

# 进行进化树构建、PCA和structure分析（都只用snp，不用其它形式变异，因为snp在基因组上的密度最高）前，会对SNP进行LD过滤，保留质量较高、多态性较高的位点。做GWAS不需进行LD过滤，只需missing和MAF过滤。

# 这里使用plink软件，vcftools软件也可以做。过滤的阈值主要看过滤后的snp数目是否够后续的分析使用，一般几百k个位点就足够了（一般做完missing和maf过滤会剩下几兆位点，再做完LD过滤后就剩下几百k了）


#--------------------------------------------------------------------------
# biallelic_missing_maf_LD_vcf.sh 
# Usage: nohup sh biallelic_missing_maf_LD_vcf.sh > biallelic_missing_maf_LD_vcf.log 2>&1 &
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2021-4-15
#--------------------------------------------------------------------------


# 输入文件为之前做过hard filter的vcf文件，并去除了染色体Mt和Pt只剩下核染色体12345
vcf=/data2/usr/LiuWW/project_2/part1.variant_calling/03.SNP_indel/step6.SNP_nuclear_only_by_vcftools_filter_result/all.filtered.snp.nuclear_only.recode.vcf



## First, before missing and maf filter, use vcftools to include only bi-allelic sites

~/miniconda3/bin/vcftools --vcf ${vcf} \
    --min-alleles 2 \
    --max-alleles 2 \
    --recode \
    --recode-INFO-all \
    --out all.biallelic
    
# --recode：重新编码为vcf文件
# --recode-INFO-all： This option can be used with the recode option. This option is used to keep all INFO values in the original file.



## filter missing and maf（速度快）
plink --vcf all.biallelic.recode.vcf `# 输入的vcf文件` \
    --geno 0.05 `# missing过滤` \
    --maf 0.05 `# maf过滤，保留多态性高的位点` \
    --out all.missing_maf `# 输出文件前缀` \
    --recode vcf-iid  `# 让输出的结果也是vcf格式（因为plink默认会将文件转成ped或bed格式）` \
    --allow-extra-chr `# 建议加，允许其他染色体编号，防止报错` \
    --set-missing-var-ids @:# `# 建议加，给每个snp赋一个id名称，后面的GWAS和后续某些过滤步骤等也需要用snp id，原本的vcf文件这一列是用.占位的（在plink里，@表示染色体，#表示snp位置）` \
    --keep-allele-order # plink软件默认读入的vcf/ped后，基于等位基因频率的高低，将ref和alt调换，因此必须加这一参数
    


## filter LD （目的：把连锁的snp过滤掉，structure分析要求是不相关的、不连锁的位点，必须过滤LD；PCA和进化树构建过滤或不过滤LD都可以）注意！！！！！！！该步千万不要重复多次，做一次就可以，否则位点数目一直减少。

### 进行LD过滤，生成需要保留的位点文件：tmp.ld.prune.in
plink --vcf  all.missing_maf.vcf  \
    --indep-pairwise 50 10 0.2 `# 通过算R^2值进行过滤，三个数字分别是：以50个snp为一个窗口，去统计窗口内部的R^2值，通过算法保证最后两两snp没有超过0.2的R^2值的snp对在，10是步长————step，指算完一个窗口后再移动10个snp计算下一个窗口，注意：50是指50个snp，如果不想以50个snp做过滤，想以物理距离做过滤，只需在50后面写kb，意思是一个窗口大小是50kb【这里的50 10 0.2就是最常用的阈值，同时要考虑剩下的snp的个数进行调整】` \
    --out tmp.ld  `# 输出文件前缀，该步生成两个文件————prune.in（要保留的snp位点）和prune.out（要丢弃的snp位点），文件内记录的是snp id信息` \
    --allow-extra-chr \
    --set-missing-var-ids @:# 

### 对vcf文件进行过滤，生成LD过滤后的vcf文件
plink --vcf  all.missing_maf.vcf  \
    --make-bed `# 生成一份bed格式（snp的一种存储格式）的文件，后续会使用到` \
    --extract tmp.ld.prune.in  `# 指定要保留的snp位点` \
    --out all.LDfilter `# 输出文件前缀` \
    --recode vcf-iid  \
    --keep-allele-order  \
    --allow-extra-chr \
    --set-missing-var-ids @:#  

