#!/bin/bash

#----------------------------------------------------------------------------------------------------------------------
# step6.SNP_nuclear_only_by_vcftools_filter.sh 
# Usage: nohup sh step6.SNP_nuclear_only_by_vcftools_filter.sh > step6.SNP_nuclear_only_by_vcftools_filter.log 2>&1 &
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2021-4-12
#----------------------------------------------------------------------------------------------------------------------

# Introduction: use vcftools to filter chromosomes Mt Pt and only reserve nuclear chromosomes 1,2,3,4,5

mkdir step6.SNP_nuclear_only_by_vcftools_filter_result

~/miniconda3/bin/vcftools --vcf step4.SNP_result/all.filtered.snp.vcf \
    --chr 1 \
    --chr 2 \
    --chr 3 \
    --chr 4 \
    --chr 5 \
    --recode \
    --recode-INFO-all \
    --out step6.SNP_nuclear_only_by_vcftools_filter_result/all.filtered.snp.nuclear_only
    
# --recode：重新编码为vcf文件
# --recode-INFO-all： This option can be used with the recode option. This option is used to keep all INFO values in the original file.
