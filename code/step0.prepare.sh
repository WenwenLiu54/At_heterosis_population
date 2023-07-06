#!/bin/bash
# 准备基础数据：对输入数据做格式转换、文件处理

#--------------------------------------------------------------------------
# step0.prepare.sh 
# Usage: nohup sh step0.prepare.sh > step0.prepare.log 2>&1 &
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2021-5-17
#--------------------------------------------------------------------------

# 输入文件为之前做过 hard filter + 去除染色体Mt和Pt只剩下核染色体12345 + only biallelic + missing和maf过滤 的vcf文件
# 链接该文件到当前目录 
ln -s /data2/usr/LiuWW/project_2/part3.polulation_genetics/00.filter/all.missing_maf.vcf

# or

# 输入文件为之前做过 hard filter + 去除染色体Mt和Pt只剩下核染色体12345 + only biallelic + missing和maf过滤 + LD过滤的vcf文件
# 链接该文件到当前目录 
# ln -s /data2/usr/LiuWW/project_2/part3.polulation_genetics/00.filter/all.LDfilter.vcf

#------------------------------------------------
# variables

vcf="./all.missing_maf.vcf"

# or

# vcf="./all.LDfilter.vcf"

#------------------------------------------------


## First, before missing and maf filter, use vcftools to include only 96 samples

~/miniconda3/bin/vcftools --vcf ${vcf} \
    --keep /data2/usr/LiuWW/project_2/part4.GWAS/06.GWAS_GEMMA_96pop/sample_list_96.txt \
    --recode \
    --recode-INFO-all \
    --out all.96pop
    
# --recode：重新编码为vcf文件
# --recode-INFO-all： This option can be used with the recode option. This option is used to keep all INFO values in the original file.


## filter missing and maf（速度快）
plink --vcf all.96pop.recode.vcf `# 输入的vcf文件` \
    --geno 0.05 `# missing过滤` \
    --maf 0.05 `# maf过滤，保留多态性高的位点` \
    --out all.96pop.missing_maf `# 输出文件前缀` \
    --recode vcf-iid  `# 让输出的结果也是vcf格式（因为plink默认会将文件转成ped或bed格式）` \
    --allow-extra-chr `# 建议加，允许其他染色体编号，防止报错` \
    --set-missing-var-ids @:# `# 建议加，给每个snp赋一个id名称，后面的GWAS和后续某些过滤步骤等也需要用snp id，原本的vcf文件这一列是用.占位的（在plink里，@表示染色体，#表示snp位置）` \
    --keep-allele-order # plink软件默认读入的vcf/ped后，基于等位基因频率的高低，将ref和alt调换，因此必须加这一参数
    


## filter LD （目的：把连锁的snp过滤掉，structure分析要求是不相关的、不连锁的位点，必须过滤LD；PCA和进化树构建过滤或不过滤LD都可以）注意！！！！！！！该步千万不要重复多次，做一次就可以，否则位点数目一直减少。

### 进行LD过滤，生成需要保留的位点文件：tmp.ld.prune.in
plink --vcf  all.96pop.missing_maf.vcf \
    --indep-pairwise 50 10 0.2 `# 通过算R^2值进行过滤，三个数字分别是：以50个snp为一个窗口，去统计窗口内部的R^2值，通过算法保证最后两两snp没有超过0.2的R^2值的snp对在，10是步长————step，指算完一个窗口后再移动10个snp计算下一个窗口，注意：50是指50个snp，如果不想以50个snp做过滤，想以物理距离做过滤，只需在50后面写kb，意思是一个窗口大小是50kb【这里的50 10 0.2就是最常用的阈值，同时要考虑剩下的snp的个数进行调整】` \
    --out tmp.ld  `# 输出文件前缀，该步生成两个文件————prune.in（要保留的snp位点）和prune.out（要丢弃的snp位点），文件内记录的是snp id信息` \
    --allow-extra-chr \
    --set-missing-var-ids @:# 

### 对vcf文件进行过滤，生成LD过滤后的vcf文件
plink --vcf  all.96pop.missing_maf.vcf  \
    --make-bed `# 生成一份bed格式（snp的一种存储格式）的文件，后续会使用到` \
    --extract tmp.ld.prune.in  `# 指定要保留的snp位点` \
    --out all.96pop.LDfilter `# 输出文件前缀` \
    --recode vcf-iid  \
    --keep-allele-order  \
    --allow-extra-chr \
    --set-missing-var-ids @:#


#-----------------------------------------------------------------------------------------------
# prepare vcf（将基因型数据vcf转换为gemma软件需要的输入格式，即plink二进制文件格式。这里使用plink软件进行。plink是一种使用最广泛的关联分析软件，其定义的ped/map文件系统，及其对应的bed/bim/fam已经成为关联分析的标准文件格式。在进行关联分析之前，首先要做的就是将其他格式的文件转换为plink对应的文件格式。）


# PLINK软件输入文件的常见格式类型：
# 1，一般格式：PED/MAP
# 2，转置格式：TPED/TFAM
# 3，二进制格式：BED/BIM/FAM
# 几种格式之间可以相互转换。推荐使用BED/BIM/FAM这种格式，读取速度快。
# bed文件包含SNP数据，是二进制格式，不能由Notepad++等文本编辑器打开。
# bim文件包括SNP位置信息。
# fam文件包括家系表型信息，这两种文件都是文本格式。


## GEMMA分析时候需要plink二进制文件：.bed .bim .fam（.fam文件有坑，要注意）。其中bed存储基因型信息，bim存储每个遗传变异（通常是SNP）的相关信息，fam存储样本信息和顺序。
## plink二进制文件分为三部分:.bed 包含基因分型的二进制文件；.fam 包含家庭号、个体号、父亲号、母亲号、性别、表型(即.ped文件前6列)；.bim 标记信息文件，类似于.map文件。


## 第1步：将vcf转成ped，该步生成的结果包括.log .map .ped .nosex，其中map和ped是主要结果

plink --vcf all.96pop.missing_maf.vcf `# 输入的vcf文件` \
    --recode `# plink1.9版本支持直接读取vcf、gen等多种文件格式，并默认会将不同的格式转换为二进制的bed文件格式，添加--recode参数将输出结果调整为ped格式` \
    --double-id `# 防plink默认的下划线分割导致出错，通过此参数使得family id和smple id保持相同，而双亲、性别默认是用0填充，表型默认是用-9填充` \
    --out snp_file `# 输出文件前缀` \
    --set-missing-var-ids @:# `# 缺失位点没有id的snp加一下snp` \
    --allow-extra-chr # 允许其他的染色体格式
    

## 第2步：将ped/map转换为fam/bed/bim，即将plink格式转化为二进制的plink格式

plink --file snp_file `# map, ped文件的前缀`\
    --make-bed \
    --out snp_bfile `# 生成的bed, bim, fam这三个文件的前缀` \
    --set-missing-var-ids @:# `# 缺失位点没有id的snp加一下snp` \
    --allow-extra-chr # 允许其他的染色体格式


#-----------------------------------------------------------------------------------------------
# prepare phenotype（对表型做样本排序，要求与基因型snp文件样品顺序一致）

## 前面还没有输入表型信息，.ped文件和.fam文件的第6列都是默认值-9。注意要按照样本顺序替换这一列的表型值。

## 这里表型数据有两种方法：（1）可以放在.fam文件第六列，用 -n 1表示。也可以放在第七列，用 -n 2表示。两者结果是完全一样的。（2）将表型数据提取为单独一列，同样按照样本顺序。

## 这里采用第（2）种方法

## 第1步：将ped/map转换为tped/tfam格式 【tfam（记录样品顺序）和tped（记录基因型信息）】

plink --file snp_file --recode --transpose --out snp_tfile \
    --set-missing-var-ids @:# `# 缺失位点没有id的snp加一下snp` \
    --allow-extra-chr # 允许其他的染色体格式


## 第2步：每种表型准备为一个文件，每次单独对一个表型做GWAS
for trait in `cat /data2/usr/LiuWW/project_2/part0.phenotype/trait_for_gwas/trait_list.txt`
do
    if [ ! -d ${trait} ];
        then mkdir ${trait}
    fi

    ### perl脚本 上一步的tfam文件 表型数据（两列：样品名 表型值） >输出的表型数据（三列：样品名称出现两次，第三列为表型值，并且样品顺序和tfam文件样品顺序一致）
    perl ../script/sort_pheno.pl snp_tfile.tfam /data2/usr/LiuWW/project_2/part0.phenotype/trait_for_gwas_pop96/${trait} > ./${trait}/${trait}.sort.txt

    awk '{print $3}' ./${trait}/${trait}.sort.txt > ./${trait}/${trait}.p.txt
done
#-----------------------------------------------------------------------------------------------
