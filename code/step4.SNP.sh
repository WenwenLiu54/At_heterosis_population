#!/bin/bash

# Usage: nohup sh step4.SNP.sh > step4.SNP.log 2>&1 &

# 把SNP单独提取出来
gatk  --java-options "-Xmx20g -Djava.io.tmpdir=./tmp" \
    SelectVariants \
    -R ../01.ref/genome.fasta \
    -V ./step3.GenotypeGVCFs_result/all.merge_raw.vcf \
    --select-type SNP \
    -O ./step4.SNP_result/all.raw.snp.vcf

# 加过滤标签：满足条件可以pass（接受），不满足的加标签‘SNP_filter’
# 关于过滤表达式————指定了在什么情况下会加标签（这里是gatk官方推荐的一组标准，||是或的关系）
# QD：位点质量值/位点深度，相当于平均质量值（因为位点覆盖度越高，gatk给的打分越高），越高越好
# MQ：read比对质量值的均方根，越高越好
# FS：考虑read比对正负链的偏向性问题，理论上应该是随机的，支持正负链的read数目应该大致差不多，差异很大时，可能是有错误
# SOR：也是过滤链特异性方面
# MQRankSum：也是过滤链特异性方面，是针对杂合位点时的过滤标准
# ReadPosRankSum：也是与杂合位点相关的过滤，若某位点为纯合则没有这个值
# --filter-name：满足任意一条过滤标准则加标签'SNP_filter'
# 在log文件里会有满屏的警告信息，没关系，无需担心，是因为有一些行没有一些统计值例如MQRankSum
gatk  --java-options "-Xmx20g -Djava.io.tmpdir=./tmp" \
    VariantFiltration \
    -R ../01.ref/genome.fasta \
    -V ./step4.SNP_result/all.raw.snp.vcf \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name 'SNP_filter' \
    -O ./step4.SNP_result/all.filter.snp.vcf

# 删除加了过滤标签的位点
gatk  --java-options "-Xmx20g -Djava.io.tmpdir=./tmp" \
    SelectVariants \
    -R ../01.ref/genome.fasta \
    -V ./step4.SNP_result/all.filter.snp.vcf \
    --exclude-filtered \
    -O ./step4.SNP_result/all.filtered.snp.vcf
