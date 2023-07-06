#!/bin/bash
#set several traps
set -e
set -u
set -o pipefail

#------------------------------------------------------------------------------------

# step01.BQSR.sh
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2020-09-23
# Usage: sh step01.BQSR.sh Per-1 /data2/usr/LiuWW/project_2/part1.variant_calling/03.SNP_indel/step0.mapq20_filter_result

#------------------------------------------------------------------------------------

##parameters required to change (based to different machines' limitation)

#------------------------------------------------------------------------------------

##Locations require to change (based to different computers) && Assign the I/O variables

sample_name=${1}    # Per-1

raw_location=${2}   # /data2/usr/LiuWW/project_2/part1.variant_calling/03.SNP_indel/step0.mapq20_filter_result


output="/data2/usr/LiuWW/project_2/part1.variant_calling/03.SNP_indel/step01.BQSR_result/${sample_name}" ## Each sample should have a uniq Dir to store files

#------------------------------------------------------------------------------------

##mkdir new directories

if [ ! -d ${output} ]
    then mkdir -p ${output}
fi

#------------------------------------------------------------------------------------

##step01.BQSR

echo "step01.BQSR begin" > ${output}/${sample_name}.log

#BaseRecalibrator碱基矫正(BQSR)包含了两个步骤：
##第一步，BaseRecalibrator，这里计算出了所有需要进行重校正的read和特征值，然后把这些信息输出为一份校准表文件。在计算的过程中，不考虑已知的变异位点的碱基质量，--known-sites指定已知变异位点对应的vcf文件。这一步对单个样本进行操作，每个样本生成一个错误模型文件。
##第二步，PrintReads，这一步利用第一步得到的校准表文件重新调整原来BAM文件中的碱基质量值，并使用这个新的质量值重新输出一份新的BAM文件。BQSR会对输入的bam文件中的碱基质量值进行替换，替换为校正之后的质量值，而原先的质量值保存在OQtag中。
##【注意，因为BQSR实际上是为了（尽可能）校正测序过程中的系统性错误，因此，在执行的时候是按照不同的测序lane或者测序文库来进行的，这个时候@RG信息（BWA比对时所设置的）就显得很重要了，算法就是通过@RG中的ID来识别各个独立的测序过程，这也是我开始强调其重要性的原因。】
gatk --java-options "-Xmx20g -Djava.io.tmpdir=./tmp"  BaseRecalibrator \
    -R ../01.ref/genome.fasta  \
    -I ${raw_location}/${sample_name}/${sample_name}.sort.rmdup.filter.bam  \
    --known-sites ../data/known1001vcf/1001genomes_snp-short-indel_only_ACGTN.vcf \
    -O ${output}/${sample_name}.sort.rmdup.filter.recal.table \
    1>${output}/${sample_name}.sort.rmdup.filter.recal.log 2>&1 

gatk --java-options "-Xmx20g -Djava.io.tmpdir=./tmp"   ApplyBQSR \
    -R ../01.ref/genome.fasta  \
    -I ${raw_location}/${sample_name}/${sample_name}.sort.rmdup.filter.bam  \
    -bqsr ${output}/${sample_name}.sort.rmdup.filter.recal.table \
    -O ${output}/${sample_name}.sort.rmdup.filter.bqsr.bam \
    1>${output}/${sample_name}.sort.rmdup.filter.bqsr.log 2>&1 

echo "step01.BQSR done!!!" >> ${output}/${sample_name}.log
#------------------------------------------------------------------------------------

