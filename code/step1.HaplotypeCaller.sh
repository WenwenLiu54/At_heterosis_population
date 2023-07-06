#!/bin/bash
#set several traps
set -e
set -u
set -o pipefail

#------------------------------------------------------------------------------------

# step1.HaplotypeCaller.sh
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2020-09-24
# Usage: sh step1.HaplotypeCaller.sh Per-1 /data2/usr/LiuWW/project_2/part1.variant_calling/03.SNP_indel/step01.BQSR_result

#------------------------------------------------------------------------------------

##parameters required to change (based to different machines' limitation)

#------------------------------------------------------------------------------------

##Locations require to change (based to different computers) && Assign the I/O variables

sample_name=${1}    # Per-1

raw_location=${2}   # /data2/usr/LiuWW/project_2/part1.variant_calling/03.SNP_indel/step01.BQSR_result


output="/data2/usr/LiuWW/project_2/part1.variant_calling/03.SNP_indel/step1.HaplotypeCaller_result/${sample_name}" ## Each sample should have a uniq Dir to store files

#------------------------------------------------------------------------------------

##mkdir new directories

if [ ! -d ${output} ]
    then mkdir -p ${output}
fi

#------------------------------------------------------------------------------------

##step1.HaplotypeCaller

echo "step1.HaplotypeCaller begin" > ${output}/${sample_name}.log

gatk --java-options "-Xmx20g -Djava.io.tmpdir=./tmp" HaplotypeCaller \
    -R ../01.ref/genome.fasta \
    -I ${raw_location}/${sample_name}/${sample_name}.sort.rmdup.filter.bqsr.bam \
    -ERC GVCF \
    -O ${output}/${sample_name}.g.vcf \
    1>${output}/${sample_name}.HC.log 2>&1

echo "step1.HaplotypeCaller done!!!" >> ${output}/${sample_name}.log
#------------------------------------------------------------------------------------

