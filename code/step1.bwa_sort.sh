#!/bin/bash
#set several traps
set -e
set -u
set -o pipefail

#------------------------------------------------------------------------------------

# step1.bwa_sort.sh
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2020-09-14
# Usage: sh step1.bwa_sort.sh Per-1 /data2/usr/LiuWW/project_2/part1.variant_calling/data/reseq_data_ln1_quality_control

#------------------------------------------------------------------------------------

##parameters required to change (based to different machines' limitation)

#------------------------------------------------------------------------------------

##Locations require to change (based to different computers) && Assign the I/O variables

sample_name=${1}    # Per-1

raw_location=${2}   # /data2/usr/LiuWW/project_2/part1.variant_calling/data/reseq_data_ln1_quality_control



output="/data2/usr/LiuWW/project_2/part1.variant_calling/02.mapping/step1.bwa_sort_result/${sample_name}" ## Each sample should have a uniq Dir to store files

#------------------------------------------------------------------------------------

##mkdir new directories

if [ ! -d ${output} ]
                then mkdir -p ${output}
fi

#------------------------------------------------------------------------------------

##step1.bwa_sort

echo "step1.bwa_sort begin" > ${output}/${sample_name}.log

bwa mem -t 8 -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:illumina" ../01.ref/genome.fasta ${raw_location}/reseq_data_ln1_${sample_name}_quality_control_result/clean.${sample_name}_R1.fq.gz ${raw_location}/reseq_data_ln1_${sample_name}_quality_control_result/clean.${sample_name}_R2.fq.gz | /data2/usr/LiuWW/software/samtools/samtools sort -@ 8 -m 1G -o ${output}/${sample_name}.sort.bam -

echo "step1.bwa_sort done!!!" >> ${output}/${sample_name}.log
#------------------------------------------------------------------------------------

