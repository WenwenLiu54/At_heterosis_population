#!/bin/bash
#set several traps
set -e
set -u
set -o pipefail

#------------------------------------------------------------------------------------

# step1.hisat2_mapping.sh
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2020-10-4
# Usage: sh step1.hisat2_mapping.sh 3-C-1 /data2/usr/LiuWW/project_3/part1.quality_control/step1.quality_control_result

#------------------------------------------------------------------------------------

##parameters required to change (based to different machines' limitation)

#------------------------------------------------------------------------------------

##Locations require to change (based to different computers) && Assign the I/O variables

sample_name=${1}    # 3-C-1
raw_location=${2}   # /data2/usr/LiuWW/project_3/part1.quality_control/step1.quality_control_result

SAMTOOLS="/data2/usr/LiuWW/software/samtools/samtools"

output="/data2/usr/LiuWW/project_3/part2.hisat2_mapping/step1.hisat2_mapping_result/${sample_name}" ## Each sample should have a uniq Dir to store files

#------------------------------------------------------------------------------------

##mkdir new directories

if [ ! -d ${output} ]
   then mkdir -p ${output}
fi

#------------------------------------------------------------------------------------

##mapping 

echo "mapping with hisat2!" > ${output}/${sample_name}_general.log

hisat2 --new-summary -t -p 10 --score-min L,0,-0.4 -x /data2/usr/LiuWW/project_3/part0.ref/genome -1 ${raw_location}/${sample_name}/clean.${sample_name}_R1.fq.gz -2 ${raw_location}/${sample_name}/clean.${sample_name}_R2.fq.gz -S ${output}/${sample_name}.sam 1>${sample_name}.log 2>&1

#p threads
#x reference genome
#S output

#------------------------------------------------------------------------------------

##sam to bam and sort

echo "sam to bam and sort!" >> ${output}/${sample_name}_general.log

${SAMTOOLS} sort -@ 8 -o ${output}/${sample_name}.bam ${output}/${sample_name}.sam

rm ${output}/${sample_name}.sam 

##get unique and bam to sam

echo "get unique and bam to sam!" >> ${output}/${sample_name}_general.log

${SAMTOOLS} view -h ${output}/${sample_name}.bam | grep -E "^@|NH:i:1$|NH:i:1[^0-9]" > ${output}/${sample_name}.unique.sam

rm ${output}/${sample_name}.bam 

##sam to bam

echo "sam to bam!" >> ${output}/${sample_name}_general.log

${SAMTOOLS} view -h -bS ${output}/${sample_name}.unique.sam > ${output}/${sample_name}.unique.bam

rm ${output}/${sample_name}.unique.sam 


##bamindex
echo "bamindex!" >> ${output}/${sample_name}_general.log
${SAMTOOLS} index ${output}/${sample_name}.unique.bam

