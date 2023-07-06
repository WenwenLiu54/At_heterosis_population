#!/bin/bash
#set several traps
set -e
set -u
set -o pipefail

#------------------------------------------------------------------------------------

# step0.mapq20_filter.sh
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2020-09-19
# Usage: sh step0.mapq20_filter.sh Per-1 /data2/usr/LiuWW/project_2/part1.variant_calling/02.mapping/step2.MarkDuplicates_result

#------------------------------------------------------------------------------------

##parameters required to change (based to different machines' limitation)

#------------------------------------------------------------------------------------

##Locations require to change (based to different computers) && Assign the I/O variables

sample_name=${1}    # Per-1

raw_location=${2}   # /data2/usr/LiuWW/project_2/part1.variant_calling/02.mapping/step2.MarkDuplicates_result


output="/data2/usr/LiuWW/project_2/part1.variant_calling/03.SNP_indel/step0.mapq20_filter_result/${sample_name}" ## Each sample should have a uniq Dir to store files

#------------------------------------------------------------------------------------

##mkdir new directories

if [ ! -d ${output} ]
    then mkdir -p ${output}
fi

#------------------------------------------------------------------------------------

##step0.mapq20_filter

echo "step0.mapq20_filter begin" > ${output}/${sample_name}.log
  

# get unique mapping reads and mapping quality > 20

/data2/usr/LiuWW/software/samtools/samtools view -h -@ 20 -q 20 -F 4 -F 256 ${raw_location}/${sample_name}/${sample_name}.sort.rmdup.bam | grep -v 'XA:Z:' | /data2/usr/LiuWW/software/samtools/samtools view -h -@ 20 -b -o  ${output}/${sample_name}.sort.rmdup.filter.bam

echo "step0.mapq20_filter done!!!" >> ${output}/${sample_name}.log
#------------------------------------------------------------------------------------

