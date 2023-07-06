#!/bin/bash
#set several traps
set -e
set -u
set -o pipefail

#------------------------------------------------------------------------------------

# step2.MarkDuplicates.sh
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2020-09-14
# Usage: sh step2.MarkDuplicates.sh Per-1 /data2/usr/LiuWW/project_2/part1.variant_calling/02.mapping/step1.bwa_sort_result

#------------------------------------------------------------------------------------

##parameters required to change (based to different machines' limitation)

#------------------------------------------------------------------------------------

##Locations require to change (based to different computers) && Assign the I/O variables

sample_name=${1}    # Per-1

raw_location=${2}   # /data2/usr/LiuWW/project_2/part1.variant_calling/02.mapping/step1.bwa_sort_result


output="/data2/usr/LiuWW/project_2/part1.variant_calling/02.mapping/step2.MarkDuplicates_result/${sample_name}" ## Each sample should have a uniq Dir to store files

#------------------------------------------------------------------------------------

##mkdir new directories

if [ ! -d ${output} ]
                then mkdir -p ${output}
fi

#------------------------------------------------------------------------------------

##step2.MarkDuplicates

echo "step2.MarkDuplicates begin" > ${output}/${sample_name}.log

picard -Xmx4g  MarkDuplicates I=${raw_location}/${sample_name}/${sample_name}.sort.bam O=${output}/${sample_name}.sort.rmdup.bam  CREATE_INDEX=true  REMOVE_DUPLICATES=true M=${output}/${sample_name}.marked_dup_metrics.txt

echo "step2.MarkDuplicates done!!!" >> ${output}/${sample_name}.log
#------------------------------------------------------------------------------------

