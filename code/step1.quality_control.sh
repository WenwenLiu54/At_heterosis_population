#!/bin/bash
#set several traps
set -e
set -u
set -o pipefail

#------------------------------------------------------------------------------------

# step1.quality_control.sh
# Author: Wenwen Liu (liuwenwen@pku.edu.cn)
# Last modified: 2020-09-27
# Usage: sh step1.quality_control.sh *1.fq.gz *2.fq.gz Per-1 /data2/usr/LiuWW/project_3/data/Cleandata_ln

#------------------------------------------------------------------------------------

##parameters required to change (based to different machines' limitation)

threads="3"

#------------------------------------------------------------------------------------

##Locations require to change (based to different computers) && Assign the I/O variables

fastqfile1=${1}     # *1.fq.gz
fastqfile2=${2}     # *2.fq.gz
sample_name=${3}    # Per-1
raw_location=${4}   # /data2/usr/LiuWW/project_3/data/Cleandata_ln

QC_output="/data2/usr/LiuWW/project_3/part1.quality_control/step1.quality_control_result/${sample_name}"

log="${QC_output}/${sample_name}_log"
QC="${QC_output}/${sample_name}_quality_control"
QC_clean="${QC_output}/${sample_name}_quality_control_clean"

#------------------------------------------------------------------------------------

##mkdir new directories
if [ ! -d ${QC_output} ]
                then mkdir -p ${QC_output}
fi

if [ ! -d ${log} ]
                then mkdir -p ${log}
fi

if [ ! -d ${QC} ]
                then mkdir -p ${QC}
fi

if [ ! -d ${QC_clean} ]
                then mkdir -p ${QC_clean}
fi

#------------------------------------------------------------------------------------

##Quality control on raw/clean sequencing data from the company by fastqc and fastp

##Quality control on raw/clean sequencing data from the company by fastqc

echo "Calculate the basic fastq info!" >> ${QC_output}/${sample_name}_general_log.txt

seqkit stats ${raw_location}/${sample_name}/${fastqfile1} > ${log}/${sample_name}_1_seq_stats.txt 
seqkit stats ${raw_location}/${sample_name}/${fastqfile2} > ${log}/${sample_name}_2_seq_stats.txt 

echo "Fastqc perform quality control on clean fastq file!" >> ${QC_output}/${sample_name}_general_log.txt

fastqc -t ${threads} -o ${QC} ${raw_location}/${sample_name}/${fastqfile1} ${raw_location}/${sample_name}/${fastqfile2}  

##Quality control on raw/clean sequencing data from the company by fastp

echo "Auto fastp fastq file!" >> ${QC_output}/${sample_name}_general_log.txt

clean1fq="${QC_output}/clean.${sample_name}_R1.fq.gz"
clean2fq="${QC_output}/clean.${sample_name}_R2.fq.gz"

fastp -w ${threads} -i ${raw_location}/${sample_name}/${fastqfile1}  -o ${clean1fq} -I ${raw_location}/${sample_name}/${fastqfile2} -O ${clean2fq} -j "${QC_output}/${sample_name}.fastp.json" -h "${QC_output}/${sample_name}.fastp.html" 

##Quality control on after_fastp_clean sequencing data by a second fastqc

echo "Fastqc perform quality control on after_fastp_clean fastq file!" >> ${QC_output}/${sample_name}_general_log.txt

fastqc -t ${threads} -o ${QC_clean} ${clean1fq} ${clean2fq} 

echo "Calculate the after_fastp_clean fastq info!" >> ${QC_output}/${sample_name}_general_log.txt

seqkit stats ${clean1fq} > ${log}/${sample_name}_clean_1_seq_stats.txt  
seqkit stats ${clean2fq} > ${log}/${sample_name}_clean_2_seq_stats.txt  

#------------------------------------------------------------------------------------

