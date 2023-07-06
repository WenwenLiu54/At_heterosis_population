#!/bin/bash

# Usage: nohup sh step2.CombineGVCFs.sh > step2.CombineGVCFs.log 2>&1 &

ls step1.HaplotypeCaller_result/*/*.g.vcf > ./step2.CombineGVCFs_result/gvcf.list

gatk  --java-options "-Xmx20g -Djava.io.tmpdir=./tmp"    CombineGVCFs \
    -R ../01.ref/genome.fasta \
    -V ./step2.CombineGVCFs_result/gvcf.list  \
    -O ./step2.CombineGVCFs_result/all.merge.g.vcf 

