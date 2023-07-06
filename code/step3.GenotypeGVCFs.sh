#!/bin/bash

# Usage: nohup sh step3.GenotypeGVCFs.sh > step3.GenotypeGVCFs.log 2>&1 &

gatk  --java-options "-Xmx20g -Djava.io.tmpdir=./tmp"   GenotypeGVCFs \
    -R ../01.ref/genome.fasta \
    --variant ./step2.CombineGVCFs_result/all.merge.g.vcf \
    -O ./step3.GenotypeGVCFs_result/all.merge_raw.vcf
