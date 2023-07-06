#!/bin/bash

# Usage:nohup sh step1.hisat_build.sh > step1.hisat_build.log 2>&1 &

hisat2-build ./genome.fasta ./genome 
