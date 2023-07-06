#!/usr/bin/Rscript
# parse parameter ---------------------------------------------------------
library(argparser, lib.loc = "/home/LiuWW/R/x86_64-redhat-linux-gnu-library/3.6", quietly=TRUE)
# Create a parser
p <- arg_parser("run featureCounts and calculate FPKM/TPM")

# Add command line arguments
p <- add_argument(p, "--bam", help="input: bam file", type="character")
p <- add_argument(p, "--gtf", help="input: gtf file", type="character")
p <- add_argument(p, "--output", help="output prefix", type="character")

# Parse the command line arguments
argv <- parse_args(p)

library(Rsubread, lib.loc = "/home/LiuWW/R/x86_64-redhat-linux-gnu-library/3.6")
library(limma, lib.loc = "/home/LiuWW/R/x86_64-redhat-linux-gnu-library/3.6")
library(edgeR, lib.loc = "/home/LiuWW/R/x86_64-redhat-linux-gnu-library/3.6")

bamFile <- argv$bam
gtfFile <- argv$gtf
nthreads <- 14
outFilePref <- argv$output

outStatsFilePath  <- paste(outFilePref, '.log',  sep = ''); 
outCountsFilePath <- paste(outFilePref, '.count', sep = ''); 

fCountsList = featureCounts(bamFile, annot.ext=gtfFile, isGTFAnnotationFile=TRUE, nthreads=nthreads, isPairedEnd=TRUE)
dgeList = DGEList(counts=fCountsList$counts, genes=fCountsList$annotation)
fpkm = rpkm(dgeList, dgeList$genes$Length)
tpm = exp(log(fpkm) - log(sum(fpkm)) + log(1e6))

write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

featureCounts = cbind(fCountsList$annotation[,1], fCountsList$counts, fpkm, tpm)
colnames(featureCounts) = c('gene_id', 'counts', 'fpkm','tpm')
write.table(featureCounts, outCountsFilePath, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

