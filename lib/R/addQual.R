#!/usr/bin/Rscript
args <- commandArgs(T)
qcFile <- args[1]
output <- args[2]
funcR <- args[3]
lociFiltFile <- args[4]

read.delim(lociFiltFile, header=F) -> lociFiltD
lociFilt <- paste(paste("chr", lociFiltD[,1], sep=""), lociFiltD[,2], sep=":")

source(funcR)
read.delim(qcFile, header=T) -> qc

#addQual(qc) -> qc.qual
addQualTransCor(qc, lociFilt) -> qc.qual
#addQualVarCor(qc) -> qc.qual

write.table(qc.qual, output, row.names=F, col.names=T, quote=F, sep="\t")
#write.table(qc[selection,], output, row.names=F, col.names=T, quote=F, sep="\t")


