library(XLConnect)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

genome <- BSgenome.Hsapiens.UCSC.hg19
basedir <- '/home/jarretinha/Dropbox/Gil'
datafile <- 'Gil_MLPA.xls' 
wb <- loadWorkbook(paste(basedir, datafile, sep='/'))
results <- readWorksheet(wb, sheet='Results', check.names=FALSE)
probes <-  readWorksheet(wb, sheet='Probes', check.names=FALSE)


