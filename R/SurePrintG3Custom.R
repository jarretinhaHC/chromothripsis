# Remember to create apropriate conf files
# source('rconf/agilent.rconf')

library(limma)
library(ArrayTV)
library(doMC)
library(foreach)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(gtools)
library(DNAcopy)

# Confs
registerDoMC(cores=24)
genome <- BSgenome.Hsapiens.UCSC.hg19
basedir <- '/data/Cytogenomics/Agilent'
outdir <- '/data/Cytogenomics/Agilent/results'
datadir <- paste(basedir, 'raw/afe', sep='/')
target_sheet <- 'targets.tsv'
targets <- readTargets(paste(basedir, target_sheet, sep='/'), row.names='Name') 

# HELPER FUNCTIONS

chr_names_to_numbers <- function(chr){

    tmp <- gsub('chr', '', chr)
    tmp <- gsub('X', '23', tmp)
    tmp <- gsub('Y', '24', tmp)
    tmp <- gsub('MT', '25', tmp)
    return(as.integer(tmp))

}

multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
    do.call(order, c(lapply(list(...),

                            function(l){
                                if(is.character(l)){

                                    factor(l, levels=mixedsort(unique(l)))

                                } else {

                                    l

                                }
                            }), list(na.last = na.last, decreasing=decreasing)))

}

### READ DATA

# Read data, annotate design, printer, names, etc.
raw <- read.maimages(targets, path=datadir, source='agilent', names=targets$Name)
raw$design <- rep(-1, dim(raw)[[2]])

# Printer layout isn't right, need additional Block column 
# Info from marray
# maNgr 1
# maNgc 1
# maNsr 1064
# maNsc 170
# maNspots 180880
# maSub 180875
raw$printer <- data.frame(ngrid.r=1,
                          ngrid.c=1,
                          nspot.r=dim(raw)[[1]],
                          nspot.c=1)

raw$ID <- targets$Name

# Quality control
limma::plotMA3by2(raw, prefix='QualityControlMA - Raw', path=outdir)

# Chromosome info and position are buried inside genes$SystematicName
# Need to split the string and retrieve chr and start coord
# Only from genes$ControlType == 0 (1 and -1 are true controls)

chr_set <- paste0('chr', c(1:22, 'X', 'Y'))
coords <- ifelse(!raw$genes$ControlType, strsplit(raw$genes$SystematicName, '[:-]'),
                 NA)

raw$genes$ID <- 1:dim(raw$genes)[[1]]
raw$genes$chr <- sapply(coords, function(x) ifelse(x[1] %in% chr_set, x[1], NA))
raw$genes$start <- sapply(coords, function(x) ifelse(x[1] %in% chr_set,
                                                        as.integer(x[2]), NA))
raw$genes$end <- sapply(coords, function(x) ifelse(x[1] %in% chr_set,
                                                        as.integer(x[3]), NA))

### CLEANING, NORMALIZATION & CORRECTIONS

# Auto is fine for background correction
# Method 'printtiploess' performed best
tmp <- backgroundCorrect(raw)
tmp <- normalizeWithinArrays(tmp)

# Get rid of controls and empty spots
# A bit tricky but works
tmp <- tmp[which(tmp$genes$chr %in% chr_set),]
tmp <- tmp[apply(tmp$M, 1, function(x) !all(is.na(x))), ]

results <- list()
# Parameter for ArrayTV
# Aparently, the bigger the window/increment, higher the TV score
maxwins <- c(1e7)
increms <- c(1e6)

for(sample in tmp$ID){

    # Genomic waves correction
    # ArrayTV works well
    # Of course, don't know how to deal with NAs
    # So, need to work on a sample basis
    NA_flag <- !is.na(tmp$M[, sample])
    val <- as.matrix(tmp$M[NA_flag, sample])
    chr <- tmp$genes$chr[NA_flag]
    start <- tmp$genes$start[NA_flag]
    end <- tmp$genes$end[NA_flag]
    probenames <- tmp$genes$ProbeName[NA_flag]

    result <- gcCorrect(object=val,
                        chr=chr,
                        starts=start,
                        increms=increms,
                        maxwins=maxwins,
                        build='hg19')$correctedVals
    cna <- smooth.CNA(CNA(result,
                          chr_names_to_numbers(chr),
                          start,
                          sampleid=sample))

    segs <- segment(cna,
                    alpha=0.01,
                    p.method='perm',
                    undo.splits='sdundo',
                    undo.SD=1.96,
                    min.width=5)

    pdf(paste(outdir, paste0(sample, '.cbs.pdf'), sep='/'))
    plot(segs, plot.type='chrombysample', ylim=c(-2, 2), xmaploc=TRUE)
    dev.off()

    segs <- segs$output
    results[[sample]] <- as(data.frame(chr=segs$chrom,
                                       start=segs$loc.start,
                                       end=segs$loc.end,
                                       nprobes=segs$num.mark,
                                       seg.mean=segs$seg.mean), 'GRanges')



}


