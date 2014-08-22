# load configuration stuff
# source('HumanCytoSNP12.conf')
# Load libs/packages

#Large file support
library(ff)
# Set packages options
options(ffcaching="ffeachflush")

library(humancytosnp12v2p1hCrlmm)
library(oligoClasses)
library(Biobase)
library(IRanges)
library(crlmm)
library(VanillaICE)
library(ArrayTV)
library(data.table)
library(doParallel)
library(foreach)

# Let's parallelize calculations.
registerDoParallel()
ocSamples(60)

# Set data folder and references
basedir <- '/data/Cytogenomics/HumanCyto12'
datadir <- paste(basedir, 'IDAT', sep='/')
sampleSheet <-read.delim(paste(basedir, 'SampleSheet/simples.tsv', sep='/'),
               as.is=TRUE)
outdir <- paste(basedir, 'results', getRversion(), sep='/')
outdir_large <- paste(basedir, 'results_large', getRversion(), sep='/')
dir.create(outdir, recursive=TRUE)

# Large file support
ldPath(outdir_large)
dir.create(outdir_large, recursive=TRUE)

# Define data, batches and related stuff
arrayNames <- file.path(datadir, sampleSheet$Sample_ID)
arrayInfo <- list(barcode=NULL, position='SentrixPosition')
cdfName <- 'humancytosnp12v2p1h'
batch <- rep('1', nrow(sampleSheet))
sns <- sampleSheet$Sample_name

# Check everything is in place
message('Everything green is OK?')
print(all(file.exists(paste(arrayNames, "_Grn.idat", sep=""))))
message('Everything red is OK?')
print(all(file.exists(paste(arrayNames, "_Red.idat", sep=""))))

# Additional checks
message('Hey dude! Did you double check everything? Calling data is costly!!!')

# Genotyping
# Fails when parallel mode in on, probably during some merge
# So, use one core a carry on

options(cores=1)
gen <- genotype.Illumina(sampleSheet=sampleSheet,
                         arrayNames=arrayNames,
                         arrayInfoCol=arrayInfo,
                         cdf=cdfName,
                         batch=batch,
                         sns=sns)
# Check data quality
snr <- cnSet$SNR[]
barcodes <- rapply(strsplit(sampleNames(cnSet), '_'), function(x) head(x, 1))
snr <- data.frame(SNR=snr, Array=barcodes)

pdf(paste('DataQualityHistogramOverall.pdf', sep='/'))
plt <- ggplot(data=snr, aes(SNR))
plt <- plt + geom_bar(binwidth=1, colour='black', alpha='0.75')
plt <- plt + ylab('Amostras')
print(plt)
dev.off()

pdf(paste(outdir, 'DataQualityHistogramArray.pdf', sep='/')
plt <- ggplot(data=snr, aes(SNR, fill=Array))
plt <- plt + geom_bar(binwidth=1, colour='black', alpha='0.75') + ylab('Amostras')
plt <- plt + facet_wrap(~Array)
print(plt)
dev.off()

# Copy number estimation
# Also fails when cores > 1. Damn!!!
crlmmCopynumber(gen)

# Downstream analysis
# Data conversion to BafLrrSetList (bls)
# Will generate some warnings
bls <- BafLrrSetList(gen)

# Illumina generates data outside normal chromosome range
# Not sure what it means. Should I ask on foruns?
# Functions will not work on this strange data
# Let's filter it out for now

cleaned_bls <- bls[-(1:1)] 

# Correct for genomic waves
# As this is a BafLrrSetList, return value is NULL
# Everything is done in place
# How do I calibrate the window size?

maxwins <- 10^(4:6)
increms <- 10^(3:5)
gcCorrect(cleaned_bls, maxwins=maxwins, increms=increms)

# Fit HMM to identify CNVs
# Takes a long time!!!
fit <- unlist(hmmBafLrrSetList2(cleaned_bls))

# Idea is to reduce ranges to a single data track and add a coverage track
for(st in 1:6){
    
    # Filter segments with same state
    ranges <- fit[state(fit) == st, ]

    # Convert to ranges with coverage
    ranges_with_coverage <- as(coverage(ranges), 'GRanges')

    # Reduce to get a idea of total length
    reduced_ranges <- reduce(ranges)

    # Annotate samples per segment
    overlaps <- findOverlaps(ranges, ranges_with_coverage)
    qidx <- queryHits(overlaps)
    sidx <- subjectHits(overlaps)
    tmp <- tapply(sampleNames(ranges)[qidx], sidx, list)
    ranges_with_coverage$samples <- NA
    ranges_with_coverage$samples[as.integer(names(tmp))] <- tmp
    ranges_with_coverage$samples <- sapply(ranges_with_coverage$samples,
                                           FUN=paste, collapse='|')

    # Print coverage data as tsv and to feed in Circos
    df <- as.data.frame(ranges_with_coverage, row.names=NULL)
    df <- as.data.table(df)
    df <- df[df$score > 0, ]

    # Table data/Coverage
    filename <- paste0('ranges_with_coverage_', st, '.tsv')
    write.table(df,
                file=paste(outdir, filename, sep='/'),
                sep='\t',
                row.names=FALSE,
                quote=FALSE)

    # Circos/Coverage
    data <- data.frame('CHR'=gsub('chr', 'hs', df$seqnames),
                       'START'=df$start,
                       'END'=df$end)

	data$SCORE <- paste0('score=', df$score)
    filename <- paste0('ranges_with_coverage_for_circos_', st, '.dat')
	write.table(data,
                file=paste(outdir, filename, sep='/'),
                sep='\t',
                row.names=FALSE,
                col.names=FALSE,
                quote=FALSE)

    # Circos/Reduced range
	df <- as.data.frame(reduced_ranges, row.names=NULL)
    df$strand <- c()
    df <- df[df$width == 1, ]

    data <- data.frame('CHR'=gsub('chr', 'hs', df$seqnames),
                       'START'=df$start,
                       'END'=df$end)
    filename <- paste0('reduced_ranges_for_circos_', st, '.dat')
    write.table(data,
                file=paste(outdir, filename, sep='/'),
                sep='\t',
                row.names=FALSE,
                col.names=FALSE,
                quote=FALSE)

}

