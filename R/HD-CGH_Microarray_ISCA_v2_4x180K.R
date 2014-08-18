# Remember to create apropriate conf files
# source('rconf/agilent.rconf')

library(doMC)
library(foreach)
library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(CGHcall)
library(reshape2)

# Confs
registerDoMC(cores=24)
genome <- BSgenome.Hsapiens.UCSC.hg19
basedir <- '/data/Cytogenomics/Agilent'
datadir <- paste(basedir, 'raw/afe', sep='/')
target_sheet <- 'targets.tsv'
targets <- readTargets(paste(basedir, target_sheet, sep='/'), row.names='Name') 

# Read target files and add some info
# Agilent data is previously dye bias normalized with aCGH.Spline
# Why the writeFE options create so badly behaved filenames?
# batch.spline(dir=datadir, writeFE=TRUE)

raw <- read.maimages(targets, path=datadir, source='agilent', names=targets$Name)
raw$design <- rep(-1, dim(raw)[[2]])
raw$ID <- targets$Name

# Chromosome info and position are buried inside genes$SystematicName
# Need to split the string and retrieve chr and start coord
# Only from genes$ControlType == 0 (1 and -1 are true controls)

chr_set <- paste0('chr', c(1:22, 'X', 'Y'))
coords <- ifelse(!raw$genes$ControlType, strsplit(raw$genes$SystematicName, '[:-]'),
                 NA)

raw$genes$ID <- 1:dim(raw$genes)[[1]]
raw$genes$chr <- sapply(coords, function(x) ifelse(x[1] %in% chr_set, x[1], NA))
raw$genes$start <- sapply(coords, function(x) ifelse(x[1] %in% chr_set,
                                                        as.numeric(x[2]), NA))
raw$genes$end <- sapply(coords, function(x) ifelse(x[1] %in% chr_set,
                                                        as.numeric(x[3]), NA))

# Manual normalization based on Chen, J. et al. (2011)
# See doi: 10.1186/gb-2011-12-8-r80 

# We don't need controls and related stuff without coords
tmp <- raw[which(raw$genes$chr %in% chr_set),]

# Remove background based on median values, threshold as you like
threshold <- 5
median_bg_G <- apply(tmp$Gb, 2, median, na.rm=TRUE)
median_bg_R <- apply(tmp$Rb, 2, median, na.rm=TRUE)

bg_flag <- (tmp$G > threshold * median_bg_G) & (tmp$R > threshold * median_bg_R)
tmp$G[!bg_flag] <- NA
tmp$R[!bg_flag] <- NA

# Clean bg data
tmp$Gb <- tmp$Rb <- c()

# Remove probes without at least one valid intensity
# A bit tricky but works
tmp$genes <- tmp$genes[apply(tmp$G, 1, function(x) !all(is.na(x))), ]
tmp$genes <- tmp$genes[apply(tmp$R, 1, function(x) !all(is.na(x))), ]

tmp$G <- tmp$G[apply(tmp$G, 1, function(x) !all(is.na(x))), ]
tmp$R <- tmp$R[apply(tmp$R, 1, function(x) !all(is.na(x))), ]

# Log intensities

tmp$G <- log2(tmp$G)
tmp$R <- log2(tmp$R)

# Calculate trimmed means
# Trim parameter is calculated with aid of quantiles/scale

calculate_trim <- function(x){

    observed <- scale(quantile(x, na.rm=TRUE, prob=seq(0, 1, 0.01)))
    mn <- mean(x, na.rm=TRUE)
    sd <- sd(x, na.rm=TRUE)
    expected <- scale(quantile(rnorm(1000, mn, sd), prob=seq(0, 1, 0.01)))
    diff <- observed - expected
    trim <- 0.01 * length(diff[which(abs(diff) > 1)])
    return(trim)

}

trimG <- calculate_trim(tmp$G)
trimR <- calculate_trim(tmp$R)

# Trimmed means per sample
tmp$Gtm <- apply(tmp$G, 2, mean, trim=trimG, na.rm=TRUE)
tmp$Rtm <- apply(tmp$R, 2, mean, trim=trimR, na.rm=TRUE)

# SD per sample
tmp$Gsd <- apply(tmp$G, 2, sd, na.rm=TRUE)
tmp$Rsd <- apply(tmp$R, 2, sd, na.rm=TRUE)

for(sample in  tmp$ID){

    tmp$G[, sample] <- (tmp$G[, sample] - tmp$Gtm[sample]) / tmp$Gsd[sample] 
    tmp$R[, sample] <- (tmp$R[, sample] - tmp$Rtm[sample]) / tmp$Rsd[sample] 

}

# Normalized ratio and average intensities
tmp$M <- tmp$G - tmp$R
tmp$A <- (tmp$G + tmp$R) / 2

# Genomic waves correction
# Parameters suggested by previous studies with ArrayTV
# Window size of about 120kbp for Agilent

ranges <- sortSeqlevels(as(tmp$genes[, c('chr', 'start', 'end')], 'GRanges'))
seqlengths(ranges) <- seqlengths(genome)[which(seqlevels(genome) %in% chr_set)]

width <- 60e4
seqs <- getSeq(genome, ranges + width)
# GC content calculations can take a long time
# Let's parallelize it!
tmp$genes$GC <- foreach(i=1:length(seqs), .combine=c) %dopar%
letterFrequency(seqs[i], letters='GC', as.prob=TRUE)

# Set GC weights
tmp$genes$w <- ifelse(tmp$genes$GC < 0.5, 1e-3, 1e3)

# Fit model, just love R!
# Residuals are the wave corrected M
# We need to watch for NAs on sample basis

# Just using tmp$M as a template
# Yeah, I know it's not safe

tmp$Mwc <- tmp$M
for(sample in tmp$ID){

    fit <- lm(tmp$M[, sample] ~ tmp$genes$GC + I(tmp$genes$GC^2),
              weights=tmp$genes$w,
              na.action=na.exclude)
    tmp$Mwc[, sample] <- resid(fit)

}

# Now, correct artifacts
# Remove medians per sample
for(sample in tmp$ID){

    tmp$Mwc[, sample] <- tmp$Mwc[, sample] - median(tmp$Mwc[, sample], na.rm=TRUE)

}
# Probe means across samples for profiling
tmp$genes$ProbeMeanMwc <- apply(tmp$Mwc, 1, mean, na.rm=TRUE)

# Calculate profile 
tmp$genes$ProbeProfile <- smooth.spline(tmp$genes$ProbeMeanMwc, all.knots=TRUE)$y
# Get the final M value!!!
tmp$Mf <- tmp$Mwc - tmp$genes$ProbeProfile

# Segmentation and copy state calling with CGHCall
segs <- make_cghRaw(data.frame(tmp$genes$ProbeName, tmp$genes[, c('chr',
                                                                 'start',
                                                                 'end')],
                              tmp$Mf, check.names=FALSE))

# Data is already cleaned, so proceed to segmentation
# Default parameters are fine
segs <- segmentData(segs)
segs <- postsegnormalized(segs)

# It's show time! Call it!
result <- CGHcall(segs)
result <- ExpandCGHcall(result, segs)

# Construct data frame from calls and locations
df <- data.frame(calls(result), featureData(result)@data, check.names=FALSE) 

# Melt it and construct ranges
df <- melt(df, id.vars=c('Chromosome', 'Start', 'End'), variable.name='ID',
           value.name='state')



