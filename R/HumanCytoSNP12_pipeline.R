# load configuration stuff
source('conf/chromothripsis.rconf')

# Let's parallelize calculations.
registerDoMC(12)
ocSamples(12)

# Genotyping
gen <- genotype.Illumina(sampleSheet=sampleSheet, arrayNames=arrayNames, arrayInfoCol=arrayInfo, cdf=cdfName, batch=batch, sns=sns, col.names=TRUE)

# Check data quality
snr <- cnSet$SNR[]
barcodes <- rapply(strsplit(sampleNames(cnSet), '_'), function(x) head(x, 1))
snr <- data.frame(SNR=snr, Array=barcodes)

pdf('DataQualityHistogramOverall.pdf')
plt <- ggplot(data=snr, aes(SNR)) + geom_bar(binwidth=1, colour='black', alpha='0.75') + ylab('Amostras')
print(plt)
dev.off()

pdf('DataQualityHistogramArray.pdf')
plt <- ggplot(data=snr, aes(SNR, fill=Array)) + geom_bar(binwidth=1, colour='black', alpha='0.75') + ylab('Amostras') + facet_wrap(~Array)
print(plt)
dev.off()

# Copy number estimation
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

window_size <- 10^(2:6)
increments <- 2 * 10^(1:5)
ArrayTV:::gcCorrect(cleaned_bls, maxwins=window_size, increms=increments)

# Fit HMM to identify CNVs
fit <- unlist(hmmBafLrrSetList2(cleaned_bls))

# Back to serial mode.
stopCluster(cluster)

# Downstream analysis
# Idea is to reduce ranges to a single data track and add a coverage track

for(st in 1:6){

    ranges <- fit[state(fit) == st, ]
    ranges_with_coverage <- as(coverage(ranges), 'GRanges')
    reduced_ranges <- reduce(ranges)
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
    df$strand <- c()
    df <- df[df$score > 0 || df$width == 1, ]

    # Table data/Coverage
    write.table(df, file=paste0('ranges_with_coverage_', st, '.tsv'), sep='\t', row.names=FALSE, quote=FALSE)

    # Circos/Coverage
    data <- data.frame('CHR'=gsub('chr', 'hs', df$seqnames),
                                'START'=df$start, 'END'=df$end)
	data$SCORE <- paste0('score=', df$score)
	write.table(data, file=paste0('ranges_with_coverage_for_circos_', st,
                                 '.dat'), sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)

    # Circos/Reduced range
	df <- as.data.frame(reduced_ranges, row.names=NULL)
    df$strand <- c()
    df <- df[df$width == 1, ]

    data <- data.frame('CHR'=gsub('chr', 'hs', df$seqnames),
                       'START'=df$start, 'END'=df$end)
	write.table(data, file=paste0('reduced_ranges_for_circos_', st,
                                 '.dat'), sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)

}


