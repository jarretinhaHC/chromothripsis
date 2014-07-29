# load configuration stuff
source('conf/chromothripsis.rconf')

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
crlmmCopyNumber(gen)

# Downstream analysis
# Data conversion to BafLrrSetList (bls)
bls <- BafLrrSetList(gen)

# Correct for genomic waves
wins <- c(1e2, 1e4, 1e6)
incs <- c(2e1, 2e3, 2e5)
ArrayTV:::gcCorrect(bls[-(1:1)], maxwins=wins, increms=incs)

# Fit HMM to identify CNVs
fit <- hmmBafLrrSetList2(bls[-(1:1)])

# Downstream analysis
ranges <- unlist(fit)

# Need to adjust chromosome names
chrs <- paste('chr', chromosome(bls), sep='')
filtered_bls <- bls[chrs %in% chromosome(ranges)]
# Filter those chromosomes that generated ranges
IDs <- unique(sampleNames(ranges))

#Find overlapping ranges by state

for(st in 1:6){

	filtered_ranges <- ranges[state(ranges) == st, ]
	tiles <- paste('all.overlapping.', st, '.ranges.dat', sep='')
	tsv <- paste('all.overlapping.', st, '.ranges.tsv', sep='')
	hist <- paste('all.overlapping.', st, '.ranges.hist', sep='')
	
	tmp <- as.data.frame(filtered_ranges, row.names=NULL)
	tmp$group <- subjectHits(findOverlaps(filtered_ranges, reduce(filtered_ranges)))
	tmp <- tmp[which(tmp$width > 1), ]
	tmp <- as.data.table(tmp)
	tmp <- tmp[, list(start=min(start), end=max(end), samples=list(unique(sample))), by=list(group, seqnames)]

	tmp$counts <- sapply(tmp$samples, FUN=length)
	tmp$samples <- sapply(tmp$samples, FUN=paste, collapse='|')
	
	data <- data.frame('CHR'=gsub('chr', 'hs', tmp$seqnames), 'START'=tmp$start, 'END'=tmp$end)
	data$STATE <- paste('state=', rep(st, dim(data)[1]), sep='')
	write.table(data, file=tiles, sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)

	data <- data.frame('CHR'=gsub('chr', 'hs', tmp$seqnames), 'START'=tmp$start, 'END'=tmp$end, COUNTS=tmp$counts * 10)
	write.table(data, file=hist, sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)


	write.table(tmp, file=tsv, sep='\t', row.names=FALSE, quote=FALSE)

}

for(ID in IDs){

	# Filter by sample ID
	filtered_ranges <- ranges[state(ranges) %in% c(1, 2, 5, 6) & sampleNames(ranges) == ID, ]
	filename <- paste(ID, '.ranges.dat', sep='')
	
	tmp <- as.data.frame(filtered_ranges, row.names=NULL)

	idx <- which(tmp$width > 1)	
	CHR <- gsub('chr', 'hs', tmp$seqnames[idx])
	START <- tmp$start[idx]
	END <- tmp$end[idx]
	STATE <- paste('state=', tmp$state[idx], sep='')
	#COLOR <- paste('color=', c('red', 'orange', 'white', 'white', 'cyan', 'blue')[STATE], sep='')

	
	data <- data.frame('CHR'=CHR, 'START'=START, 'END'=END, 'STATE'=STATE)
	#data <- data.frame('CHR'=CHR, 'START'=START, 'END'=END, 'STATE'=STATE, 'COLOR'=COLOR)
	
	write.table(data, file=filename, sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)

}

for(ID in IDs){

	# Filter by sample ID
	filtered_ranges <- ranges[sampleNames(ranges) == ID, ]
	tmp <- filtered_bls[ , match(ID, sampleNames(filtered_bls))]

	dataLRR <- data.frame('CHR'=NA, 'START'=NA, 'END'=NA, 'LRR'=NA)
	dataBAF <- data.frame('CHR'=NA, 'START'=NA, 'END'=NA, 'BAF'=NA)
	for(i in 1:length(tmp)){

		CHR <- unlist(paste('hs', chromosome(tmp[i]), sep=''))
		START <- unlist(position(tmp[i]))
		END <- unlist(position(tmp[i]))
		LRR <- unlist(lrr(tmp[i]))/100
		BAF <- unlist(baf(tmp[i]))/1000
		#COLOR <- ifelse(LRR == 0, 'color=blue', 'color=red')
		dflrr <- data.frame('CHR'=CHR, 'START'=START, 'END'=END, 'LRR'=LRR)
		dfbaf <- data.frame('CHR'=CHR, 'START'=START, 'END'=END, 'BAF'=BAF)
		dataLRR <- rbind(dataLRR, dflrr)
		dataBAF <- rbind(dataBAF, dfbaf)
		
	}

	filename <- paste(ID, '.LRR.dat', sep='')
	write.table(dataLRR, file=filename, sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)
	
	filename <- paste(ID, '.BAF.dat', sep='')
	write.table(dataBAF, file=filename, sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)

}

