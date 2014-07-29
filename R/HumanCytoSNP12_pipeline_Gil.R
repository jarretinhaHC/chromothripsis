library(ff)
library(humancytosnp12v2p1hCrlmm)
library(oligoClasses)
library(Biobase)
library(IRanges)
library(crlmm)
library(VanillaICE)
library(ArrayTV)
library(data.table)

options(ffcaching="ffeachflush")

#Set results folder
outdir <- paste('/data/Cytogenomics/HumanCyto12/scratch/R-', getRversion(), '/r00', sep='')
ldPath(outdir)
dir.create(outdir, recursive=TRUE)

# Set data folder and sheet
datadir <- '/data/Cytogenomics/HumanCyto12/IDAT'
sampleSheet <- read.delim('/data/Cytogenomics/HumanCyto12/SampleSheet/simples.tsv', as.is=TRUE)

arrayNames <- file.path(datadir, sampleSheet$Sample_ID)
arrayInfo <- list(barcode=NULL, position='SentrixPosition')
cdfName <- 'humancytosnp12v2p1h'
batch <- rep("1", nrow(sampleSheet))

print(all(file.exists(paste(arrayNames, "_Grn.idat", sep=""))))
print(all(file.exists(paste(arrayNames, "_Red.idat", sep=""))))

# Genotyping
gen <- genotype.Illumina(sampleSheet=sampleSheet, arrayNames=arrayNames, arrayInfoCol=arrayInfo, cdf=cdfName, batch=batch)

# Copy number estimation
crlmmCopynumber(gen)

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
IDs <- c('7685541145_R01C01', '7685541145_R03C01', '7685541145_R04C01', '7685541145_R05C02', '7685541057_R01C01')

#Find overlapping ranges by state

for(st in c(1, 2, 3, 4, 5, 6)){

	filtered_ranges <- ranges[state(ranges) == st & sampleNames(ranges) %in% IDs, ]
	
	if(length(ranges(filtered_ranges)) == 0) next
	
	tiles <- paste('Gil.overlapping.', st, '.ranges.dat', sep='')
	tsv <- paste('Gil.overlapping.', st, '.ranges.tsv', sep='')
	hist <- paste('Gil.overlapping.', st, '.ranges.hist', sep='')
	
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
	filtered_ranges <- ranges[state(ranges) %in% c(1, 2, 3, 4, 5, 6) & sampleNames(ranges) == ID, ]
	
	if(length(ranges(filtered_ranges)) == 0) next
	
	filename <- paste('Gil.', ID, '.ranges.dat', sep='')
	
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
	
	if(length(ranges(filtered_ranges)) == 0) next
	
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

	filename <- paste('Gil.', ID, '.LRR.dat', sep='')
	write.table(dataLRR, file=filename, sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)
	
	filename <- paste('Gil.', ID, '.BAF.dat', sep='')
	write.table(dataBAF, file=filename, sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)

}

