library(ShortRead)
library(Rsubread)

# Some notes
# Rsubread align/subjunc appear to have a 32 threads limit
# Only part of their internal pipeline really uses threads
# I'm struggling to get proper template folder trees

# Set BiocParallel and get control over it
nthreads <- 32
smpParam <- MulticoreParam(workers=nthreads)
snowParam <- SnowParam(workers=nthreads, type='SOCK')
register(smpParam, default=TRUE)
register(snowParam)

# Initially, doing things on a sample basis
# Set paths and related stuff
baseDir <- '/data/Cytogenomics/TruSeq'
refDir <- paste(baseDir, 'genomes/UCSC/hg19', sep='/')
samplesDir <- paste(baseDir, user, 'samples', sep='/')
outDir <- paste(baseDir, user, 'results', sep='/')

# Create a place for results
dir.create(outDir)

# Read reference data
user <- 'Marília'
sampleSheetFilename <- 'sampleSheetTruSeq.csv'
sampleSheet <- read.csv(paste(baseDir, user, sampleSheetFilename, sep='/'),
                        as.is=TRUE)

# Quality control
qcDir <- paste(outDir, 'QC', sep='/')
dir.create(qcDir)
for(sample in sampleSheet$SampleID){

    QC <- qa(dirPath=samplesDir, pattern=paste0(sample, '_*'), type='fastq')
    tmp <- paste(qcDir, sample, sep='/')
    report(QC, dest=tmp)

}

# It seems that all samples are OK
# Let's map them using hg19 and Rsubread!!!
# As of 2014, trimming/filtering are unnecessary
# Well, I expect a lot of CIGAR alignment tweaking

idxName <- 'hg19.Rsubread.idx' 
idxDir <- paste(refDir, 'Rsubread', idxName, sep='/')
genomeFile <- paste(refDir, 'fa/hg19.fa', sep='/')
variantFile <- paste(refDir, 'vcf/00-All.vcf', sep='/')

bamDir <- paste(outDir, 'bam', sep='/')
dir.create(bamDir)

vcfDir <- paste(outDir, 'vcf', sep='/')
dir.create(vcfDir)

# Map read, takes about 1 hour per sample
for(sample in sampleSheet$SampleID){

    # Align reads for DE
    sampleFile <- Sys.glob(paste(samplesDir, paste0(sample, '*'), sep='/'))
    alnFile <- paste0(sample, '.aln.bam')
    align(index=idxDir,
          readfile1=sampleFile,
          input_format='gzfastq',
          output_file=paste(bamDir, alnFile, sep='/'),
          tieBreakQS=TRUE,
          nthreads=nthreads)

}

# Find splicing junctions, takes 2-3 hours per sample
for(sample in sampleSheet$SampleID){

    # Align reads for splicing junction discovery
    subjuncFile <- paste0(sample, '.subjunc.bam')
    subjunc(index=idxDir, 
            readfile1=sampleFile,
            input_format='gzfastq',
            output_file=paste(bamDir, subjuncFile, sep='/'),
            tieBreakQS=TRUE,
            reportAllJunctions=TRUE,
            nthreads=nthreads)

}

# For some reason, exactSNP isn't limited
# So, let's use more threads!!!
# Takes about 10-15 min per sample
for(sample in sampleSheet$SampleID){

    # SNP discovery
    snpFile <- paste0(sample, '.snp.vcf')
    exactSNP(readFile=sampleFile,
             refGenomeFile=genomeFile,
             SNPAnnotationFile=variantFile,
             isBAM=TRUE,
             outputFile=paste(vcfDir, snpFile, sep='/'),
             nthreads=60)

    # Remove residuals bin files
    unlink(Sys.glob('./*.bin'))

}


