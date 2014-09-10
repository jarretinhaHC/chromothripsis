library(Rsubread)
library(DESeq2)

# Set BiocParallel and get control over it
nthreads <- 32
smpParam <- MulticoreParam(workers=nthreads)
snowParam <- SnowParam(workers=nthreads, type='SOCK')
register(smpParam, default=TRUE)
register(snowParam)

# Initially, doing things on a sample basis
# Set paths and related stuff

user <- 'Marília'
sampleSheetFilename <- 'sampleSheetTruSeq.csv'

baseDir <- '/data/Cytogenomics/TruSeq'
resDir <- paste(baseDir, user, 'results', sep='/')
bamDir <- paste(resDir, 'bam', sep='/')
deDir <- paste(resDir, 'DE', sep='/')
dir.create(deDir)

# Read reference data
sampleSheet <- read.csv(paste(baseDir, user, sampleSheetFilename, sep='/'),
                        as.is=TRUE)
 
# Let's count features
alnFiles <- Sys.glob(paste(bamDir, '*.aln.bam', sep='/'))

# Features grouped by gene models
fcGene <- featureCounts(files=alnFiles,
                       annot.inbuilt='hg19',
                       nthreads=nthreads)

# Let's construct DESeq2 data container
# Control will be the reference level (Group)
# In this particular case rownames will be SampleID
countData <- fcGene$counts
colnames(countData) <- sampleSheet$SampleID
rownames(sampleSheet) <- sampleSheet$SampleID
ddsGene <- DESeqDataSetFromMatrix(countData, sampleSheet, ~Group)
ddsGene$Group <- relevel(ddsGene$Group, 'Control')

# Run analysis pipeline
ddsGene <- DESeq(ddsGene)


