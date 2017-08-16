# FILE
#   makefeatures.R
# DESCRIPTION
#   Creates feature files for chemicals and cell lines in folders:
#   CelllineFeatures, ChemicalFeatures
# more specifically, we have:
#   CelllineFeatures/covariates.txt: a matrix 3 X ncells with columns of attributes: sex, population and batch
#   CelllineFeatures/RNASeqnormalized_withNA.txt and CelllineFeatures/RNASeqnormalized_noNA.txt: 46256 X ncells normalized counts and there is only RNA-seq data for 337 cells
#   ChemicalFeatures/cdk.txt and ChemicalFeatures/sirms.txt, NB: + the snp matrix see above in GenotypeDistance folder
# NOTE
#   Each feature file should be a tab-separated file
#   First row = name of the chemicals or cell lines
#   Following rows : name of the attribute, then attributes for all chemicals/cell lines

# Directory to store chemicals' features
chempath <- 'ChemicalFeatures/'
system(paste("mkdir -p", chempath))
chemfile <- function(filename) { paste(chempath,filename,sep="") }

# Directory to store cell lines' features
cellpath <- 'CelllineFeatures/'
system(paste("mkdir -p", cellpath))
cellfile <- function(filename) { paste(cellpath,filename,sep="") }


### Cell lines

# Covariates provided by DREAM8
x <- read.table(file='Covariates/ToxChallenge_Covariates.txt', header=TRUE, row.names=1)
write.table(t(x), file=cellfile('covariates.txt'), quote=FALSE, sep='\t')

# RNA-seq
x <- read.table(file='RNASeq/ToxChallenge_RNASeq_counts_norm.txt', header=TRUE, row.names=1)
write.table(x,file=cellfile('RNASeqnormalized_withNA.txt'), quote=FALSE, sep="\t")
z=apply(x, 2, sum)
indna <- which(is.na(z))
xx <- x[,-indna]
write.table(xx,file=cellfile('RNASeqnormalized_noNA.txt'), quote=FALSE, sep="\t")


### Chemicals

# CDK descriptors (provided by DREAM)
x <- read.table('ChemicalAttributes/ToxChallenge_ChemicalAttributes_CDK.txt',header=TRUE)
write.table(t(x),file=chemfile('cdk.txt'),quote=FALSE,sep="\t")

# SiRMS descriptors (provided by DREAM)
x <- read.table('ChemicalAttributes/ToxChallenge_ChemicalAttributes_SiRMS.txt',header=TRUE)
write.table(t(x),file=chemfile('sirms.txt'),quote=FALSE,sep="\t")
