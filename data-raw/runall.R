# FILE
#   runall.R
# DESCRIPTION
#   download dream8toxicogenetics data and some hand-tailored chemical descriptors from own server, 
#   then generate features and kernels as saved in data/

# setwd("data-raw")

# download data, descriptors and generate features and kernels ------------
# data to be saved locally in each folder

# system("runall_data_kernels.sh")

# process features --------------------------------------------------------

# toxicity - response - sizeof ncell x nchem
tox1 = read.table("EC10/ToxChallenge_CytotoxicityData_Train_Subchal1_Extended.txt", sep="\t", stringsAsFactors=F)
tox2 = read.table("EC10/ToxChallenge_Subchall1_Test_Data.txt", sep="\t", stringsAsFactors=F)
tox = rbind(tox1, tox2)
CellNames = rownames(tox)
ChemNames = colnames(tox)
tox = as.matrix(tox)
cat("Number of total chemicals = ", length(ChemNames))
cat("Number of total cell lines = ", length(CellNames))
# remove non predictable compounds [given challenge paper]
toremove = kmr4toxicogenetics:::nontoxic()
tox = tox[,!(colnames(tox) %in% toremove)]
keepchem = colnames(tox)
cat("Number of in-use chemicals = ", length(keepchem))
# keep training or full set of cell lines
keepcell = CellNames
tox = tox[keepcell,]
cat("Number of in-use cell lines = ", length(keepcell))
# toxicity data dim
cat("Toxicity response dim = ", dim(tox))

# cell lines - features - sizeof ncell x p
# 1) covariates
x  = read.table(file='CelllineFeatures/covariates.txt',header=TRUE, row.names=1)
g = as.data.frame(t(x))
covariates = cbind( model.matrix(~ 0 + Sex , g) , model.matrix( ~ 0 + Population, g), model.matrix(~ 0 + Cytotoxicity_Batch , g) )
# 2) SNP gram
SNPgram = read.table(file='GenotypeDistance/SNP_eu_distance_all.txt', stringsAsFactors=F)
rownames(SNPgram) = CellNames
colnames(SNPgram) = paste0("SNPgram",1:ncol(SNPgram))
# 3) RNA gram
RNA_seq = read.table(file='CelllineFeatures/RNASeqnormalized_withNA.txt', header = TRUE, row.names=1)
RNAgram = crossprod(x=as.matrix(RNA_seq), y=NULL)
ind.na <- which(is.na(RNAgram[1,]))
RNAfix <- RNAgram[,-ind.na]
RNAfix <- randomForest::na.roughfix(RNAfix)
colnames(RNAfix) = paste0("RNAgram",1:ncol(RNAfix))
# combine
design = cbind(covariates, SNPgram, RNAfix)
stopifnot(all(rownames(design) == CellNames))
# resize
design = design[keepcell,]
# cell line design matrix size
cat("Cell line design matrix dim = ", dim(design))
# different feature categories
cat("Number of dim in each feature category = ")
table(as.vector( sapply(colnames(design), function(a) if (!grepl("gram",a)) {return("covariate")} else {return(paste(strsplit(a,split="")[[1]][1:7],collapse=""))}) ))

# save up
devtools::use_data(design, tox)

# process kernels ---------------------------------------------------------

# Read kernel files for cell lines
pathcellkernel <- "CelllineKernels/"
kcellname <- list.files(path= pathcellkernel ,pattern="*.txt")
kcell <- lapply(kcellname,function(n) as.matrix(read.table(paste(pathcellkernel , n,sep=''),header=TRUE,row.names=1)))
kcellname <- gsub("[.]txt$","",kcellname)
names(kcell) <- kcellname
kcell <- lapply(kcell,function(k) k[keepcell,keepcell])
ncelllines <- length(keepcell)
# Get idx of RNA-seq cell lines only
id.rna <- keepcell[apply(kcell[[grep("rna", kcellname)[1]]][keepcell,],1,sum)>1]
cat("Number of in-use RNA-seq only cell lines = ", length(id.rna))

# Read kernel files for chemicals
pathchemkernel <- "ChemicalKernels/"
kchemname <- list.files(path= pathchemkernel ,pattern="*.txt")
kchem <- lapply(kchemname,function(n) as.matrix(read.table(paste(pathchemkernel , n,sep=''),header=TRUE,row.names=1)))
kchemname <- gsub("[.]txt$","",kchemname)
names(kchem) <- kchemname
kchem <- lapply(kchem,function(k) k[keepchem,keepchem])
nchemicals <- length(keepchem)

# save up
devtools::use_data(kcell, kchem, id.rna)
