# FILE
#   makekernels.R
# DESCRIPTION
#   Creates kernel files for chemicals and cell lines in folders:
#   CelllineKernels, ChemicalKernels

# Directory to store chemicals' features
chempath <- 'ChemicalKernels/'
chemfile <- function(filename) {paste(chempath,filename,sep="")}

# Directory to store cell lines' features
cellpath <- 'CelllineKernels/'
cellfile <- function(filename) {paste(cellpath,filename,sep="")}

# Bandwidth of Gaussian kernels (as fraction of the maximum distances between points)
siglist <- c(1/20 , 1/10 , 1/5 , 2/5 , 3/5 , 4/5 , 1, 2, 4, 8)

### Cell lines


# Covariates
cat('Covariate kernel...')
x <- read.table(file='CelllineFeatures/covariates.txt',header=TRUE)
g <- as.data.frame(t(x))
a <- cbind(model.matrix(~ 0 + Sex , g) , model.matrix( ~ 0 + Population, g), model.matrix(~ 0 + Cytotoxicity_Batch , g) )
write.table(a%*%t(a)/3,cellfile('Kcovariates.txt'),quote=FALSE,sep='\t')
a <- model.matrix(~ 0 + Sex , g)
write.table(a%*%t(a),cellfile('KcovariatesSex.txt'),quote=FALSE,sep='\t')
a <- model.matrix(~ 0 + Population , g)
write.table(a%*%t(a),cellfile('KcovariatesPopulation.txt'),quote=FALSE,sep='\t')
a <- model.matrix(~ 0 + Cytotoxicity_Batch , g)
write.table(a%*%t(a),cellfile('KcovariatesBatch.txt'),quote=FALSE,sep='\t')
cat('done\n')


# RNA-seq
cat('RNAseq kernels...')
x <- read.table(file='CelllineFeatures/RNASeqnormalized_withNA.txt',header=TRUE)
allcell <- colnames(x)
ncell <- length(allcell)
valid <- (apply(is.na(x),2,sum)==0)
x <- x[,valid]
# log transform
v <- log10(x + min(x[x>0]))
d <- dist(t(v))
dm <- as.matrix(d)
# Hack : Since distances range between 50 and 100, we subtract the smallest value to create the Gaussian kernel and will later add back to the diagonal of the kernel a scalar large enough to make is positive semidefinite.
dd <- dm^2 - min(dm[dm>0])^2
diag(dd) <- 0
slist <- max(dd) * siglist
kall <- matrix(0,nrow=ncell,ncol=ncell,dimnames=list(allcell,allcell))
for (i in seq(length(slist))) {
	k <- exp(-dd/slist[i])
	# make it pd by removing the smallest eigenvalue from the diagonal, then scale the diagonal to 1
	k <- k/(1-min(eigen(k)$values))
	kall[valid,valid] <- k
	diag(kall) <- 1
	write.table(kall , cellfile(paste('KrnaseqRbf',i,'.txt',sep='')),quote=FALSE,sep='\t')
}
cat('done\n')

# Genotype
cat('SNP kernel...')
dm <- read.table(file='GenotypeDistance/SNP_eu_distance_all.txt',header=FALSE)
# Hack : Since distances range between 293k and 491k, we subtract the smallest value to create the Gaussian kernel and will later add back to the diagonal of the kernel a scalar large enough to make is positive semidefinite.
dd <- dm-min(dm[dm>0])
diag(dd) <- 0
slist <- max(dd) * siglist
for (i in seq(length(slist))) {
	kall <- exp(-dd/slist[i])
	dimnames(kall) <- list(allcell,allcell)
	# make it pd by removing the smallest eigenvalue from the diagonal, then scale the diagonal to 1
	kall <- kall/(1-min(eigen(kall)$values))
	diag(kall) <- 1
	write.table(kall , cellfile(paste('KsnpRbf',i,'.txt',sep='')),quote=FALSE,sep='\t')
}
cat('done\n')


### Chemicals

# CDK
cat('CDK kernel...')
x <- read.table(file='ChemicalFeatures/cdk.txt',header=TRUE)
a <- scale(t(x))

dd <- as.matrix(dist(a))
slist <- max(dd) * siglist
for (i in seq(length(slist))) {
	write.table(exp(-dd/slist[i]) , chemfile(paste('KcdkRbf',i,'.txt',sep='')),quote=FALSE,sep='\t')
}
cat('done\n')

# SIRMS
cat('SIRMS kernel...')
x <- read.table(file='ChemicalFeatures/sirms.txt',header=TRUE)
a <- scale(t(x))

dd <- as.matrix(dist(a))
slist <- max(dd) * siglist
for (i in seq(length(slist))) {
	write.table(exp(-dd/slist[i]) , chemfile(paste('KsirmsRbf',i,'.txt',sep='')),quote=FALSE,sep='\t')
}
cat('done\n')

# Predicted targets?
cat('Predicted target kernel...')
x <- read.table("ChemicalFeatures/pred_matrix_vero.csv",header=TRUE,row.names=1)
row.names(x) <- make.names(row.names(x))
a <- scale(x)
dd <- as.matrix(dist(a))
slist <- max(dd) * siglist
for (i in seq(length(slist))) {
	write.table(exp(-dd/slist[i]) , chemfile(paste('KpredtargetRbf',i,'.txt',sep='')),quote=FALSE,sep='\t')
}
cat('done\n')


# Substructures
cat('Substructure kernel...')
x <- read.table("ChemicalFeatures/sub_strucure_features.csv",header=TRUE,row.names=1)
row.names(x) <- make.names(row.names(x))
k <- 1 - as.matrix(dist(x,"binary"))
write.table(k , chemfile("Ksubstructure.txt"),quote=FALSE,sep='\t')
cat('done\n')

# ChemCPP
cat('ChemCPP kernel...')
x <- read.table("ChemicalFeatures/sd2gram.csv",header=TRUE,row.names=1)
rownames(x) <- make.names(row.names(x)) 
write.table(x , chemfile("Kchemcpp.txt"),quote=FALSE,sep='\t')
cat('done\n')


