#!/bin/bash

mkdir -p CelllineFeatures

# Download original data 
echo Downloading DREAM data and making features ...
./getdata.sh
./getAllEC10.sh
echo Done !

# Download SNP data
echo Downloading SNP euclidean distances ...
./getsnpdist.sh
echo Done !

# Download chemical descriptors
echo Downloading chemical descriptors ...
./getchemdescriptors.sh
echo Done ! 

mkdir -p CelllineKernels
mkdir -p ChemicalKernels

# Create some formatted feature files
echo Creating some formatted features files ...
Rscript makefeatures.R
echo Done ! 

# Make kernel matrices
echo Making kernels ...
Rscript makekernels.R
echo Done !
