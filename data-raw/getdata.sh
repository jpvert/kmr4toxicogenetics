## FILE
##   getdata.sh
## DESCRIPTION
##   creates the following folders containing raw data downloaded from Dream 8 toxicogenetics challenge 
##   ChemicalAttributes, Covariates, Genotype, RNASeq

#!/bin/bash

SERVER=http://members.cbio.ensmp.fr/~twalter/data/tox_challenge_data/
DATA="ChemicalAttributes Covariates Genotype RNASeq"
#DATA="ChemicalAttributes Covariates EC10 Genotype RNASeq"

# Download and extract data
for d in ${DATA}; do
    if [ ! -d ${d} ]; then
        echo Downloading ${d} from ${SERVER}
	curl -O ${SERVER}${d}.tar.gz
        tar xvzf ${d}.tar.gz
    fi
done
