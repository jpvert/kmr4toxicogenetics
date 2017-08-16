## FILE
##   getchemdescriptors.sh
## DESCRIPTION
##   creates the ChemicalFeatures folder containing hand-tailored chemical features

#/bin/bash

SERVER=http://members.cbio.ensmp.fr/~twalter/data/tox_challenge_data/
#DATA="convert_prot_gene additional_chemical_descriptors"
#WARNING: the convert_prot_name download does not work

DATA="additional_chemical_descriptors"

# Download and extract data
for d in ${DATA}; do
    if [ ! -d ${d} ]; then
        echo Downloading ${d}.tgz from ${SERVER}
	curl -O ${SERVER}${d}.tgz
        mkdir -p ./ChemicalFeatures
        mv ${d}.tgz ./ChemicalFeatures/
        cd ./ChemicalFeatures
        tar xvzf ${d}.tgz
	rm ${d}.tgz
    fi
done

