## FILE
##   getsnpdist.sh
## DESCRIPTION
##   creates the GenotypeDistance folder with a ncells X ncells matrix of snp euclidean distance

#!/bin/bash

SERVER=http://members.cbio.ensmp.fr/~twalter/data/tox_challenge_data/
DATA="SNP_distance_all"

# Download and extract data
for d in ${DATA}; do
    if [ ! -d ${d} ]; then
        echo Downloading ${d} from ${SERVER}
	curl -O ${SERVER}${d}.tgz
        tar xvzf ${d}.tgz
    fi
done

