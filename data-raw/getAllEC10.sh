## FILE
##   getAllEC10.sh
## DESCRIPTION
##   creates EC10 folder with response (EC10) available from Dream 8 toxicogenetics challenge

#!/bin/bash

SERVER=http://members.cbio.ensmp.fr/~ebernard/dream8_tox_data/EC10/
DATA="ToxChallenge_CytotoxicityData_Train_Subchal1_Extended.txt ToxChallenge_CytotoxicityData_Train_Subchal2_Extended.txt ToxSubchallenge_1_Test_Submission_File_Format.txt ToxChallenge_Subchall1_Test_Data.txt"

SERVER2=http://members.cbio.mines-paristech.fr/~jvert/svn/dream8toxicogenetics/data/
DATA2="EC10.tar.gz"

mkdir -p EC10

# Download and extract data
for d in ${DATA}; do
    if [ ! -d EC10/${d} ]; then
        echo Downloading ${d} from ${SERVER}
	curl -o EC10/${d} ${SERVER}${d}
    fi
done

for d in ${DATA2}; do
    if [ ! -d ${d} ]; then
        echo Downloading ${d} from ${SERVER2}
	curl -o ${d} ${SERVER2}${d}
	tar -xvzf ${d}
    fi
done
