#/bin/bash

#Create reference file for test cases.

for FILE in *.c; do
    export FILENAME="${FILE%.*}";

    #Check for no convergence
    grep -q "#" ${FILENAME}/log
    if [ $? == 0 ]; then
	echo "Convergence problem in "${FILENAME}".c";
    fi
    
    # Check for difference with .ref
    diff ${FILENAME}/log ${FILENAME}.ref
    if [ $? -ne 0 ]; then
	echo ${FILENAME};
    fi
done
