#/bin/bash

#Create reference file for test cases.

for FILE in *.c; do
    export FILENAME="${FILE%.*}";
    cp ${FILENAME}/log ${FILENAME}.ref
done
