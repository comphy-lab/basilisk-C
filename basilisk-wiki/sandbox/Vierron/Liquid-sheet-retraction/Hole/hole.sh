#!/bin/bash
MAXLEVEL=10
process_id=$!
for Oh in '0.1' '0.01' 
do
    for H in '1.'
    do
	mkdir hole$Oh$H
	cp hole RP.py `pwd`/hole$Oh$H/
	cd hole$Oh$H
	mpirun -np 30 ./hole $Oh $H $MAXLEVEL &
	cd ..
    done
done
wait $process_id



