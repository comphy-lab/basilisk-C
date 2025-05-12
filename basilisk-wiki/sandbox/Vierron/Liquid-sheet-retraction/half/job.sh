#!/bin/bash                
#SBATCH -n 4                           
#SBATCH -p batch576                    

MAXLEVEL=10
process_id=$!
for Oh in '0.0001' '0.0005' '0.001' '0.005' '0.01' '0.05' '0.1'
do
    for H in '1.' '0.9' '0.8' '0.7' '0.6' '0.5'
    do
	mkdir half$Oh$H
	cp half h.py RP.py `pwd`/half$Oh$H/
	cd half$Oh$H
	srun --mpi=pmi2 -J half${Oh}${H} ./half $Oh $H $MAXLEVEL &
	avconv -r 10 -vcodec ppm -f image2pipe -i Hole${Oh}.mp4.ppm -vcodec mpeg
4 -vb 40M -r 24 Hole.mp4 &
	avconv -r 10 -vcodec ppm -f image2pipe -i grid${Oh}.mp4.ppm -vcodec mpeg
4 -vb 40M -r 24 grid.mp4 &
	rm *.ppm &
	cd ..
    done
done
wait $process_id
