set -x
mkdir helix2_base/
v v v v v v v
for level in 8 ; do

# mkdir helix2_base/level${level}
# mkdir helix2_base/level${level}/000
# make clean
# CFLAGS='-march=native -fopenmp -DMAXLEVEL='"${level}"' -D_INIT -DMTRACE=2 -DTRACE=2' make helix2_base.tst
# mv -v ./helix2_base/dump0 ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/000/
=============
for level in 8 9 ; do

mkdir helix2_base/level${level}
mkdir helix2_base/level${level}/000
*************
for level in 9 10 11 12 ; do

mkdir helix2_base/level${level}
mkdir helix2_base/level${level}/000
^ ^ ^ ^ ^ ^ ^
make clean
v v v v v v v


mkdir helix2_base/level${level}/001
cp ./helix2_base/level${level}/000/dump0 ./helix2_base/dump0
CFLAGS='-march=native -fopenmp -D_TMAX=2.0 -DMAXLEVEL='"${level}"' -DMTRACE=2 -DTRACE=2 ' make helix2_base.tst
mv -v ./helix2_base/dump-* ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/001/
=============
CFLAGS='-march=native -fopenmp -DMAXLEVEL='"${level}"' -D_INIT -DMTRACE=2 -DTRACE=2' make helix2_base.tst
#CFLAGS='-march=native -fopenmp -DMAXLEVEL='"${level}"' -DMTRACE=2 -D_TMAX=0.1 -DTRACE=2' make helix2_base.tst
mv -v ./helix2_base/dump0 ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/000/
cp ./helix2_base/level${level}/000/dump0 ./helix2_base/level${level}/dump0
*************
CFLAGS='-march=native -fopenmp -DMAXLEVEL='"${level}"' -D_INIT -DMTRACE=2 -DTRACE=2' make helix2_base.tst
mv -v ./helix2_base/dump0 ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/000/
cp ./helix2_base/level${level}/000/dump0 ./helix2_base/level${level}/
^ ^ ^ ^ ^ ^ ^
make clean

v v v v v v v
mkdir helix2_base/level${level}/002
cp ./helix2_base/level${level}/001/dump-2 ./helix2_base/dump0
CFLAGS='-march=native -fopenmp -D_TMAX=4.0 -DMAXLEVEL='"${level}"' -DMTRACE=2 -DTRACE=2 ' make helix2_base.tst
mv -v ./helix2_base/dump-* ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/002/
make clean

mkdir helix2_base/level${level}/003
cp ./helix2_base/level${level}/002/dump-4 ./helix2_base/dump0
CFLAGS='-march=native -fopenmp -D_TMAX=6.0 -DMAXLEVEL='"${level}"' -DMTRACE=2 -DTRACE=2 ' make helix2_base.tst
mv -v ./helix2_base/dump-* ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/003/
make clean
=============

# mkdir helix2_base/level${level}/001
# cp ./helix2_base/level${level}/000/dump0 ./helix2_base/dump0
# CFLAGS='-march=native -fopenmp -D_TMAX=2.0 -DMAXLEVEL='"${level}"' -DMTRACE=2 -DTRACE=2 ' make helix2_base.tst
# mv -v ./helix2_base/dump-* ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/001/
# make clean
#
# mkdir helix2_base/level${level}/002
# cp ./helix2_base/level${level}/001/dump-2 ./helix2_base/dump0
# CFLAGS='-march=native -fopenmp -D_TMAX=4.0 -DMAXLEVEL='"${level}"' -DMTRACE=2 -DTRACE=2 ' make helix2_base.tst
# mv -v ./helix2_base/dump-* ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/002/
# make clean
#
# mkdir helix2_base/level${level}/003
# cp ./helix2_base/level${level}/002/dump-4 ./helix2_base/dump0
# CFLAGS='-march=native -fopenmp -D_TMAX=6.0 -DMAXLEVEL='"${level}"' -DMTRACE=2 -DTRACE=2 ' make helix2_base.tst
# mv -v ./helix2_base/dump-* ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/003/
# make clean
*************

#mkdir helix2_base/level${level}/001
#cp ./helix2_base/level${level}/000/dump0 ./helix2_base/dump0
#CFLAGS='-march=native -fopenmp -D_TMAX=2.0 -DMAXLEVEL='"${level}"' -DMTRACE=2 -DTRACE=2 ' make helix2_base.tst
#mv -v ./helix2_base/dump-* ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/001/
#make clean

#mkdir helix2_base/level${level}/002
#cp ./helix2_base/level${level}/001/dump-2 ./helix2_base/dump0
#CFLAGS='-march=native -fopenmp -D_TMAX=4.0 -DMAXLEVEL='"${level}"' -DMTRACE=2 -DTRACE=2 ' make helix2_base.tst
#mv -v ./helix2_base/dump-* ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/002/
#make clean

#mkdir helix2_base/level${level}/003
#cp ./helix2_base/level${level}/002/dump-4 ./helix2_base/dump0
#CFLAGS='-march=native -fopenmp -D_TMAX=6.0 -DMAXLEVEL='"${level}"' -DMTRACE=2 -DTRACE=2 ' make helix2_base.tst
#mv -v ./helix2_base/dump-* ./helix2_base/*.vts ./helix2_base/*.asc ./helix2_base/*.mp4 ./helix2_base/*.png ./helix2_base/log ./helix2_base/out ./helix2_base/level${level}/003/
#make clean
^ ^ ^ ^ ^ ^ ^
done
