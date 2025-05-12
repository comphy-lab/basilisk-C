include $(BASILISK)/Makefile.defs
export BASILISK_INCLUDE_PATH=/home/castillo/Workspace/wiki/sandbox/Antoonvh/


CC99 = gcc
CFLAGS += -Wall -std=c99 -O3
#CFLAGS += -Wall -std=c99 -fopenmp -O3
OPENGLIBS=-lglutils -lfb_osmesa -lGLU -lOSMesa -lgsl -lgslcblas
#OPENGLIBS=-lfb_glx -lGLU -lGLEW -lGL -lX11

export PROF_PATH=/home/castillo/Workspace/wiki/sandbox/acastillo/filaments/helical/

helical_relax:
	export OMP_NUM_THREADS=20;
	export KMP_NUM_THREADS=20;
	qcc -fopenmp -D_RELAX -source helical.c ;
	icc -fopenmp $(CFLAGS) _helical.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o helical/helical_icc -lm ;
	cd helical ; ./helical_icc > out 2> log ; cd ..

helical_frame:
	export OMP_NUM_THREADS=20;
	export KMP_NUM_THREADS=20;
	qcc -fopenmp -D_RELAX -D_FRAME -source helical.c ;
	icc -fopenmp  $(CFLAGS) _helical.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o helical/helical_icc -lm ;
	cd helical ; ./helical_icc > out 2> log ; cd ..

helical_relax_mpi:
	export OMP_NUM_THREADS=1;
	export KMP_NUM_THREADS=1;
	qcc -D_MPI=20 -D_RELAX -source helical.c ;
	mpiicc $(CFLAGS) _helical.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o helical/helical_mpiicc -lm ;
	cd helical ; mpirun -np 20 -genv I_MPI_PIN_PROCESSOR_LIST=allcores ./helical_mpiicc > out 2> log ; cd ..

helical_frame_mpi:
	export OMP_NUM_THREADS=1;
	export KMP_NUM_THREADS=1;
	qcc -D_MPI=20 -D_RELAX -D_FRAME -source helical.c ;
	mpiicc $(CFLAGS) _helical.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o helical/helical_mpiicc -lm ;
	cd helical ; mpirun -np 20 -genv I_MPI_PIN_PROCESSOR_LIST=allcores ./helical_mpiicc > out 2> log ; cd ..

twix_relax:
	export OMP_NUM_THREADS=16;
	export KMP_NUM_THREADS=16;
	qcc -fopenmp $(QCFLAGS) -source twin_helices.c ;
	icc -fopenmp $(CFLAGS) _twin_helices.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o twin_helices/twin_helices -lm ;
	cd twin_helices ; ./twin_helices > out 2> log ; cd ..

twix_mpi:
		export OMP_NUM_THREADS=1;
		export KMP_NUM_THREADS=1;
		qcc -D_MPI=16 $(QCFLAGS) -source twin_helices.c ;
		mpiicc $(CFLAGS) _twin_helices.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o twin_helices/twin_helices -lm ;
		cd twin_helices ; mpirun -np 16 -genv I_MPI_PIN_PROCESSOR_LIST=allcores ./twin_helices > out 2> log ; cd ..

# helical_frame_mpi:
# 	export OMP_NUM_THREADS=1;
# 	export KMP_NUM_THREADS=1;
# 	qcc -D_MPI=20 -source helical_linearised.c ;
# 	mpiicc $(CFLAGS) _helical_linearised.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o helical_linearised/helical_linearised -lm ;
# 	cd helical_linearised ; mpirun -np 20 -genv I_MPI_PIN_PROCESSOR_LIST=allcores ./helical_linearised > out 2> log ; cd ..

	#icc $(CFLAGS) -prof-gen -prof-dir$(PROF_PATH) _helical.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o helical/helical_icc -lm ;
	#cd helical ; ./helical_icc > out 2> log ; cd ..
	#icc $(CFLAGS) -prof-use -ipo -prof-dir$(PROF_PATH) _helical.c -I$(BASILISK) -L$(BASILISK)/gl $(OPENGLIBS) -o helical/helical_icc -lm ;
	#cd helical ; ./helical_icc > out 2> log ; cd ..
