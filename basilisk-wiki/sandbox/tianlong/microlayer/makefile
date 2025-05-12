#make tempinit mpi=1 solid=1 to simulate nucleation process
#make init solid=1 to generate an initial temperature field to procedd bubble growth
#make expbucci mpi=1 solid=1 to simulate the case in (bures and sato, JFM, 2022)
#make longexpbucci to simulate the whole bubble cycle

TARGET = aboiling

MCC = qcc
MACRO :=

ifeq ($(dbg), 1)
	MACRO += -DDEBUGTL=1
	FLAGS += -g  -lm
else
	FLAGS += -O2 -lm
endif

#For MPI
ifeq ($(mpi), 1)
	MCC = CC99='mpicc -std=c99' qcc
	MACRO += -D_MPI=1
endif

ifeq ($(omp),1)
	MACRO += -DUSEOPENMP
	FLAGS = -O2 -lm -fopenmp
endif

ifeq ($(solid),1)
MACRO += -DUSE_MY_SOLID=1
endif

FLAGS += -autolink -Isrc_pc


longexpbucci: data
	$(MCC) $(MACRO) expbucci_longtime.c -o $(TARGET) $(FLAGS)

expbucci: data
	$(MCC) $(MACRO) expbucci.c -o $(TARGET) $(FLAGS)

tempinit:data
	$(MCC) $(MACRO) main_temp_init.c -o $(TARGET) $(FLAGS)

microinit:data
	$(MCC) $(MACRO) main_init.c -o $(TARGET) $(FLAGS)


data:
	@mkdir data

cleandata:
	rm -f data/*

clean:
	rm -rf *.o $(TARGET) data 
