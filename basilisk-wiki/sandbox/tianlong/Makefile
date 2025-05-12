# The makefile for compiling the code.
# To use MPI, add 'mpi=1'.
# For non-cubic domains, add 'solid=1'.
# For example, to compile the bubble rising case, use the command:
# make pcrising mpi=1 solid=1


TARGET = semupc

SEMU_DIR = semushin

MCC = qcc
MACRO := -I$(SEMU_DIR)

ifeq ($(dbg), 1)
	MACRO +=
	FLAGS += -g  -lm
else
	FLAGS += -O2 -lm
endif

#For MPI
ifeq ($(mpi), 1)
	MCC = CC99='mpicc -std=c99' qcc
	MACRO += -D_MPI=1
endif

ifeq ($(solid),1)
MACRO += -DUSE_MY_SOLID=1
endif

FLAGS += -autolink

SEMUSHIN_H = $(SEMU_DIR)/ebit_utils.h $(SEMU_DIR)/fitting.h $(SEMU_DIR)/semushin.h $(SEMU_DIR)/semushin_two-phase.h $(SEMU_DIR)/advection.h

all: data $(SEMUSHIN_H)
	$(MCC) $(MACRO) main.c -o $(TARGET) $(FLAGS)

stefan:data $(SEMUSHIN_H)
	$(MCC) $(MACRO) stefanproblem.c -o $(TARGET) $(FLAGS)

sucking:data $(SEMUSHIN_H)
	$(MCC) $(MACRO) suckingproblem.c -o $(TARGET) $(FLAGS)

scriven:data $(SEMUSHIN_H)
	$(MCC) $(MACRO) scrivenproblem.c -o $(TARGET) $(FLAGS) -lgsl -lgslcblas

pcrising:data $(SEMUSHIN_H)
	$(MCC) $(MACRO) bubblerising.c -o $(TARGET) $(FLAGS) -lgsl -lgslcblas

fixedflux:data $(SEMUSHIN_H)
	$(MCC) $(MACRO) fixedflux.c -o $(TARGET) $(FLAGS)

filmboiling:data $(SEMUSHIN_H)
	$(MCC) $(MACRO) filmboiling.c -o $(TARGET) $(FLAGS)

data:
	@mkdir -p data

cleandata:
	rm -f data/*

clean:
	rm -rf *.o $(TARGET) data 
