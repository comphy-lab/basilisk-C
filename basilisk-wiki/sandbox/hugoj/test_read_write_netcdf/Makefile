# exec=breaking
#
# CC=$(BASILISK)/qcc
# # CFLAGS="-autolink -disable-dimensions -grid=multigrid -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp"
# CFLAGS= -autolink -grid=multigrid -g -Wall -pipe -O2 -D_FORTIFY_SOURCE=2 -fopenmp
# LIBGL= -L$(BASILISK)/gl -lglutils
# INCLUDE= $(DOCUMENT_ROOT)/sandbox
# LIBS= -lm -lfb_tiny -lnetcdf
# DEPS= 
# SRC := $(exec).c 
#


# $(exec): $(SRC) $(DEPS)
# 	mkdir -p $(exec)
# 	qcc -I$(INCLUDE) $(CFLAGS) $(EVENTS) -o $(exec)/$(exec) $(SRC) $(LIBGL) $(LIBS)
#
#

CFLAGS += -O2 -I$(DOCUMENT_ROOT)/sandbox 
CFLAGS += -disable-dimensions # skip dimensions check for now

all: breaking.tst

breaking.tst: CC = mpicc -D_MPI=16

clean:
	rm $(exec)/$(exec)


include $(BASILISK)/Makefile.defs
