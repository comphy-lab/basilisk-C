exec=main

CC=$(BASILISK)/qcc
# CFLAGS="-autolink -disable-dimensions -grid=multigrid -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp"
CFLAGS= -autolink -grid=multigrid -g -Wall -pipe -D_FORTIFY_SOURCE=2 -O2 -fopenmp
LIBGL= -L$(BASILISK)/gl -lglutils
INCLUDE= $(DOCUMENT_ROOT)/sandbox # "$(SANDBOX)"
OPENGLIBS= -lfb_tiny
MATHLIB= -lm
DEPS= 
SRC := $(exec).c 


$(exec): $(SRC) $(DEPS)
	mkdir -p $(exec)
	cp params.in $(exec)/params.in
	qcc -I$(INCLUDE) $(CFLAGS) $(EVENTS) -o $(exec)/$(exec) $(SRC) $(LIBGL) $(OPENGLIBS) $(MATHLIB) 





