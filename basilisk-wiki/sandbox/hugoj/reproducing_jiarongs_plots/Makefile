# #SET FILES OF INTEREST 
# SRC_FILES = $(addsuffix .c, $(OUTPUT))
# MPI_FILES = $(addprefix _, $(SRC_FILES))
# OBJ_FILES = $(addsuffix .o,$(basename $(SRC_FILES)))
# OUT_FOLDER = .
#
# #SET C FLAGS AND BASILISK COMPILER, LINKING FLAGS
# BAS=qcc 
# CFLAGS += -disable-dimensions -autolink -D_MPI=1 -DPARAVIEW=0 -source
# MPICC = mpicc
# MPICCFLAGS += -std=c99 -O2 -g -Wall -D_MPI=1
# LDFLAGS = -lgfortran -L${BASILISK}/ppr -lppr -lm
# OPENGLINC = -I${HOME}/local/include
# OPENGLIBS = -L${BASILISK}/gl -lfb_tiny
# #OPENGLIBS = -L${BASILISK}/gl -lglutils -lfb_osmesa -L${HOME}/local/lib -lGLU -L${HOME}/local/lib -lOSMesa 
#
# #MAKE COMMAND: TAKES ARGUMENT AND MAKES EXECUTABLE IN DESIRED FOLDER
# $(OUTPUT):
# #	$(BAS) $(CFLAGS) $(SRC_FILES)
# #	$(MPICC) $(MPICCFLAGS) $(MPI_FILES) -o $(OUT_FOLDER)/$(OUTPUT) -I$(BASILISK) ${OPENGLINC} ${OPENGLIBS} $(LDFLAGS)
	#CC99='mpicc -std=c99' qcc -disable-dimensions -autolink -D_MPI=1 -DPARAVIEW=0 -O2 -g -Wall -o $(OUT_FOLDER)/$(OUTPUT) $(SRC_FILES) -I$(BASILISK) -I${HOME}/local/include -L${BASILISK}/gl -lfb_tiny -lgfortran -L${BASILISK}/ppr -lppr -lm

CFLAGS += -O2 -I$(DOCUMENT_ROOT)/sandbox 
CFLAGS += -disable-dimensions # skip dimensions check for now

field_PM.tst: field_PMmpi.tst

# Jiarong's way of compiling
# field_PM/field_PM: field_PM.c
# 	mkdir -p field_PM
# 	qcc $(CFLAGS) -autolink -DPARAVIEW=0 -O2 -g -Wall -o field_PM/field_PM field_PM.c -I$(BASILISK) -L${BASILISK}/gl -lfb_tiny -lgfortran -L${BASILISK}/ppr -lppr -lm
#
# field_PMmpi/field_PM: field_PM.c
# 	mkdir -p field_PMmpi
# 	CC99='mpicc -std=c99' qcc $(CFLAGS)-autolink -DPARAVIEW=0 -O2 -g -Wall -o field_PM/field_PM field_PM.c -I$(BASILISK) -L${BASILISK}/gl -lfb_tiny -lgfortran -L${BASILISK}/ppr -lppr -lm


field_PMmpi.c: field_PM.c
	ln -sf field_PM.c field_PMmpi.c

field_PMmpi.tst: CC = mpicc -D_MPI=16


include $(BASILISK)/Makefile.defs


