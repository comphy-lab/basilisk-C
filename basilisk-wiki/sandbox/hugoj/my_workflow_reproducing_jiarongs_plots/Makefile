CFLAGS += -O2 -I$(DOCUMENT_ROOT)/sandbox 

CFLAGS += -disable-dimensions # skip dimensions check for now
LIBS += -lgsl

EXEC=ml_breaking
NAMELIST=namelist


all: gen_ini/${NAMELIST}.toml gen_ini.tst gen_ini/plots


gen_ini/${NAMELIST}.toml: ${NAMELIST}.toml 
	mkdir -p gen_ini
	cp ${NAMELIST}.toml gen_ini/

	









# $(EXEC)/F_kxky $(EXEC).tst 
#
# # initial fields
# F_kxky: gen_ini.c
# 	qcc -autolink ${CFLAGS} -Wall -o gen_ini gen_ini.c -lm
# 	./gen_ini
#
#
#
# $(EXEC)/F_kxky: 
# 	mkdir -p $(EXEC)
# 	mkdir -p $(EXEC)/surface
# 	mkdir -p $(EXEC)/field
# 	cp F_kxky kx ky ./$(EXEC)
#
#
#
# # running in // using mpi
# $(EXEC).tst: CC = mpicc -D_MPI=16
#
# # generate source code for HPC
# _$(EXEC).c: CFLAGS += -D_MPI=1



include $(BASILISK)/Makefile.defs

