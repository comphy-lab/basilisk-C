CFLAGS += -I$(DOCUMENT_ROOT)/sandbox #-O2
CFLAGS += -disable-dimensions # skip dimensions check for now
#LIBS += -lgsl

EXEC=ml_breaking
NAMELIST=namelist
INI=gen_ini

CFLAGS += -DUSE_PPR=1

all: initial_condition 

sim: ${EXEC}/${NAMELIST}.toml ${EXEC}.tst ${EXEC}/plots

sim.gpu: ${EXEC}.gpu/${NAMELIST}.toml ${EXEC}.gpu.tst ${EXEC}.gpu/plots

${EXEC}.tst: CC = mpicc -D_MPI=16

initial_condition: ${INI}/${NAMELIST}.toml ${INI}.tst ${INI}/plots

namelist: ${INI}/${NAMELIST}.toml ${EXEC}/${NAMELIST}.toml


${INI}/${NAMELIST}.toml: ${NAMELIST}.toml 
	mkdir -p ${INI}
	ln -sf ../${NAMELIST}.toml ${INI}/

${EXEC}/${NAMELIST}.toml: ${NAMELIST}.toml 
	mkdir -p ${EXEC}
	ln -sf ../${NAMELIST}.toml ${EXEC}/
	ln -sf ../${INI}/F_kxky ${EXEC}/
	ln -sf ../${INI}/ky ${EXEC}/
	ln -sf ../${INI}/kx ${EXEC}/
	ln -sf ../${INI}/omega ${EXEC}/
	ln -sf ../${INI}/phase ${EXEC}/

${EXEC}.gpu/${NAMELIST}.toml: ${NAMELIST}.toml 
	mkdir -p ${EXEC}.gpu
	ln -sf ../${NAMELIST}.toml ${EXEC}.gpu/
	ln -sf ../${INI}/F_kxky ${EXEC}.gpu/
	ln -sf ../${INI}/ky ${EXEC}.gpu/
	ln -sf ../${INI}/kx ${EXEC}.gpu/
	ln -sf ../${INI}/omega ${EXEC}.gpu/
	ln -sf ../${INI}/phase ${EXEC}.gpu/

_${EXEC}.c: CFLAGS += -D_MPI=1







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

