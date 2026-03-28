CFLAGS += -I$(DOCUMENT_ROOT)/sandbox #-O2
CFLAGS += -disable-dimensions # skip dimensions check for now
#LIBS += -lgsl

EXEC=ml_breaking
NAMELIST=namelist
INI=gen_ini

#CFLAGS += -DUSE_PPR=1

# all: initial_condition sim
#
# sim: ${EXEC}/${NAMELIST}.toml ${EXEC}.tst ${EXEC}/plots
#
# sim.gpu: ${EXEC}.gpu/${NAMELIST}.toml ${EXEC}.gpu.tst ${EXEC}.gpu/plots
#
# ${EXEC}.tst: CC = mpicc -D_MPI=16
#
# initial_condition: ${INI}/${NAMELIST}.toml ${INI}.tst ${INI}/plots
#
# namelist: ${INI}/${NAMELIST}.toml ${EXEC}/${NAMELIST}.toml
#
#
# ${INI}/${NAMELIST}.toml: ${NAMELIST}.toml 
# 	mkdir -p ${INI}
# 	ln -sf ../${NAMELIST}.toml ${INI}/
#
# ${EXEC}/${NAMELIST}.toml: ${NAMELIST}.toml 
# 	mkdir -p ${EXEC}
# 	ln -sf ../${NAMELIST}.toml ${EXEC}/
# 	ln -sf ../${INI}/F_kxky ${EXEC}/
# 	ln -sf ../${INI}/ky ${EXEC}/
# 	ln -sf ../${INI}/kx ${EXEC}/
# 	ln -sf ../${INI}/omega ${EXEC}/
# 	ln -sf ../${INI}/phase ${EXEC}/
#
# ${EXEC}.gpu/${NAMELIST}.toml: ${NAMELIST}.toml 
# 	mkdir -p ${EXEC}.gpu
# 	ln -sf ../${NAMELIST}.toml ${EXEC}.gpu/
# 	ln -sf ../${INI}/F_kxky ${EXEC}.gpu/
# 	ln -sf ../${INI}/ky ${EXEC}.gpu/
# 	ln -sf ../${INI}/kx ${EXEC}.gpu/
# 	ln -sf ../${INI}/omega ${EXEC}.gpu/
# 	ln -sf ../${INI}/phase ${EXEC}.gpu/
#
# _${EXEC}.c: CFLAGS += -D_MPI=1
#
#

INI_FILES = F_kxky ky kx omega phase
define setup_exec_dir
	mkdir -p $(@D)
	cp $(NAMELIST).toml $(@D)/
	$(foreach f,$(INI_FILES), ln -sf ../$(INI)/$(f) $(@D)/;)
endef

# General rule
all: ${EXEC}.tst

# We want to run on local machine

# > link ini files 
${EXEC}/${NAMELIST}.toml ${EXEC}_pprf90/${NAMELIST}.toml ${EXEC}.gpu/${NAMELIST}.toml: ${NAMELIST}.toml
	$(setup_exec_dir)

# > ppr_C, mpirun
${EXEC}.tst: _prep_lowres
${EXEX}.tst: ${EXEC}/${NAMELIST}.toml
${EXEC}.tst: CC = mpicc -D_MPI=16

# > ppr_f90, mpirun
${EXEC}_pprf90.c: ${EXEC}.c
	ln -sf ${EXEC}.c ${EXEC}_pprf90.c
${EXEC}_pprf90.tst: _prep_lowres
${EXEC}_pprf90.tst: ${EXEC}_pprf90/${NAMELIST}.toml
${EXEC}_pprf90.tst: CC = mpicc -D_MPI=16 -DUSE_PPR=1


# > gpu (ppr_C)
# dont forget the env variables to select nvidia gpu !
# __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia make ml_breaking.gpu.tst
${EXEC}.gpu.c: ${EXEC}.c
	ln -sf ${EXEC}.c ${EXEC}.gpu.c
${EXEC}.gpu.tst: _prep_midres
${EXEC}.gpu.tst: ${EXEC}.gpu/${NAMELIST}.toml

# We want to run on HPC
_${EXEC}.c: CFLAGS += -D_MPI=1
_${EXEC}.c: _prep_highres

hpc_run: _${EXEC}.c
	mkdir -p $@
	mv _${EXEC}.c $@/
	mv ${NAMELIST}.toml $@/


.PHONY: _prep_highres _prep_lowres, _prep_midres
_prep_highres:
	cp ${NAMELIST}_highres.toml ${NAMELIST}.toml

_prep_midres:
	cp ${NAMELIST}_midres.toml ${NAMELIST}.toml


_prep_lowres:
	cp ${NAMELIST}_lowres.toml ${NAMELIST}.toml



include $(BASILISK)/Makefile.defs

