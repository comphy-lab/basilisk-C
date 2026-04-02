CFLAGS += -I$(DOCUMENT_ROOT)/sandbox
CFLAGS += -disable-dimensions # skip dimensions check for now
CFLAGS := $(filter-out -g, $(CFLAGS)) # remove debug mode
CFLAGS += -O2 # compilation optmisation

CFLAGS += -DTRACE=2 	# toggle tracing
#CFLAGS += -DMTRACE=1 # memory profiling

EXEC=ml_convection
NAMELIST=namelist
INI=gen_ini

define setup_exec_dir
	mkdir -p $(@D)
	cp $(NAMELIST).toml $(@D)/
endef

# General rule
all: ${EXEC}.tst

# We want to run on local machine

# > link ini files 
${EXEC}/${NAMELIST}.toml ${EXEC}.gpu/${NAMELIST}.toml: ${NAMELIST}.toml
	$(setup_exec_dir)

# > ppr_C, mpirun
${EXEC}.tst: _prep_lowres
${EXEC}.tst: ${EXEC}/${NAMELIST}.toml
#${EXEC}.tst: CC = mpicc -D_MPI=16

# > gpu (ppr_C)
# dont forget the env variables to select nvidia gpu !
# __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia make ml_convection.gpu.tst
# nohup __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia make ml_convection.gpu.tst > nohup.out 2>&1 &
${EXEC}.gpu.c: ${EXEC}.c
	ln -sf ${EXEC}.c ${EXEC}.gpu.c
${EXEC}.gpu.tst: _prep_midres
${EXEC}.gpu.tst: ${EXEC}.gpu/${NAMELIST}.toml

# We want to run on HPC
_${EXEC}.c: CFLAGS += -D_MPI=1
_${EXEC}.c: _prep_highres

hpc: _${EXEC}.c
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

