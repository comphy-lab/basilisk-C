CFLAGS += -I$(DOCUMENT_ROOT)/sandbox
CFLAGS += -disable-dimensions # skip dimensions check for now
#CFLAGS := $(filter-out -g, $(CFLAGS)) # remove debug mode
#CFLAGS += -O2 # compilation optmisation

#CFLAGS += -DTRACE=2 	# toggle tracing
#CFLAGS += -DMTRACE=1 # memory profiling

all: gen_ini.tst ml_breaking_strat.tst
#all: _prep_low_res gen_ini.tst ml_breaking_strat.tst

gen_ini.tst: _prep_ini

# > ppr_C, mpirun
ml_breaking_strat.tst: _prep_local
ml_breaking_strat.tst: CC = mpicc -D_MPI=16

# GPU run
ml_breaking_strat.gpu.c: ml_breaking_strat.c
	ln -sf ml_breaking_strat.c ml_breaking_strat.gpu.c
ml_breaking_strat.gpu.tst: _prep_gpu

# We want to run on HPC
_ml_breaking_strat.c: CFLAGS += -D_MPI=1
hpc_run: _ml_breaking_strat.c
	mkdir -p ml_breaking_strat.hpc
	mv _ml_breaking_strat.c ml_breaking_strat.hpc/
	cp namelist.toml ml_breaking_strat.hpc/
	$(foreach f,$(INI_FILES), cp gen_ini/$(f) ml_breaking_strat.hpc/;)

_prep_ini:
	ln -sf ../namelist.toml gen_ini/namelist.toml
INI_FILES = F_kxky ky kx omega phase

_prep_local:
	mkdir -p ml_breaking_strat
	ln -sf ../namelist.toml ml_breaking_strat/namelist.toml
	$(foreach f,$(INI_FILES), ln -sf ../gen_ini/$(f) ml_breaking_strat/;)

_prep_gpu:
	mkdir -p ml_breaking_strat.gpu
	ln -sf ../namelist.toml ml_breaking_strat.gpu/namelist.toml
	$(foreach f,$(INI_FILES), ln -sf ../gen_ini/$(f) ml_breaking_strat.gpu/;)


include $(BASILISK)/Makefile.defs


















