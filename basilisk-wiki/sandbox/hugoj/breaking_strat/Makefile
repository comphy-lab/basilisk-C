CFLAGS += -I$(DOCUMENT_ROOT)/sandbox
CFLAGS += -O2 -disable-dimensions # skip dimensions check for now

#CFLAGS += -DTRACE=2 	# toggle tracing
#CFLAGS += -DMTRACE=1 # memory profiling

INI_FILES = F_kxky ky kx omega phase


all: gen_ini.tst ml_breaking_strat.gpu.tst
#all: _prep_low_res gen_ini.tst ml_breaking_strat.tst

gen_ini.tst: _prep_ini


# > ppr_C, mpirun
ml_breaking_strat.tst: _prep_local
ml_breaking_strat.tst: CC = mpicc -D_MPI=16

# GPU run
ml_breaking_strat.gpu.tst: _prep_gpu
ml_breaking_strat.gpu.tst: CFLAGS+="-fopenmp"

ml_breaking_strat.cuda.tst: _prep_cuda
#ml_breaking_strat.gpu.tst: CFLAGS+="-fopenmp"

# We want to run on HPC
_ml_breaking_strat.c: CFLAGS +="-D_MPI=1 -DDUMPS"
#_ml_breaking_strat.c: CFLAGS += -grid=cuda/multigrid
hpc_run: _ml_breaking_strat.c
	mkdir -p ml_breaking_strat.hpc
	mv _ml_breaking_strat.c ml_breaking_strat.hpc/
	cp namelist.toml ml_breaking_strat.hpc/namelist.toml
	$(foreach f,$(INI_FILES), cp gen_ini/$(f) ml_breaking_strat.hpc/;)

_prep_ini:
	mkdir -p gen_ini
	ln -sf ../namelist.toml gen_ini/namelist.toml


_prep_local:
	mkdir -p ml_breaking_strat
	ln -sf ../namelist.toml ml_breaking_strat/namelist.toml
	$(foreach f,$(INI_FILES), ln -sf ../gen_ini/$(f) ml_breaking_strat/;)

_prep_gpu:
	mkdir -p ml_breaking_strat.gpu
	ln -sf ../namelist.toml ml_breaking_strat.gpu/namelist.toml
	$(foreach f,$(INI_FILES), ln -sf ../gen_ini/$(f) ml_breaking_strat.gpu/;)

_prep_cuda:
	mkdir -p ml_breaking_strat.cuda
	ln -sf ../namelist.toml ml_breaking_strat.cuda/namelist.toml
	$(foreach f,$(INI_FILES), ln -sf ../gen_ini/$(f) ml_breaking_strat.cuda/;)


include $(BASILISK)/Makefile.defs


















