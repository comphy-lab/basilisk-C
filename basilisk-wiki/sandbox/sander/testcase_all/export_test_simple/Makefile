INCL=-I../../
all: export_test_simple export_test_simple_mpi

export_test_simple_mpi: export_test_simple.c
	qcc -O2 -Wall $(INCL) $(CONF) -D_MPI=1 -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -o $@ $< -lm -lmpi
export_test_simple: export_test_simple.c
	qcc -O2 -Wall $(INCL) $(CONF) -o $@ $<  -lm

clean:
	rm -rf htg/
	rm -rf vti/
	rm -rf vtp/
	rm -f export_test_simple
	rm -f export_test_simple_mpi
	rm -rf header/
	rm -f *.pvd
	rm -f perfs
