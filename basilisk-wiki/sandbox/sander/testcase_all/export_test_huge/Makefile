CONF=
INCL=-I../../
all: export_test_huge export_test_huge_mpi

export_test_huge_mpi: export_test_huge.c
	CC99='mpicc -std=c99' qcc -O2 -Wall $(INCL) $(CONF) -D_MPI=1 -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -o $@ $< -lm -lmpi
export_test_huge: export_test_huge.c
	qcc -O2 -Wall $(INCL) $(CONF) -o $@ $<  -lm

clean:
	rm -rf htg/
	rm -rf vti/
	rm -rf vtp/
	rm -f export_test_huge
	rm -f export_test_huge_mpi
	rm -rf header/
	rm -f *.pvd
	rm -f perfs
