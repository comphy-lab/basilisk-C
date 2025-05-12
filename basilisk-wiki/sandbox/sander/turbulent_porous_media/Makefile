#CC = CC99='mpicc -std=c99'
CC = 
#S=-source
#S=
CFLAGS = -Wall -O2 -I../
LIB = -lm
TARGET=turbulent_flow_porous_media_sq_cylinder.c
MPI = -D_MPI=1 -I/usr/lib/x86_64-linux-gnu/openmpi/include/ 
LIB_MPI = -lmpi

all: turbulent_flow_porous_media_sq_cylinder

turbulent_flow_porous_media_sq_cylinder: $(TARGET)
	$(CC) qcc $S $(CFLAGS) $(MPI) -o $@ $< $(LIB) $(LIB_MPI)

clean:
	rm -rf htg
	rm -f turbulent_flow_porous_media_sq_cylinder perfs
	rm -f *.pvd
