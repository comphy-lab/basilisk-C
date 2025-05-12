#if _MPI

typedef struct {
  double x;
  double y;
  double z;
  double val;
} particle;

#define MPI_DIST_INIT \
   int out_mpi = 0;\
   particle *mpi_array = malloc(sizeof(particle) * mall_size);
 
#define MPI_DIST_(scalar, val, SUM, err) \
       if (is_LOCAL_and_ACTIVE) {\
           if (SUM){\
              scalar[] += val;\
           }else{\
              scalar[] = val;\
           }\
       } else {\
           if (out_mpi >= mall_size-1) { \
               fprintf(stdout,"out_mpi (near_interface) %d - increase mall_size %d (mpi_distribute.h)\n", out_mpi, mall_size);\
               fprintf(stderr,"out_mpi (near_interface) %d - increase mall_size %d (mpi_distribute.h)\n", out_mpi, mall_size);\
               STOP(err);\
           }\
           particle pp = {x, y, z, val};\
           mpi_array[out_mpi] = pp;\
           out_mpi++;\
       }
 
#define MPI_DIST_FINISH_(scalar, SUM) \
   mpi_boundary_distribute(scalar, mpi_array, out_mpi, SUM);\
   out_mpi = 0;\
   free(mpi_array);

#define MPI_DIST_FINISH(scalar) MPI_DIST_FINISH_(scalar, false)
#define MPI_DIST(scalar, val, err) MPI_DIST_(scalar, val, false, err)
// sumation versions
#define MPI_DIST_FINISH_SUM(scalar) MPI_DIST_FINISH_(scalar, true)
#define MPI_DIST_SUM(scalar, val, err) MPI_DIST_(scalar, val, true, err)


static void set_cell_func(Point point, particle P, scalar P_scalar, bool SUM) {
  if (SUM) {  // for m_dot
    P_scalar[] += P.val;
  } else {  // for interface distance
    if (fabs(P_scalar[]) < 1.0e-10) {
      // no closer value
      P_scalar[] = P.val;
    } else {
      if (fabs(P_scalar[]) > fabs(P.val)) {
        // only update if closer
        P_scalar[] = P.val;
      }
    }
  }
  return;
}

void mpi_boundary_distribute(scalar P_scalar, particle *send_array, int out,
                             bool SUM) {
  particle *senddata = malloc(sizeof(particle) * out);
  for (int i = 0; i < out; i++) {
    senddata[i] = send_array[i];
  }
  // MPI DISTRIBUTE
  /// Gather lost particles among threads:
  // First, count all of them
  int outa[npe()];
  int outat[npe()];
  outat[0] = 0;

  MPI_Allgather(&out, 1, MPI_INT, &outa[0], 1, MPI_INT, MPI_COMM_WORLD);
  //
  // Compute displacements
  for (int j = 1; j < npe(); j++) {
    outat[j] = outa[j - 1] + outat[j - 1];
  }
  int outt = outat[npe() - 1] + outa[npe() - 1];

  // Allocate recieve buffer and gather
  particle *recdata = malloc(sizeof(particle) * outt);
  for (int j = 0; j < npe(); j++) {
    outat[j] *= sizeof(particle);
    outa[j] *= sizeof(particle);
  }
  // send and recieve data
  MPI_Allgatherv(&senddata[0], outa[pid()], MPI_BYTE, &recdata[0], outa,
                 outat, MPI_BYTE, MPI_COMM_WORLD);

  int in = 0;
  for (int i = 0; i < outt; i++) {
    particle P = recdata[i];
    if (locate(P.x, P.y, P.z).level > 0) {
      in++;
    }
  }
  // Collect new particles from `recdata`
  if (in > 0) {
    int indi[in];
    int m = 0;
    for (int i = 0; i < outt; i++) {
      particle P = recdata[i];
      if (locate(P.x, P.y, P.z).level > 0) {
        indi[m] = i;
        m++;
      }
    }
    for (int i = 0; i < in; i++) {
      int index = indi[i];
      particle P = recdata[index];
      // locate point
      Point point = locate(P.x, P.y, P.z);
      set_cell_func(point, P, P_scalar, SUM);
    }
  }
  
  // clean the mess
  free(senddata);
  free(recdata);

  return;
}
#endif
