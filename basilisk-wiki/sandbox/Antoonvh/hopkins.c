/**
This is not basilisk C, but it is part of the methods I use. It relies on some external libraries that are not under basilisk/src.  
A manual will follow some time.  
Database host request; *"Do not use too many parallel threads. Upto 10 orso."*.  
The code is self explainatory known that data is read in XY-slices.
*/

#include <mpi.h> // the sandbox does not complain....
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "turblib.h"


int main(int argc, char *argv[]) {
  MPI_Init(0,0);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Status status1;
  MPI_Status status2;
  MPI_Status status3;
  printf("thread with Rank %d starts working\n",world_rank);

  soapinit();
  turblibSetExitOnError(1);
  
  char * authtoken = "Your_Token_Goes_Here";
  char * dataset = "isotropic1024coarse";
  float time = 0.364;
  int Nx = 1024; // x-width of slice per call
  int Ny = 1024; // y-width of slice per call
  int dbNz = 1024;
  int Zw=1; // slices per call
  int Nz=10; //Total number of slices
  
  //CLA 
  if (argc>1)
    Nz = atoi(argv[1]);
  if (argc>2)
    Zw = atoi(argv[2]);
  if (argc>3)
    Ny = atoi(argv[3]);
  if (argc>4)
    Nx = atoi(argv[4]);
  
  int K=Nz/(Zw*world_size)+0.5; //number of server calls per thread
  if (K*world_size*Zw > dbNz)
    printf("You might obtain zero-padded data\n");
  if (K*world_size*Zw == dbNz)
     printf("Downloading full snapshot\n");
  if (K*world_size*Zw < dbNz)
    printf("Downloading fraction of the database, Nz = %d, K = %d\n",Nz,K);
   
  MPI_File fpu;
  MPI_File_open(MPI_COMM_WORLD,"datau",MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fpu);
  MPI_File fpv;
  MPI_File_open(MPI_COMM_WORLD,"datav",MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fpv);
  MPI_File fpw;
  MPI_File_open(MPI_COMM_WORLD,"dataw",MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fpw);
  
  int components = 3; // Velocity components per cell
  int amount = Nx*Ny*Zw; //Cells per server query 
  
  float * rawdata  = (float*) malloc(amount*sizeof(float)*components);
  float * ud = (float*) malloc(amount*sizeof(float));
  float * vd = (float*) malloc(amount*sizeof(float));
  float * wd = (float*) malloc(amount*sizeof(float));
  if(rawdata == NULL||ud==NULL|| vd==NULL||wd==NULL) {
    printf("malloc failed on thread %d!\n",world_rank);   
    return 1;
  }
  
  for (int it = 0;it<K;it++){
    int Z = world_rank+world_size*it;
    getRawVelocity(authtoken, dataset, time, 0, 0, Z, Nx, Ny, Zw, (char*)rawdata);
    for (int l = 0; l < amount; l++) {
      ud[l]=rawdata[3*l];
      vd[l]=rawdata[3*l+1];
      wd[l]=rawdata[3*l+2];
    }
    // Guide file pointers and write the data
    int g = Z*amount*sizeof(float);
    MPI_Offset offset = g;
    MPI_File_seek(fpu,offset,MPI_SEEK_SET);
    MPI_File_seek(fpv,offset,MPI_SEEK_SET);
    MPI_File_seek(fpw,offset,MPI_SEEK_SET);
    MPI_File_write(fpu,&ud[0],amount,MPI_FLOAT,&status1);
    MPI_File_write(fpv,&vd[0],amount,MPI_FLOAT,&status2);
    MPI_File_write(fpw,&wd[0],amount,MPI_FLOAT,&status3);
    fprintf(stdout,"thread #: %d, Z= %d. So Progress is at %g %%\n",world_rank,Z,100.*(double)(it+1)/(double)K);
  }
  
  // Finish the session
  free(rawdata);
  free(ud);
  free(vd);
  free(wd);
  soapdestroy();
  MPI_Barrier(MPI_COMM_WORLD); // Is this really needed?
  MPI_File_close(&fpu);
  MPI_File_close(&fpv);
  MPI_File_close(&fpw);
  MPI_Finalize();
  return 0;
}
