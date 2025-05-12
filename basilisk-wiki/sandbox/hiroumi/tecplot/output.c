/**
This code demonstrates the output of selected scalars and vectors that are defined at cell center in Tecplot binary format. The compilation should be done with linking a TecIO library (tecio) to the code.

For instance, 

in case with openmp : qcc -O2 -Wall -fopenmp output.c -I/(Path to the directory of tecio source files, i.e. teciosrc) -L/(Path to the directory of libtecio.a, e.g. teciosrc) -ltecio -lstdc++ -lpthread -lm

in case with MPI : CC99='mpicc -std=c99' qcc -O2 -Wall -D_MPI=3 output.c -I/(Path to the directory of tecio source files, i.e. teciosrc) -L/(Path to the directory of libtecio.a, e.g. teciosrc) -ltecio -lstdc++ -lpthread  -lm

The TecIO library can be download from a Tecplot company website for free.
*/

#include "grid/octree.h"
#include "run.h"
#if dimension>1
#include "tecplot.h"
#endif

//Primary variables
  scalar T[], P[];
  vector vel[];
  scalar * list_Q = {T,P,vel};
//Dependent variables
  scalar RHO[];
//New boundary
  bid ball;

int main(){
  N=64;
  L0 = 1.0;
  origin(-L0/2,-L0/2.,-L0/2.);
  dt = 0.25;
  run();
}

event init(i=0){
#if TREE
   refine((sq(x)+sq(y)+sq(z)<sq(0.25))&& level<8 );
#ifndef _MPI
   mask (sq(x)+sq(y)+sq(z)<sq(0.2) ? ball:
        none);
#endif
#endif
  foreach(){
    double rad = sqrt(sq(x)+sq(y)+sq(z));
    T[]=296.0+20*rad;
    P[]=1.01325e+5;
    vel.x[]=10.0*sin(2*rad*pi);
    vel.y[]=10.0*cos(2*rad*pi);
#if dimension>2
    vel.z[]=10.0*tan(2*rad*pi);
#endif
    RHO[] = P[]/(296.0*T[]);
  }
  boundary(list_Q);
  boundary({RHO});
}

event step(i++){
  dt = dtnext(dt);
}

event output(i=0;i+=10;i<=100){
  struct OutputTec tec;
//Give primary variables to the tec_cc list. "cc" means "cell center".
  tec.tec_cc = list_copy(list_Q);
//Give the name of variables to varname.
#if dimension>1
  sprintf(tec.varname,"x y T P vel.x vel.y");
#endif
#if dimension>2
  sprintf(tec.varname,"x y z T P vel.x vel.y vel.z");
#endif

//Extra variables can be output by adding them to the list of scalar and name.
  tec.tec_cc = list_concat(tec.tec_cc,{RHO});
  sprintf(tec.extname," RHO");

  output_tec(tec,i);
  free(tec.tec_cc);
}

