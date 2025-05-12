/**

*This is "not a bug" because fields are indeed recognised by their
names in dump files and if you give the same name to different
fields, they will not be distinguished as different in the dump
file. This is intended behaviour.*

The restore() doesn't read more than one scalar data with the same
name. In this sample program, an attribute is defined at the begining,
and scalar fields for an attribute are allocated with the same or
different name. When the same name, qx in this exmaple, is given to
scalar fields using the list iteration, the restore() doesn't
distinguish T.qx, P.qx, and U.qx of dump files. The data of U.qx in a
dump file is restored in the field of T.qx at the restart in some
reason. On the other hand, the restore() can read T.qx, P.qx, and U.qx
separetly, when they are allocated with different names. This issue
didn't happen with the header version of 161020.

For a test, please run a.out from the scratch. If you look into
plot.txt1, it's obvious that T.qx, P.qx, and U.qx store expected
values at this moment. Then, please restart with dump0. You will see
that T.qx, P.qx, and U.qx in plot.txt1 change. */

#include "grid/bitree.h"
#include "run.h"

scalar T[],P[],U[];
scalar * list_Q ={T,P,U};
attribute {
 scalar qx;
}

int main(){
  init_grid(8);
  L0 = 1.0; 
  origin(0,0);
  run(); 
}

event allocation(i=0){
<<<<<<< Inappropriate for header.version 170901
  for (scalar q in list_Q){
    scalar qx = new scalar;
    q.qx = qx;
    boundary({qx});
  }
======
  scalar Tx = new scalar; T.qx = Tx;
  scalar Px = new scalar; P.qx = Px;
  scalar Ux = new scalar; U.qx = Ux;
  boundary({Tx,Px,Ux});
>>>>>>> Compatible with header.version 170901
}

event init(i=0){
  dt =1e-6;
  if(!restore(file="dump")){
      printf("No dump file.\n");
      scalar Tx = T.qx;
      scalar Px = P.qx;
      scalar Ux = U.qx;
      foreach(){
         Tx[]=100.; 
         Px[]=10.; 
         Ux[]=1.; 
      }
  }
}

event output(i=0;i+=1;i<=1){
  char filename[20]="plot.txt",num[5];
  sprintf(num,"%d",i);
  strcat(filename,num);
  FILE * fpx = fopen(filename, "w");
  scalar Tx = T.qx;
  scalar Px = P.qx;
  scalar Ux = U.qx;
  foreach(){
    fprintf (fpx,"%5.3e %5.3e %5.3e %5.3e\n",
                   x,  Tx[],Px[],Ux[]);
  }
  fclose(fpx);
}

event dumpout(i=0;i+=1;i<=2){
  char dumpname[20]="dump",num[100];
  sprintf(num,"%d",i);
  strcat(dumpname,num);
  dump(file=dumpname);
}


