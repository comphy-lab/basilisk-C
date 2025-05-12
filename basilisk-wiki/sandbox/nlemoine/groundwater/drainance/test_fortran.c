#include "nlemoine/groundwater/drainance/libexpokit.h"
#include "run.h"
#pragma autolink -L. -lexpokit -lgfortran

int main (int argc, char * argv[])
{
  run();
  return(0);
}

event init (i=0)
{
   int m = 3, ldh = 3;
   int iflag;
/*   double H[3][3] = {{-1.0/86400.,0.,0.},
                     {0.2/86400./0.05,-(1.0e-7)/0.05,+(1.0e-7)/0.05},
                     {0.,+(1.0e-7)/0.001,-(1.0e-7)/0.001}};*/
   double H[3][3] = {{-1.0/86400.,0.2/86400./0.05,0.},
                     {0.,-(1.0e-7)/0.05,+(1.0e-7)/0.001},
                     {0.,+(1.0e-7)/0.05,-(1.0e-7)/0.001}};

   double Dtexp = 3600.;
   double yvec[3] = {0.0,0.0,1.0};
   double wsp[30]; // 2*m*(m+2)
   int iwsp[3];

/*                        A3 = [ -1/tauperc,0,0 ; ..
                               [+wsoil/tauperc,-Cd,+Cd]/wsapro ; ..
                               [0,+Cd,-Cd]/Sfrac ];

Fortran arrays are stored in column-major order: A(3,2)
A(1,1) A(2,1) A(3,1) A(1,2) A(2,2) A(3,2)

C arrays are stored in row-major order: A[3][2]
A[0][0] A[0][1] A[1][0] A[1][1] A[2][0] A[2][1]

      integer          m, ldh, iflag, iwsp(m)
      double precision t, H(ldh,m), y(m)
      complex*16       wsp(m*(m+2))
*/

   // Test de l'exponentiation matricielle

   (void) dgchbv_ ( &m, &Dtexp, &H[0][0], &ldh, &yvec[0], &wsp[0], &iwsp[0], &iflag );

   fprintf(stderr,"iflag = %d\n",iflag);
   fprintf(stderr,"yvec = [%.4f,%.4f,%.4f]\n",yvec[0],yvec[1],yvec[2]);

}