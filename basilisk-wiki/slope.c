#include "grid/cartesian1D.h"
scalar dth[];
#include "saint-venant.h"
double a,tmax;
scalar d2h[],d3h[];
scalar ho[];
u.n[left] = neumann(0);
h[left] = h[];
event init (i = 0)
{
  foreach(){
    zb[] =   (x>10)*(x-10)/(25);
    u.x[]=a*exp(-(x+12)*(x+12)) ;
    h[]=fmax(1+u.x[]-zb[],0);
    }
}
event field (t<tmax;t+=1) {
    printf("p[-15:][-.1:1.1]'-' u 1:3 t'u'w l,''u 1:4 t'z' w l,''u 1:($2+$4) t'eta' w l\n");
    foreach()
    printf (" %g %g %g %g %g \n", x, h[], u.x[], zb[], t);
    printf ("e\n\n");
}

event end (t = tmax) {
    printf("p[-15:][-.1:1.1]sin(x)\n\n\n");
     
}
int main()
   {  
    X0 = -15.;
    L0 = 60.;
    G = 1;
    N = 1024; 
    a=0.01;
    tmax=150;	
    run();
    }
