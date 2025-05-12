/**
# kind of hydrodynamic internal dam break
 
 `to be cleand internaldam_2`
 
 Two layer Shallow Water equations.
The layer 1 is the upper layer, the layer 2 is the lower layer. The bottom is flat here. The mass conservation for each layer (no exchange of mass)
 $$
 \frac{\partial}{\partial t}(h_1) +\frac{\partial}{\partial x}(h_1u_1 )=0$$
 $$
 \frac{\partial}{\partial t}(h_2) +\frac{\partial}{\partial x}(h_2u_2 )= 0$$
Momentum:
$$
 \frac{\partial}{\partial t}(h_1u_1) +\frac{\partial}{\partial x}(h_1u_1^2+\frac{gh_1^2}{2})= - g h_1 \frac{\partial}{\partial x}h_2$$
 $$
 \frac{\partial}{\partial t}(h_2u_2) +\frac{\partial}{\partial x}(h_2u_2^2+\frac{gh_2^2}{2})= - g \frac{\rho_1}{\rho_2} h_2 \frac{\partial}{\partial x}h_1$$
 
note that the ratio $\frac{\rho_1}{\rho_2}$ is lesser than one
 
 */
#include "grid/multigrid.h"
/**
 Using Emily two layer implementation
 */
#include "two_layer_emily.h"
double tmax,x0;
scalar dh[];

int main(){
    RHO1 = .98;   //layer1 over layer2
    RHO2 = 1.;	  //lower layer
    X0 = -5.;
    L0 = 10.;
    N = 128*2;
    x0 = 0;
    G = 1;
    tmax = 1;
#define DA
    RHO1 =0.7;
    DT = HUGE;
    tmax=1*sqrt(10);
    run();
}
/** 
 simple infinite domain BC
*/
h1[left] = neumann(0);
h1[right] = neumann(0);
u1.n[left] = neumann(0);
u1.n[right] = neumann(0) ;
h2[left] = neumann(0);
h2[right] = neumann(0);
u2.n[left] = neumann(0);
u2.n[right] = neumann(0) ;
/**
 first case, the mooving internal jump with no perturbation of the total interface
 initial case taken from Abgrall (and Pares)
 */
#ifdef RE
event init (i = 0)
{
    foreach(){
        zb[]=0.;
        h2[]=0.5*(x<x0) + 0.45*(x>=x0);
        u2.x[]=(1.25*(x<x0) + 1.125*(x>=x0))/h2[];
        h1[]=1-h2[]-zb[];
        u1.x[]=(1.25*(x<x0) + 1.375*(x>=x0))/h1[];
    }
    boundary({zb,h1,u1,h2,u2});
}
#endif
/**
 second case, the internal dam break
 initial case taken from Berthon et al.
*/
#ifdef DA
event init (i = 0)
{
    foreach(){
        zb[]=0.;
        h2[]=0.2*(x<x0) + 1.8*(x>=x0);
        u2.x[]=0;
        h1[]=2-h2[]-zb[];
        u1.x[]=0;
    }
    boundary({zb,h1,u1,h2,u2});
}
#endif
/**
 save profiles along a line
*/
event printprofile  (t<=tmax;t+=0.1){
    char name[100];
    FILE * fp;
    sprintf(name,"profil.OUT");
    fp=fopen(name,"w");
    double dx=L0/N;
    for(double x=X0;  x<X0+L0;x+=dx)
    {
        fprintf (fp,"%g %g %g %g %g \n",
                 x,interpolate (h1, x, 0),interpolate (h2, x, 0),interpolate (u1.x, x, 0),interpolate (u2.x, x, 0));}
    // foreach()
    // fprintf(fp,"%g %g %g %g %g\n",x,h1[],h2[],u1.x[],u2.x[]);
    fclose(fp);
}
/**
 plot on the fly
*/
event plot (t<=tmax;t+=0.1) {
    printf("set title 't= %.2lf '\n",t);
    printf("p[%lf:%lf][-.5:2.5]  '-' u 1:($2+$3) t'free surface' w l lt 3,'-' u 1:3 t'internal' w l linec 5,'-' u 1:4 t'zb' w l lt -1 \n",
           X0,X0+L0);
    foreach()
    printf (" %g %g %g %g  \n", x, h1[], h2[],zb[]);
    printf ("e \n");
}
/**
 
# Run
 
 compilation and execution
 
~~~
 qcc internaldam.c -o internaldam -lm
 ./internaldam | gnuplot
~~~
 
# Results
 
 The heavy thick layer falls on it self to the left, it creates a surface wave going to the left and a depression wave goieng to the right
 
 
~~~gnuplot plot of profiles
 
 p [:][-.5:3]'profil.OUT' u 1:($2+$3+.5) t'surf' w l,'profil.OUT'u 1:(($3+.5)) t 'int.' w l,\
 ''u 1:4 t'u surf' w l,''u 1:5 t'u int' w l,0 not lc-1
 
~~~
 
 We superpose to Berthon et al. simulation
 
~~~gnuplot (Comparison with Figure 6.3 from Berthon)
 
 set output 'fig1.png'
 X0 = 103
 X1 = 758
 Y0 = 128
 Y1 = 604
 unset tics
 plot [:][:] '../Img/berthon.png' binary filetype=png with rgbimage not,\
 'profil.OUT' u (X0+($1+5)/10*(X1-X0)):(($2+$3+.5)/3*(Y1-Y0)+Y0)t'surf' w l,\
 'profil.OUT'u (X0+($1+5)/10*(X1-X0)):(($3+.5)/3*(Y1-Y0)+Y0) t 'int' w l lc 1
 
~~~
 
 and observe the early stage of the generated wave: slow wave of media 2, fast wave of media 1.
 

 
# Bibliography
 
 *  Abgrall and Karni [The two-layer shallow water ...](https://www.math.u-bordeaux.fr/~rabgrall/mes_papiers/SISC_31_3_2009.pdf)
 
 * Berthon Foucher Morales [.. for two-layer shallow water...](http://www.math.sciences.univ-nantes.fr/~berthon/articles/two_layers.pdf)
 
 */
