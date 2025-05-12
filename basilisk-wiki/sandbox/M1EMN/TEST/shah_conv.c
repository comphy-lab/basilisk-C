/**
# Avalanche of two newtonian fluids of different density and viscosity in 1D

To address the issue of shock formation in transport equations, we propose an improved approach as a preferable alternative to the one used in [one dimensional two layer avalanche](http://basilisk.fr/_edit/sandbox/M1EMN/TEST/shah_flux.c). Instead of initially introducing artificial viscosity, our idea involves adding a correction term to the advective part of the flux. 

$$ \frac{\partial}{\partial \bar t} \bar h_u  + \frac{\partial}{\partial \bar x} \bar Q_u = 0$$
$$ \frac{\partial}{\partial \bar t} \bar h_l  + \frac{\partial}{\partial \bar x} \bar Q_l = 0$$

solved as $\bar Q_u=\bar Q_{u1}+\bar Q_{u2}$ and 
$\bar Q_l=\bar Q_{l1}+\bar Q_{l2}$
where subscript 1 refers to the "advective" part :

$$\frac{\partial}{\partial \bar t} \bar h_u  +  \bar c_u \frac{\partial}{\partial \bar x} \bar h_u = 0 \;\;\text{ with }\;\; \bar c_u = \partial \bar Q_{u1}/\partial \bar h_u$$

$$\frac{\partial}{\partial \bar t} \bar h_l  + \bar c_l \frac{\partial}{\partial \bar x} \bar h_l = 0 \;\;\text{ with }\;\; \bar c_l = \partial \bar Q_{l1}/\partial \bar h_l$$
and subscript 2 refers to the "diffusive" part :
$$\frac{\partial}{\partial \bar t} \bar h_l  = -   \frac{\partial}{\partial \bar x} \bar Q_{l2} $$
$$\frac{\partial}{\partial \bar t} \bar h_u  = -   \frac{\partial}{\partial \bar x} \bar Q_{u2}$$
## Code
Mandatory declarations
 
*/
#include "grid/cartesian1D.h"
#include "run.h"

scalar hu[];
scalar hl[];
scalar hus[];
scalar hls[];
scalar h[];
scalar zb[];

face vector  hux[];
face vector  hlx[];
face vector Qu[];
face vector Ql[];
face vector Qud[];
face vector Qld[];
face vector cmu[];
face vector cml[];
face vector Quc[];
face vector Qlc[];

double alpha;
double X0,L0;
double M,R;
double hl0,hu0,Ql0,Qu0;
double DT;
double DTprt=20,tmax = 5000;

/**
We create a file that contains the solutions obtained for different mesh sizes, enabling us to analyze the convergence of the solution.

*/

FILE * convergence;
/**
Main with definition of parameters :

*/


int main() 
 {
    
  L0 = 200; 
  X0 = 0; 
  alpha = 1;
  M = 10;
  R = .5;
  hl0=.5;
  hu0=.5;
  Qu0=(R*(M * pow(hu0,3)/3+pow(hu0,2)*hl0)*alpha)+((hu0*pow(hl0,2)/2)* alpha);
  Ql0=(pow(hl0,3)/3)*alpha+(R*hu0*pow(hl0,2)/2)*alpha;
    
  for (int n=500; n<=500*4; n *=2)
   {
    printf("\n\n N = %d\n\n",n);
    N = n;
    DT = (L0/N)*(L0/N)/20;
    run();
 
   }
}
/**
Boundary conditions

*/

h[left] =  neumann (0);
hl[left] = dirichlet(hl0);
hu[left] = dirichlet(hu0);

hlx.n[left] = neumann (0);
hux.n[left] = neumann (0);

Ql.x[left] = dirichlet(Ql0);
Qu.x[left] = dirichlet(Qu0);
/**
Initial elevation: a "double square" of surface 2 of both fluids. \
As a simple case, we assume a constant slope, i.e. $z_b = -\alpha x$, for  the bed shape.

*/
event init (t = 0) 
 {
  foreach()
   {
    zb[] = -(x-X0)*alpha; 
    hu[] =  (fabs(x)<1)*hu0*0;
    hl[] =  (fabs(x)<1)*hl0*0 ;
      
    h[] = hu[] + hl[]; 
   }
    
  foreach_face()
   {   
    Qu.x[]=0;
    Ql.x[]=0;
        
    cmu.x[]=1;
    cml.x[]=1;
  }
}
/**
print data in stdout

*/
event printdata (t += DTprt; t <= tmax) 
 {
  if(t>100) DTprt=100;
  
  foreach()
    fprintf (stdout, "%g %g %g %g %g %g %g %g \n", x, h[], hu[], hl[], Qu.x[], Ql.x[], zb[], t);
  
  fprintf (stdout, "\n");
}
/**
Integrating the domaine over time, updating the elevation at each time step.

*/
event integration (i++) 
 {
    double dt = DT;
   
    
/**
finding the good next time step

*/
    dt = dtnext (dt);
/**
  the flux $Q(h)$ : separating in two independant flows : $Q_u$ upper fluid and $Q_l$ lower fluid
 
update of the height with mass conservation for each layer in two steps as below:
 
1- The advective part: gives a estimation of the height in the next time step $h^*$
 
$$ \frac{h^* - h^n}{\Delta t} = -\frac{1}{\Delta x} (Q^n_{1j} - Q^n_{1j-1})$$

where $Q_1$ is adjusted, for each layer separately, to prevent the shock formation.

$$ Q_{1j} = \frac{Q_{1j} + Q_{1j-1}}{2}  - c \frac{h_j - h_{j-1}}{2}$$
with $$ c = \frac{\partial Q_{1i}}{\partial h_i}$$

*/
    foreach_face()
      {
      double hlm = ((hl[0,0] + hl[-1,0])/2);
      double hum = ((hu[0,0] + hu[-1,0])/2);

            
      Qu.x[] = R * (M * pow(hum,3)/3 + pow(hum,2)*hlm) + (hum*pow(hlm,2)/2);
      Ql.x[] = (pow(hlm,3)/3) + (R * hum * pow(hlm,2)/2);
       
      cml.x[] = pow(hlm,2)+ R* hlm * hum;
      cmu.x[] = R * (M * pow(hum,2) + 2*hum*hlm) + (pow(hlm,2)/2);
          
      Quc.x[] = (Qu.x[]+Qu.x[-1])/2.  - (cmu.x[] *(hu[]-hu[-1])/2);
      Qlc.x[] = (Ql.x[]+Ql.x[-1])/2.  - (cml.x[] *(hl[]-hl[-1])/2);
        
      }
  
    foreach()
      {

      hls[] = hl[] - dt*( Qlc.x[1,0] - Qlc.x[0,0] )/Delta;
      hus[] = hu[] - dt*( Quc.x[1,0] - Quc.x[0,0] )/Delta;
              
    }
    
  
/**

2- The diffusion part : calculate the next step elevation h^{n+1}

$$ \frac{h^{n+1} - h^*}{\Delta t} = - -\frac{1}{\Delta x} (Q^*_{2j} - Q^*_{2j-1})$$

`i = u, l`: the layer index, `j`: the spatial index and `n`: the temporal index.

*/


 
    
    foreach_face()
     {
         
      hlx.x[] = ((hls[0,0] - hls[-1,0] )/Delta);
      hux.x[] = ((hus[0,0] - hus[-1,0] )/Delta);
    
      double hlsm = ((hls[0,0] + hls[-1,0])/2);
      double husm = ((hus[0,0] + hus[-1,0])/2);
          
      Qld.x[] = -((pow(hlsm,3)/3) * (R*hux.x[] + hlx.x[])) - ((R * husm * pow(hlsm,2)/2) * (hux.x[] + hlx.x[]));
      Qud.x[] = -(R * (M * pow(husm,3)/3 + pow(husm,2)*hlsm) * (hux.x[] + hlx.x[])) - ((husm*pow(hlsm,2)/2) * (R*hux.x[] + hlx.x[]));
    }
    
    
    foreach()
      {

      hl[] = hls[] - dt*( Qld.x[1,0] - Qld.x[0,0] )/Delta;
      hu[] = hus[] - dt*( Qud.x[1,0] - Qud.x[0,0] )/Delta;
    
      h[] = hl[] + hu[];
          
    }
}



event conv (t=end,last)
 {
    char name[30];
    sprintf(name,"hfinal-n%d",N);
    convergence= fopen(name,"w");
    foreach()
      fprintf(convergence,"%g %g %g %g \n", x, hl[], hu[], zb[]);
  
    fclose(convergence);
}

/**
## Plots 

Wave solution of form : $$x-ct$$

~~~gnuplot
reset
set xlabel "x-ct"
set size ratio .25
p[-30:1][] 'out' u ($1-(.163*$8)):3 w l
p[:5][] 'out' u ($1-(.146*$8)):4 w l
~~~





## Plots of convergence

At different mesh number N


~~~gnuplot
reset
set xlabel "x"
set size ratio .25
p[270:300][0:1 ]'hfinal-n500' u 1:2 w l,\
 'hfinal-n1000' u 1:2 w l,\
 'hfinal-n2000' u 1:2 w l,\
 'hfinal-n4000' u 1:2 w l,\
 'hfinal-n8000' u 1:2 w l
~~~



on se rapproche des solutions précises

~~~gnuplot
reset
set xlabel "x"
set size ratio .25
p[290:293][0:1 ]'hfinal-n500' u 1:2 w l,\
 'hfinal-n1000' u 1:2 w l,\
 'hfinal-n2000' u 1:2 w l,\
 'hfinal-n4000' u 1:2 w l,\
 'hfinal-n8000' u 1:2 w l,\
 'hfinal-n16000' u 1:2 w lp,\
 'hfinal-n32000' u 1:2 w l
~~~

 
# Links
 * [http://basilisk.fr/sandbox/M1EMN/TEST/shah.c]()
 
 
 
# Bibliography


*  Kasturi S. Shah, Samuel S. Pegler and Brent M. Minchew
Two-layer fluid flows on inclined surfaces  
J. Fluid Mech. (2021), vol. 917, A54, doi:10.1017/jfm.2021.273
 

*/