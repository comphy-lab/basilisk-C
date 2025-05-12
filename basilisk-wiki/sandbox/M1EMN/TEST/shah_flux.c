/**
# Avalanche of two newtonian fluids of different density and viscosity in 1D

This is Shah 2021 first problem of a two layer gravity driven flow (with two different densities and viscosities of ratio $R$ and $M$) over an inclined flat plate. The model is based on thin layer approximation and is an explicit resolution of depth-integrated continuity equations (diffusive wave).

$$\frac{\partial h_i}{\partial t}+  \frac{\partial Q_i(h)}{\partial x}=0$$
with $i= u, l$ corresponding to the upper and the lower fluid respectively.

The flux expressions for each layer are:

$$Q_u = - R \left(M \frac{h_u^3}{3} + h_u^2h_l\right)\left(\frac{\partial h_u}{\partial x} + \frac{\partial h_l}{\partial x} + \frac{\partial z_b}{\partial x}\right) - \frac{h_l^2 h_u}{2}\left(R \frac{\partial h_u}{\partial x} + \frac{\partial h_l}{\partial x} + \frac{\partial z_b}{\partial x}\right)$$

$$Q_l = -\frac{h_l^3}{3}\left(\frac{R \partial h_u}{\partial x} + \frac{\partial h_l}{\partial x} + \frac{\partial z_b}{\partial x}\right) - R\frac{h_l^2 h_u}{2}\left(\frac{\partial h_u}{\partial x} + \frac{\partial h_l}{\partial x} + \frac{\partial z_b}{\partial x}\right)$$

with  $M = \frac{\mu_l}{\mu_u}$ and $R = \frac{\rho_u}{\rho_l}$

We aim to  track the evolutions of the flow fronts.


For initial layer $h_{u0}, h_{l0}$ has associated fluxes 

$Q_{l0}=\frac{h_{l0}^3}{3} \alpha+Rh_{u0}\frac{h_{l0}^2}{2}\alpha$ 

$Q_{u0}=(R(M  h_{u0}^3/3+h_{u0}^2*hl0)\alpha)+((h_{u0}(h_{l0}^2)/2)\alpha)$



## Code
Mandatory declarations
 
*/

#include "grid/cartesian1D.h"
#include "run.h"

scalar hu[];
scalar hl[];
scalar h[];
scalar zb[];
face vector  hux[];
face vector  hlx[];
face vector  zbx[];
face vector Qu[];
face vector Ql[];

double alpha;
double X0,L0;
double M,R;
double hl0,hu0,Ql0,Qu0;
double DT;
double DTprt=20,tmax = 500;

/**
Main with definition of parameters :

*/
int main() {
    
    L0 = 40;
    X0 = -5;
    alpha = 1;
    M = .1;
    R = .5;
    hl0=.25;
    hu0=.75;
    Ql0=(pow(hl0,3)/3)*alpha+(R*hu0*pow(hl0,2)/2)*alpha;
    Qu0=(R*(M * pow(hu0,3)/3+pow(hu0,2)*hl0)*alpha)+((hu0*pow(hl0,2)/2)* alpha);
  
    N =  500;
    DT = (L0/N)*(L0/N)/10;
    run();
}
/**
Boundary conditions :  constant injection at the left boundary
*/

h[left] =  neumann (0);
hl[left] = dirichlet(hl0); 
hlx.n[left] = neumann (0);
hu[left] = dirichlet(hu0); 
hux.n[left] = neumann (0);
Ql.x[left] = dirichlet(Ql0);
Qu.x[left] = dirichlet(Qu0);
/**
No initial elevation. \
As a simple case, we assume a constant slope, i.e. $z_b = -\alpha x$, for  the bed shape.
*/
event init (t = 0) {
  foreach(){
      
    zb[] = -(x-X0)*alpha; //bottom topography
    hl[] =   0 ;
    hu[] =   0;
      
    h[] = hu[] + hl[];   //total height 
    }
    foreach_face(){
      Qu.x[]=0;
      Ql.x[]=0;
  }
  boundary ({hu,Qu});   
  boundary ({hl,Ql});
  }
/**
print data in stdout
*/
event printdata (t += DTprt; t <= tmax) {   //500
  if(t>100) DTprt=100;
  foreach()
    fprintf (stdout, "%g %g %g %g %g %g %g %g \n", x, h[], hu[], hl[], Qu.x[], Ql.x[], zb[], t);
  fprintf (stdout, "\n");
}
/**
Integrating the domaine over time, updating the elevation at each time step.
*/
event integration (i++) {
    double dt = DT;
    
/**
finding the good next time step
*/
    dt = dtnext (dt);
/**
  the flux $Q(h)$ : separating in two independant flows
*/
    foreach_face(){
       
     hlx.x[] = ((hl[0,0] - hl[-1,0])/Delta);
     hux.x[] = ((hu[0,0] - hu[-1,0])/Delta);
        
     zbx.x[] = ((zb[0,0] - zb[-1,0])/Delta);
    
     double hlm = ((hl[0,0] + hl[-1,0])/2);
     double hum = ((hu[0,0] + hu[-1,0])/2);

            
     Ql.x[] = -((pow(hlm,3)/3) * (alpha*(R*hux.x[] + hlx.x[]) + zbx.x[]))
        - ((R * hum * pow(hlm,2)/2) * ( alpha*(hux.x[] + hlx.x[]) + zbx.x[]));

     Qu.x[] = -(R*(M*pow(hum,3)/3+pow(hum,2)*hlm)*(alpha*(hux.x[]+ hlx.x[]) + zbx.x[])) 
       - ((hum*pow(hlm,2)/2) * (alpha*(R*hux.x[] + hlx.x[]) + zbx.x[]));
  }    
    boundary ({Qu});
    boundary ({Ql});   
    
/** 
 dirty trick : use artificial viscosity to prevent oscillations
 this corresponds to the flux limiter 
 
  $\frac{\partial h_i}{\partial t}=-\frac{\partial Q_i}{\partial x}+\nu_{art}\frac{\partial^2 h_i}{\partial x^2}$
 
*/
   foreach_face(){
     double nuart=Delta/100;
     Qu.x[] += - nuart*hux.x[];
     Ql.x[] += - nuart*hlx.x[];
   }
    boundary ({Qu});
    boundary ({Ql});
    
/**
update of the height with mass conservation for each layer:

 $\frac{\partial h_u}{\partial t}= - \frac{\partial Q_u}{\partial x}$
and  $\frac{\partial h_l}{\partial t}= - \frac{\partial Q_l}{\partial x}$
*/
    foreach(){     
      hl[] -=  dt*( Ql.x[1,0] - Ql.x[0,0] )/Delta;
      hu[] -=  dt*( Qu.x[1,0] - Qu.x[0,0] )/Delta;  

      h[] = hl[] + hu[];
  }
        boundary ({hl});
        boundary ({hu});
    
}

/**
## Results

note the small oscillations (artificial visocity has to be tuned)


early times

~~~gnuplot
reset
set xlabel "x"
set size ratio .25
p[-1.5:7][ ]'out' u 1:(($8<200)?$2+$7 :NaN) t'h+z_b' w l,'' u 1:($7) t'z_b' w l 
~~~


later times

~~~gnuplot
reset
set xlabel "x"
set size ratio .25
tt=10
p[][0: ]'out' u 1:($8>=tt?$2:NaN) t'total h'w l,'' u 1:($8>=tt?$4:NaN)  t'low' w l,''u 1:($8>=tt?$3:NaN)  t'up'w l
~~~

tracking the front, it is moving at constant velocity:

~~~gnuplot
reset
set xlabel "x"
set size ratio .25
tt=10
p[-10:0][0: ]'out' u (5.5+$1-$8/500*40):($8>=tt?$2:NaN) t'total h'w l 
~~~


## Links

* [http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c]()
* [http://basilisk.fr/sandbox/M1EMN/TEST/shah.c]()




## References

*  Kasturi S. Shah, Samuel S. Pegler and Brent M. Minchew
Two-layer fluid flows on inclined surfaces  
J. Fluid Mech. (2021), vol. 917, A54, doi:10.1017/jfm.2021.273


*/ 
