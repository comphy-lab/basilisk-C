/**
# collapse of a rectangular viscous column containing two newtonian fluids of different density and viscosity  in 1D


## Problem
Inspired by Shah et al. 2021, we consider the collapse of a two layer (with two different densities and viscosities of ratio $R$ and $M$) gravity driven flow over an inclined flat bed. The model is based on thin layer approximation and is an explicit resolution of depth-integrated continuity equations, this is a "diffusive wave problem".

$$\frac{\partial h_i}{\partial t}+  \frac{\partial Q_i(h)}{\partial x}=0$$
with $i= u, l$ corresponding to the upper and the lower fluid respectively.

The flux expressions for each layer are:

$$Q_u = - R \left(M \frac{h_u^3}{3} + h_u^2h_l\right)\left(\frac{\partial h_u}{\partial x} + \frac{\partial h_l}{\partial x} + \frac{\partial z_b}{\partial x}\right) - \frac{h_l^2 h_u}{2}\left(R \frac{\partial h_u}{\partial x} + \frac{\partial h_l}{\partial x} + \frac{\partial z_b}{\partial x}\right)$$

$$Q_l = -\frac{h_l^3}{3}\left(\frac{R \partial h_u}{\partial x} + \frac{\partial h_l}{\partial x} + \frac{\partial z_b}{\partial x}\right) - R\frac{h_l^2 h_u}{2}\left(\frac{\partial h_u}{\partial x} + \frac{\partial h_l}{\partial x} + \frac{\partial z_b}{\partial x}\right)$$

with  $M = \frac{\mu_l}{\mu_u}$ and $R = \frac{\rho_u}{\rho_l}$

In case of constant mass realised at time $t=0$, the problem is self similar.


## Code
Mandatory declarations
 
*/

#include "grid/cartesian1D.h"
#include "run.h"

scalar zb[];
scalar hl[];
scalar hu[];
scalar h[];

scalar hus[];
scalar hls[];

face vector hux[];
face vector hlx[];
face vector zbx[];

face vector Qu[];
face vector Ql[];
face vector Qud[];
face vector Qld[];

face vector Quc[];
face vector Qlc[];

face vector cmu[];
face vector cml[];

double alpha;
double L0;
double X0;
double M;
double R;
double DT;
double DTprt=10,tmax = 4000;  // @Mina:tmax changé 

/**
Main with definition of parameters :
*/
int main() 
 {
    
  L0 = 25;   // @Mina: domaine plus court
  X0 = -5.;  // @Mina: qui ne commence pas à 0 pour éviter de toucher le bord
  alpha = 1;
  M = .1;
  R = .5;
  N =  500;    // @Mina: less points
  DT = (L0/N)*(L0/N)/2;   // @Mina: sur 20 trop sévère 
  run();
}
/**
Boundary conditions

*/

h[left] = neumann (0);
hl[left] = neumann (0);
hlx.n[left] = neumann (0);
hu[left] = neumann (0);
hux.n[left] = neumann (0);
Ql.x[left] = dirichlet(0);
Qu.x[left] = dirichlet(0);
/**
Initial elevation: a "double square" of unit surface of both fluids. \
As a simple case, we assume a constant slope, i.e. $z_b = -\alpha (x-X_0)$, for  the bed shape.
*/
event init (t = 0) 
  {
  foreach()
   {  
    zb[] = -(x-X0)*alpha; 
    hl[] =  (fabs(x-.5)<.5)*.5; 
    hu[] =  (fabs(x-.5)<.5)*.5;  // longueur 1 de  x=0 à 1
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
event integration (i++) {
    double dt = DT;
    
/**
finding the good next time step
*/
    dt = dtnext (dt);
 /** 
It is common for shock waves to form in transport equations. In order to tackle this issue, our initial approach involved incorporating artificial viscosity. 
- Dirty trick : use artificial viscosity to prevent oscillations
this corresponds to the flux limiter 
 
  $\frac{\partial h_i}{\partial t}=-\frac{\partial Q_i}{\partial x}+\nu_{art}\frac{\partial^2 h_i}{\partial x^2}$
 
*/

   //foreach_face()
    //{
     //double nuart=Delta/500.;
     //Qu.x[] += - nuart*hux.x[];
     //Ql.x[] += - nuart*hlx.x[];
   //}

/**

However, we now propose a more sophisticated numerical correction technic to smooth out the solution.

Splitting method: to calculate the evolution of elevation, we update of the height using mass conservation, for each layer, in two steps as explained below:
 
1- The advective part: gives a estimation of the height in the next time step $h^*$
 
$$ \frac{h^* - h^n}{\Delta t} = -\frac{1}{\Delta x} (Q^n_{1j} - Q^n_{1j-1})$$

where $Q_1$ is adjusted, for each layer separately, to prevent the shock formation.

$$ Q_{1j} = \frac{Q_{1j} + Q_{1j-1}}{2}  - c \frac{h_j - h_{j-1}}{2}$$
with $$ c = \frac{\partial Q_{1i}}{\partial h_i}$$
 
 The flux $Q(h)$ is separated into advective and diffusive components for two independent flows: one in the upper layer and another in the lower layer.
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
   
/**
Finally update of the height with the corrected flux expression using mass conservation for each layer:

 $\frac{\partial h_u}{\partial t}= - \frac{\partial Q_u}{\partial x}$
and  $\frac{\partial h_l}{\partial t}= - \frac{\partial Q_l}{\partial x}$
*/   
   foreach()
     {

     hl[] = hls[] - dt*( Qld.x[1,0] - Qld.x[0,0] )/Delta;
     hu[] = hus[] - dt*( Qud.x[1,0] - Qud.x[0,0] )/Delta;
   
     h[] = hl[] + hu[];
     
         
   }
}



    
/**

## Results 

### Early times ($t<100$)

Collapse along the slope
~~~gnuplot
reset
set xlabel "x"
set size ratio .25
p[0:7][ ]'out' u 1:(($8<100)?$2+$7 :NaN) t'h+z_b' w l,'' u 1:($7) t'z_b' w l 
~~~


### Later times ($t>500°

@Mina : mettre "\\" pour saut de ligne dans gnuplot

~~~gnuplot
reset
set xlabel "x"
set size ratio .25
p[-2:15][0:]'out' u 1:($8>=500?$2:NaN) t'total h'w l,\
'' u 1:($8>=100?$4:NaN)  t'low' w l,''u 1:($8>=100?$3:NaN)  t'up'w l
~~~

##Similarity solution

Our problems admet a similarity solution defined as:

$$ h(\eta) = t^{-1/3}H(\eta)$$ with $$ \eta = \frac{x}{t^{1/3}}$$


The slope component becomes dominant at long times. we derive then the flux expressions in this limit by adapting the reduced equations in terms of similarity variable. Integrating these latter equations considering the flux conditions, i.e. $Q_u(x=0,t) = Q_l(x=0,t) = 0$, we ultimately obtain:

- Upper layer

$$ \frac{1}{3} \eta = R(\frac{M}{3} H_u^2 + H_u H_l) + \frac{1}{2} H_l^2$$

- Lower layer 

$$ \frac{1}{3} \eta = \frac{1}{3}H_l^2+ \frac{R}{2}H_l H_u$$

The lighter fluid gradually spills ahead of the heavier fluid. This suggests that over time, the fluids will separate into two distinct single-layer regions. This separation allows us to find the evolution law for the thickness of each layer as they move towards their respective fronts.

$$H_u = \sqrt{\frac{\eta}{RM}}\;\;\text{ for }\;\; \eta_l<\eta<\eta_u$$ 
$$H_l = \sqrt{\eta}\;\;\text{ for }\;\; 0<\eta<\eta_l$$ 

where $\eta_l$ and $\eta_u$ are the front of each layer.
  
  
- Self similarity of total $h$ 

~~~gnuplot
reset
set xlabel "x/t^{1/3}"
  set ylabel "h t^{1/3}"
set size ratio .25
p[-1:2][0: ]'out' u ($1/$8**(1./3)):($8>=600?($2*$8**(1./3)):NaN) t'total h'w l 
~~~

- Self similarity of $h_u$

~~~gnuplot
reset
set xlabel "x/t^{1/3}"
  set ylabel "h_u t^{1/3}"
set size ratio .25
p[-1:2][0: ]'out' u ($1/$8**(1./3)):($8>=600?($3*$8**(1./3)):NaN) t'h_u'w l,\
sqrt(x/(.05)) t'Anal._{upper}' w l
~~~

- Self similarity of $h_l$

~~~gnuplot
reset
set xlabel "x/t^{1/3}"
  set ylabel "h_l t^{1/3}"
set size ratio .25
p[-1:2][0: ]'out' u ($1/$8**(1./3)):($8>=600?($4*$8**(1./3)):NaN) t'h_l'w l,\
 sqrt(x) t'Anal._{lower}' w l
~~~

see figure 11 of Shah (with some oscillations...), the two phases are completely separated:

~~~gnuplot
reset
set xlabel "x/t^{1/3}"
set ylabel "h_u t^{1/3},h_lt^{1/3}"
set size ratio .25
p[-.5:2][0: ]'out' u ($1/$8**(1./3)):($8>=500?($4*$8**(1./3)):NaN) t'lower'w l,\
'out' u ($1/$8**(1./3)):($8>=500?($3*$8**(1./3)):NaN) t'upper'w l 
~~~


$$H_u = \sqrt{\frac{\eta}{RM}}\;\;\text{ for }\;\; \eta_l<\eta<\eta_u$$ 
$$H_l = \sqrt{\eta}\;\;\text{ for }\;\; 0<\eta<\eta_l$$ 

where $\eta_l$ and $\eta_u$ are the front of each layer.
 
  
@Mina:  set key left et "\\" pour les sauts de ligne dans gnuplot

~~~gnuplot
set key left
set xlabel "x/t^{1/3}"
set ylabel "t^{1/3}h_u, t^{1/3}h_l" 
tt=1500
p[-0.5:1.5][] 'out' u ($1/$8**(1./3)):($8>=tt?(($3+$4)*$8**(1./3)) : NaN) t'Comp._{upper}'w l,\
'' u ($1/$8**(1./3)):($8>=tt?($4*$8**(1./3)) : NaN) t'Comp._{lower}'w l,\
sqrt(x/(.05)) t'Anal._{upper}' w l, sqrt(x) t'Anal._{lower}' w l
~~~

 



## Links

* [http://basilisk.fr/sandbox/M1EMN/Exemples/viscous_collapse_noSV.c]()
* [http://basilisk.fr/sandbox/M1EMN/TEST/shah_flux.c]()



## References

*  Kasturi S. Shah, Samuel S. Pegler and Brent M. Minchew
Two-layer fluid flows on inclined surfaces  
J. Fluid Mech. (2021), vol. 917, A54, doi:10.1017/jfm.2021.273
 
*/ 


