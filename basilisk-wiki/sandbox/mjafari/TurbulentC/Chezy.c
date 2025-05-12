/**
##Problem

Our objective is to simulate flood currents over large spatial and temporal scales. The goal is to analyze the characteristics of these flows over long periods to validate the self-similarity of their solutions.

##Model

This problem is addressed using flood wave equations, also known as kinematic wave equations. These are a simplified, long-wave approximation of the Saint-Venant equations. At this scale, the effects of inertia and pressure variations are negligible compared to the dominant contributions from gravity and friction (refer to Whitham, p. 82; Fowler, p. 76; or Chanson). As a result, the waves propagate only downstream, and the mass and momentum conservation equations simplify to:

$$\frac{\partial h}{\partial t} +  \frac{\partial  Q}{\partial x}  =0$$
$$0 =  - \rho g h \frac{\partial z_b}{\partial x}   - \tau_0 $$
￼
In these equations:

- $z_b$ represents the bed profile, with $\partial_x z_b = - S$ being its slope, assumed constant in this example.

- $\tau_0$ is the bottom friction term, which reflects the nature and rheology of the flow.

For flood currents, turbulent friction is significant and is modeled using historical formulations. These models typically express friction as proportional to the square of the flow velocity ($\rho u^2$), with empirically defined friction coefficients. Notable contributions to these models include those by Chezy, Fanning, Manning, and Strickler.

###Simplified Turbulent Friction Model

In this study, we solve the flood wave equations using Chezy’s friction model, which is one of the simplest forms. The bottom friction is given by:

$$\tau_0 = \rho C u^2$$ 
￼
where $C$ is a constant dependent on the bed roughness. Using this friction model, the flow velocity $u$ can be approximated as:

$$u = \left( \frac{g S}{C}\right)^{1/2} h^{1/2}$$

Note that he real Chézy coefficient is $C_C=(\frac{g }{C})^{1/2}$ such as $u= C_C \sqrt{S h}$ by definition, 
and the flux $Q$ is expressed as:

$$Q = ( \frac{g S}{C})^{1/2} h^{3/2}$$￼

###Governing Equation

The main equation to solve is an advection equation, in a dimensionless form, that updates the flux at each time step:

$$\frac{\partial \bar h}{\partial \bar t} +  \frac{\partial \bar Q}{\partial \bar x}  =0$$
where:

$$ \bar Q = \beta \bar h^{3/2}$$

and $\beta$￼is the proportionality coefficient, defined as:

$\beta = \left( \frac{g S}{C}\right)^{1/2}$￼

### Code
*/


#include "grid/cartesian1D.h"
#include "run.h"

/** definition of the Parameters
 */
 
scalar h[];
scalar Q[];
double S, C;
double DTprt = 10., tmax = 500.;
char s[80];

h[left] = neumann(0);
h[right] = neumann(0);

/**
Function that computes the flux, with $n$ the power of $h$, and $\beta$ the propostionality constant.

*/
double n=3./2;  // change for tests, m=1 advection, m=2 Burgers
//double beta = pow((S/C), (1./ 2.)); 

double flux(double z, double S, double C, double n)
{
    return pow((S/C), (1./ 2.))*pow(fabs(z),n);
}
/**
 the velocity $\bar c = \partial \bar Q/\partial \bar h$ is then calculated in order to enhance the stability of the numerical scheme, using a corrective approach. (see [advection](http://basilisk.fr/sandbox/M1EMN/BASIC/advecte1.c))
*/
double celerity(double z, double S, double C, double n)
{
    return fmin(pow((S/C), (1./ 2.))*n*pow(fabs(z),n-1.),10000);
}
/**
 Main with definition of parameters, note that $S$ is the normalized slope and $C$ the friction coefficient.
 */
int main() {
    L0 = 80.;
    X0 = -1.;
    N = 1024;
    S = 1.;
    C = 0.5;
    DT = (L0/N)/16/4;

    run(); 
}
/**
 initial elevation: a constant rectangle that we release over the slope.
 */
event init (t = 0) {
    foreach()
      h[] = (x>=0 && x<=.5)*.5;
}

/**
Saving data; also some other parameters such as the tail and front positions; the maximum height:
*/

event printdata (t += DTprt; t <= tmax) {
  double xf=0, xe=0, hmax=0, xstop=0, he=0;

     foreach(){
      xf = h[] > 1e-5 ?  max(xf,x) :  xf ;
     
      hmax = max(h[], hmax);

      xstop = h[] == hmax ? max(xstop,x) : xstop;
    
      
      xe = h[] > 1e-5 ?  min(xe,x) :  xe ;
      if (xe>=0 && x<=xstop) {
        xe = (h[] > 1e-5 && h[] < 1e-4) ? max(xe,x) :  xe ;
      }
      if (x == xe){
        he =  h[];
      }
  
    }
 // sprintf (s, "shape-N%d-S%.2f-C%.2f-T%.2f.txt",N, S, C, tmax);
  sprintf (s, "shape.txt");
  FILE * fp = fopen(s, "a");
  foreach()
    fprintf(fp, "%g %g %g %g %g %g %g %g\n", 
            x, h[], t, xe, xf, hmax, xstop, he);
  fprintf(fp, "\n");
  fclose(fp);
}


/**
 integration step, at each time step
 */
event integration (i++) {
    double dt = DT;
    double cDelta = 1;
    /**
     finding the good next time step
     */
    dt = dtnext(dt);
    /**
     the algorithm is based on the flux.
     Approximation of the numerical flux taking into account
     $$Q_i = \frac{h_i^{3/2}+h_{i-1}^{3/2}}2 - c \frac{(h_i-h_{i-1})}2 $$
     */
    foreach()
    {
        cDelta = (celerity(h[], S, C, n)+celerity(h[-1], S, C, n))/2.;
        Q[] = ((flux(h[], S, C, n)+flux(h[-1], S, C, n))/2.)*S  
          - cDelta *(h[]-h[-1])/2;
    }

    foreach()
      h[] +=  - dt* ( Q[1] - Q[] )/Delta;
}


/**
####Self-Similarity
By studying the governing equations, we find stationnary solutions in the self-similar domain which are expressed as:

$$  \hat{h} \sim \hat{t}^{-2/3} \mathcal{H}(\eta), \quad \mathrm{where} \quad \eta = \hat{x}\hat{t}^{-2/3}, \quad \mathcal{H}(\eta) = \left( \frac{2 \eta}{3 \beta} \right)^2.$$
*/

/**
#plot
     

~~~gnuplot  

 set xlabel "x " 
 set ylabel "h " 
 
 p  'shape.txt'  w l
~~~



~~~gnuplot  

 set xlabel "xt^{-2/3}" 
 set ylabel "ht^{2/3}" 
 tt=1
 S=1.
 C=0.5
 V0=0.25
 etaf = (V0*S/C*(27./4))**(1./3)
 H(x)=  C/S*((2./3)*(x))**2
 Hmax=H(etaf)
 p[-0.1:2][0:1.] 'shape.txt' u ($1/($3**(2./3.))):($3>tt?($2*($3**(2./3.))):NaN) t'h' w l,etaf,\
  H(x) t ' ss' w l, 20*(x-etaf),Hmax
~~~

*/