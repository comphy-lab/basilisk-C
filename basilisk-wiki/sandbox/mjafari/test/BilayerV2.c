/**
# A bi-layer gravity driven flow

![animation  ](bilayer/animate-B0.00.gif)

## Colapse of a two layer heap

We use the Saint-Venats equations to solve the evolution of a two layer flow on a flat or an inclined bed. 
The model is based on thin layer approximation and is an explicit resolution of depth-integrated continuity equations (diffusive wave).

$$\frac{\partial h_i}{\partial t}+  \frac{\partial Q_i(h)}{\partial x}=0$$

$$Q_0 = (\frac{Y^3}{3} + h_0^2Y - Y^2h_0)(S-\partial_xp_0) + ((h_0 - \frac{Y}{2})Yh_1)(RS - \partial_xp_1) - (h_0 - \frac{Y}{2})YB.$$
     $$Q_1 = (h_0 - \frac{Y}{2})Yh_1 (S-\partial_xp_0) + (M\frac{h_1^3}{3} +Yh_1^2)(RS - \partial_xp_1) - Yh_1B.$$
 
## Transition between models

**Mono-layer models**

_Huppert_: The yield stress is zero. We investigate the collapse of a viscous heap over an inclined and a horizontal bed (i.e., $S=1$ and $S=0$ respectively).
$$\boxed{B=0 \Rightarrow Y=h_0}$$

_Balmforth_ : Collapse of a single-layered Bingham heap, with yield stress $B=0.5$.
$$\boxed{B\neq 0 \Rightarrow Y = h_0 - \frac{B}{(S -\partial_xp_0)}}$$

**Bi-layer models**

_Shah_ : Collapse of a bi-layer viscous heap
$$\boxed{B=0 \Rightarrow Y = h_0 +\frac{(RS - \partial_xp_1)}{(S -\partial_xp_0)}h_1}$$

_Roya_ : Collapse of a viscous-bingham heap
$$\boxed{B\neq 0 \Rightarrow Y = h_0 +\frac{(RS - \partial_xp_1)}{(S -\partial_xp_0)}h_1 - \frac{B}{(S -\partial_xp_0)}}$$

## Codes

Mandatory declarations:
 */
#include "grid/cartesian1D.h"
#include "run.h"

scalar h0[], h1[], h1s[];
scalar hx0[], hx1[];
scalar hxc0[], hxc1[];
scalar p0[], p1[];
scalar px0[], px1[];
scalar pxc0[], pxc1[];
scalar Y[], tau[];
scalar nu0[], nu01[];
scalar nu1[], nu10[], nu1B[];
scalar Q0[], Q0star[], Q1[], Q1star[];
scalar dQ0[], dQ1[];
scalar c0[], c1[], zb[], beta[];
double dt,S,B,inct=0.0001;
double DTprt=30,tmax = 600;
double M, R,eps,l,i;
char s[80];
/**
 Main with definition of parameters
 */
int main() {
    L0 = 2.*8;
    X0 = -1;
    M = 5.;
    R = .2;
    N = 128*4;
    DT = (L0/N)*(L0/N)/5;
    eps = pow(10,-10);
    for (i=0.; i<=1.; i+=1){
        for (B=0.; B<=.5; B+=.5){
            for (S=0.; S<=1.; S+=1.){
                l=i;
                sprintf (s, "shape-l%.2f-B%.2f-S%.2f.txt",l, B, S);
                FILE * fp = fopen (s, "a");
                fclose(fp);
                run();
            }
        }
    }
}


event init (t = 0) {
    foreach(){
        h0[] = (fabs(x)<.25);
        h1[] = l*(fabs(x)<.25);
        zb[] = -(x-X0)*S;
        Q0[]=0;
        Q1[]=0;
    }
}

    /**
     Save the hight the flux and the yield surface as a function of time
     */

event printdata (t += DTprt; t <= tmax){
        sprintf (s, "shape-l%.2f-B%.2f-S%.2f.txt",l, B, S);
        FILE * fp = fopen (s, "a");
        foreach()
        fprintf (fp, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g \n", x, h0[], h1[], Q0[], Q1[], Y[], tau[], t, zb[], beta[], px0[], px1[], pxc0[], pxc1[]);
        fprintf (fp, "\n");
        fclose(fp);
    }

    /**
     The fluxes are written as: $$Q_0= FQ_0 (S-\partial_xp_0) + FQ_{01}\tau_{01} - FQ_{01}B$$
     $$Q_1= FQ_10 (S-\partial_xp_0) + FQ_{1}M(RS - \partial_xp_1) + FQ_{1B}\tau_{01} - FQ_{1B}B$$
     where $\tau_{01}$ is the shear at the interface between the two flows.
     $$\tau_{01} = (RS - \partial_xp_1)h_1.$$
     */
    
double FQ0(double h, double Y)
{
    return 1./3.*pow(Y,3) + pow(h,2)*Y - pow(Y,2)*h;
        
}

double FQ01(double h,double Y)
{
    return (h - 1./2.*Y)*Y;
}
    
double FQ1(double h)
{
    return 1./3.*pow(h,3);
}
    
double FQ10(double hl, double hu, double Y)
{
    return (hl - 1./2.*Y)*Y*hu;
}
    
double FQB1(double hu, double Y)
{
    return Y*hu;
}
    
event integration (i++) {
    double dt = DT;
        
    /**
    Finding the good next time step
    */
        
    dt = dtnext (dt);
        
    /**
    $O(\Delta)$ down stream derivative
    */
        
    foreach(){
        hx0[] =  ( h0[0,0] - h0[-1,0] )/Delta;
        hx1[] =  ( h1[0,0] - h1[-1,0] )/Delta;
    }
        
    /**
    And the derivative of the pressure in each layer
    */
        
    foreach(){
        px0[] = R*hx1[] + hx0[];
        px1[] = R*(hx1[] + hx0[]);
    }
        
   /**
    Centered derivative for $h_i$ used in the yield criteria for the yield surface $Y$
    */
        
    foreach(){
        hxc0[] =  ( h0[1,0] - h0[-1,0] )/2/Delta;
        hxc1[] =  ( h1[1,0] - h1[-1,0] )/2/Delta;
    }
        
    foreach(){
        pxc0[] = R*hxc1[] + hxc0[];
        pxc1[] = R*(hxc1[] + hxc0[]);
    }
        
    /**
    Yield surface can be writtern as:
    $$Y = h_0 + \underbrace{\frac{(RS - \partial_xp_1)}{(S -\partial_xp_0)}}_\textrm{$\beta$} h _1 - \frac{B}{(S -\partial_xp_0)}$$
      
    $\beta$ represent derivative of the yield surface $Y$ in respect to  $h_1$ used in the calculation of the advection celerity $c_i$
    */
        
    foreach(){
        tau[] = (R*S - px1[])*h1[];
        beta[] = (R*S - px1[])/(fabs(S -  pxc0[])+eps);
        Y[] =  max(h0[] + (tau[] - B)/(fabs(S -  pxc0[])+eps)  ,0);
       
    }
      
    /**
    The flux is taken with the mean value with the next cell, $\nu$ is a kind of viscosity (that is why we center)
    */
        
    foreach(){
          
        nu0[] = (FQ0(h0[-1,0],Y[-1,0]) + FQ0(h0[0,0],Y[0,0]))/2;
        nu01[] = (FQ01(h0[-1,0],Y[-1,0]) + FQ01(h0[0,0],Y[0,0]))/2;
            
        nu1[] = (FQ1(h1[-1,0]) + FQ1(h1[0,0]))/2;
        nu10[] = (FQ10(h0[-1,0],h1[-1,0],Y[-1,0]) + FQ10(h0[0,0],h1[0,0],Y[0,0]))/2;
        nu1B[] = (FQB1(h1[-1,0],Y[-1,0]) + FQB1(h1[0,0],Y[0,0]))/2;
    }
        
    /**
    The flux is decomposed into its diffusive and advective components. With the help of the advection celerity, which is  expressed as $c_i = \frac{\partial Q_i}{\partial h_i}$, we introduce corrections to the expressions of each fluxe to enhance the stability of the scheme, as:
    $$Q_i^* = Q_i - c_i \frac{(h_i^j - h_i^{j-1})}{2}$$
    with
    $$ c_0 = h_0^2S + h_0h_1RS - h_0B.$$
    $$c_1 =( h_0 + 2\beta h_1 - \frac{Y}{2})YS + (h_1^2(M+\beta) + 2h_1Y)RS - (Y + \beta h_1)B.$$
    */
        
    foreach(){
            
        c0[] = pow(h0[],2)*S + h0[]*h1[]*R*S - h0[]*B;
        c1[] = (h0[] + 2*beta[]*h1[] - Y[]/2.)*Y[]*S + (pow(h1[],2)*(M+beta[])+ 2*h1[]*Y[])*R*S - (Y[] + beta[]*h1[]) * B;
            
        Q0[] = - nu0[]*(-S + px0[]) - nu01[] * h1[]*(-R*S + px1[]) - nu01[]*B;
        Q0star[] = Q0[] - c0[]*(h0[] - h0[-1])/2.;
            
        Q1[] = - nu10[]*(-S + px0[]) - M*nu1[]*(-R*S + px1[])  - nu1B[]*h1[]*(-R*S + px1[]) - nu1B[]*B;
        Q1star[] = Q1[] - c1[]*(h1[] - h1[-1])/2.;
        }
        
    /**
    derivative $O(\Delta)$ up stream like for heat equation
    */
    foreach(){
        dQ0[] =  ( Q0[1,0] - Q0[0,0] )/Delta;
        dQ1[] =  ( Q1[1,0] - Q1[0,0] )/Delta;
    }
    /**
    update $h$
    */
    foreach(){
        h0[] +=  - dt*dQ0[];
        h1[] +=  - dt*dQ1[];
    }
        
}

/**
#plot
     
~~~gnuplot Huppert's first problem

 set xlabel "x"
 set ylabel "h(x,t)" 
 set title 'B=0, S=0'
 p[-1:][] 'shape-l0.00-B0.00-S0.00.txt' u 1:2 w l 
 ~~~

~~~gnuplot Auto-similarity
 set xlabel "xt^{-1/5}"
 set ylabel "ht^{1/5}" 
 b= 0.493108
 p[0:1.5] 'shape-l0.00-B0.00-S0.00.txt' u ($1/($8**.2)):($2*($8**.2)) t 'num' w l,(9./10*(b*b-x*x))**(1/3.) t'Self Sim'
 ~~~

~~~gnuplot Huppert's second problem
 set xlabel "x"
 set ylabel "h(x,t)" 
 set title 'B=0, S=1'
 p[-1:][] 'shape-l0.00-B0.00-S1.00.txt' u 1:2 w l 
 ~~~

 ~~~gnuplot analytical solution
 set xlabel "xt^{-1/3}" 
 set ylabel "ht^{1/3}" 
 p[-1:8][0:] 'shape-l0.00-B0.00-S1.00.txt' u ($1/($8**(1/3))):($2*($8**(1/3))) t'comp.' w l, sqrt(x/(980)+0.0005) t'anal'' w l
~~~

~~~gnuplot collapse of a one layer Bingham heap
 set xlabel "x"
 set ylabel "h(x,t)" 
 set title 'B=0.5, S=0'
 p[-1:1][] 'shape-l0.00-B0.50-S0.00.txt' u 1:2 w l 
~~~
~~~gnuplot
 set xlabel "x"
 set ylabel "h(x,t)" 
 set title 'B=0.5, S=1'
 p[-1:1.5][] 'shape-l0.00-B0.50-S1.00.txt' u 1:2 w l 
~~~

~~~gnuplot exact solution
 set xlabel "h + B /S Log(B-(S h)/B)"
 set ylabel "h(x,t)" 
 B=0.5
 p[-1:2][:1]'shape-l0.00-B0.50-S1.00.txt' u 1:($8>500?$2:NaN) t'B=0.50' w l, ''u (($2) + B*log((B-$2)/B) )+.88:($8>500?$2:NaN) t'exact' w l

~~~
 
*/
