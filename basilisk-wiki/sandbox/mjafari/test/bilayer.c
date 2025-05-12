/**
# A bi-layer gravity driven flow

![animation  ](bilayer/animate-B0.00.gif)

## Colapse of a two layer heap

We use the Saint-Venats equations to solve the evolution of a two layer flow on a flat or an inclined bed. 
The model is based on thin layer approximation and is an explicit resolution of depth-integrated continuity equations (diffusive wave).

$$\frac{\partial h_i}{\partial t}+  \frac{\partial Q_i(h)}{\partial x}=0$$

$$Q_0 = (\frac{Y^3}{3} + h_0^2Y - Y^2h_0)(S-\partial_xp_0) + ((h_0 - \frac{Y}{2})Yh_1)(RS - \partial_xp_1) - (h_0 - \frac{Y}{2})YB.$$
     $$Q_1 = (h_0 - \frac{Y}{2})Yh_1 (S-\partial_xp_0) + M\frac{h_1^3}{3}(RS - \partial_xp_1) - Yh_1B.$$
 

## Codes

Mandatory declarations:
 */
#include "grid/cartesian1D.h"
#include "run.h"

scalar h0[], h1[];
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
double DTprt=.1,tmax = 20;
double M, R;
char s[80];
/**
 Main with definition of parameters
 */
int main() {
    L0 = 2.*4;
    X0 = -1;
    M = 1.;
    R = 1.;
    N = 128*8;
    DT = (L0/N)*(L0/N)/5;
    B = 0.;
    S = 0;
    run();
    S = 1;
}
event init (t = 0) {
    foreach(){
        h0[] = (fabs(x)<.25);
        h1[] = 0*(fabs(x)<.25);
        Q0[]=0;
        Q1[]=0;
    }
}
     /**
         Tracking the front and the end of the heap
         */
        
event output (t += inct; t < tmax) {
  double  x0f=0,x0e=0;
  double  x1f=0,x1e=0;
        inct = inct*2;
        inct = min(1,inct);
      
       
        foreach(){
            x0f = h0[] > 1e-4 ?  max(x0f,x) :  x0f;
            x0e = h0[] > 1e-4 ?  min(x0e,x) :  x0e;
            
            x1f = h1[] > 1e-4 ?  max(x1f,x) :  x1f;
            x1e = h1[] > 1e-4 ?  min(x1e,x) :  x1e;
        }
        
        sprintf (s, "front-%.2f.txt", B);
        FILE * f = fopen (s, "w");
        foreach()
        fprintf (f, "%g %g %g %g %g %g \n", fmin((x-x0f),0), h0[], x0e-x0f, fmin((x-x1f),0), h1[], x1e-x1f);
        fclose(f);
    
        sprintf (s, "x-%.2f.txt", B);
        FILE * fp = fopen (s, "a");
        fprintf (fp, "%g %g %g \n", t, x0f, x1f);
        fclose(fp);
    }

    /**
     Save the hight the flux and the yield surface as a function of time
     */

event printdata (t += DTprt; t <= tmax){
        sprintf (s, "shape-B%.2f-S%.2f.txt", B, S);
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
        // tau[] = (R*S - pxc1[])*h1[];
        //beta[] = (R*S - px1[])/(fabs(S -  pxc0[]));
        beta[] = 0.;
        //Y[] =  max(h0[] + beta[]*h1[] - B/fabs(S -  pxc0[])  ,0);
        Y[] =  h0[];
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
            
        Q0[] = - nu0[]*(-S + px0[]) + nu01[] * h1[]*(R*S - pxc1[]) - nu01[]*B;
        Q0star[] = Q0[] - c0[]*(h0[] - h0[-1])/2.;
            
        Q1[] = - nu10[]*(-S + px0[]) - M*nu1[]*(R*S - px1[])  + nu1B[]*h1[]*(R*S - pxc1[]) - nu1B[]*B;
        Q1star[] = Q1[] - c1[]*(h1[] - h1[-1])/2.;
        }
        
    /**
    derivative $O(\Delta)$ up stream like for heat equation
    */
    foreach(){
        dQ0[] =  ( Q0star[1,0] - Q0star[0,0] )/Delta;
        dQ1[] =  ( Q1star[1,0] - Q1star[0,0] )/Delta;
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
animation
*/
event animatedplot (t+=.1) {
    static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
    if(t==0) fprintf (fp,"set term gif animate;set output 'animate-B%.2f.gif';set size ratio .333333\n", B);
    fprintf (fp,"\nset grid\n");
    fprintf (fp,"set title 'Collapse en 1D --- t= %.2lf '\n"
      "t= %.2lf  ; "
      "p[-1:1][-1:1]  '-' u 1:($2+$3+$4) t'free surface' w l lt 3,"
      "'' u 1:($2+$4) t'lower layer' w l lt 4,\\\n"
      "'' u 1:4 t'topo/L0' w l lt -1\\\n",
           t,t);
    foreach()
    fprintf (fp,"%g %g %g %g %g\n", x, h0[], h1[], zb[]/L0, t);
    fprintf (fp,"e\n\n");
    fflush (fp);
    if(t==tmax){
        fprintf (fp,"! cp animate.gif a2.gif \n");
    }
}

/**
#plot
     
~~~gnuplot
 p[-1:1][] 'shape-B0.00-S0.00.txt' u ($1/$8**0.2):($2*$8**0.2) w l
~~~
 
*/

