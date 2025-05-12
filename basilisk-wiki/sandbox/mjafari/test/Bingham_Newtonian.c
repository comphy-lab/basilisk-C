/**
# A bi-layer gravity driven flow


## Colapse of a two layer heap

We use the Saint-Venats equations to solve the evolution of a two layer flow on a flat or an inclined bed. 
The model is based on thin layer approximation and is an explicit resolution of depth-integrated continuity equations (diffusive wave).

$$\frac{\partial h_i}{\partial t}+  \frac{\partial Q_i(h)}{\partial x}=0$$

$$Q_0 = (\frac{Y^3}{3} + h_0^2Y - Y^2h_0)(S-\partial_xp_0) + (h_0 - \frac{Y}{2})Y(\tau_{10}-B).$$
     $$Q_1 = (h_0 - \frac{Y}{2})Yh_1 (S-\partial_xp_0) + m\frac{h_1^3}{3}(RS - \partial_xp_1) + Yh_1(\tau_{10}-B).$$
with $\tau_{01}$ the shear at the interface of the two layers and $Y$ represents the boundary separating the shearing flow from the plug flow.
$$\tau_{01} = h_1(rS -\partial_{x}p_1)$$
$$Y = h_0 +\frac{(\tau_{10}-B)}{\lvert(S -\partial_xp_0)\lvert}$$

## Code

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
scalar nu0[], nu0B[];
scalar nu1[], nu10[], nu1B[];
scalar Q0[], Q0star[], Q1[], Q1star[];
scalar dQ0[], dQ1[];
scalar c0[], c1[], zb[], beta[], Bt[];
double dt,S,B,inct=0.0001;
double DTprt=10,tmax = 600;
double M, R,eps,l;
int i,n;
char s[80];
/**
 Main with definition of parameters
 */
int main() {
    L0 = 2.*8;
    X0 = -1.;//-L0/2.;
    M = 20.;
    R = .2;
    N = 128*4;
    DT = (L0/N)*(L0/N)/5;
    eps = pow(10,-10);
    S = 1.;
    n = 0.;
    for (i=0; i<=1; i+=1){
        for (B=0.; B<=.5; B+=.5){  
            l=i;
            sprintf (s, "shape-l%.2f-B%.2f.txt",l, B);
            FILE * fp = fopen (s, "a");
            fclose(fp);
            run();
        
        }
    }

}


event init (t = 0) {
    foreach(){
        h0[] = (fabs(x)<.25);
        h1[] = l*(fabs(x)<.25);
        zb[] = -(x-X0)*S;
        c0[]=0;
        c1[]=0;
        Q0[]=0;
        Q1[]=0;
    }
}

    /**
     Save the hight the flux and the yield surface as a function of time
     */

event printdata (t += DTprt; t <= tmax){
    sprintf (s, "shape-l%.2f-B%.2f.txt",l, B);
    FILE * fp = fopen (s, "a");
    foreach()
    fprintf (fp, "%g %g %g %g %g %g %g %g %g \n", x, h0[], h1[], Q0[], Q1[], Y[], tau[], t, zb[]);
    fprintf (fp, "\n");
    fclose(fp);
    }

     /**
     Tracking the front of the flow
     */

event outputfront (t += DTprt, t<= tmax ) {
    double  x0f=0,x0e=0, x1f=0,x1e=0;
    foreach(){
        x0f = h0[] > 1e-5 ?  max(x0f,x) :  x0f ;
        x0e = h0[] > 1e-5 ?  min(x0e,x) :  x0e ;
        x1f = h1[] > 1e-5 ?  max(x1f,x) :  x1f ;
        x1e = h1[] > 1e-5 ?  min(x1e,x) :  x1e ;
    }

    sprintf (s, "front-l%.2f-B%.2f.txt",l, B);
    FILE * f = fopen (s, "w");
    foreach()
        fprintf (f, "%g %g %g %g %g %g \n", fmin((x-x0f),0), h0[], x0e-x0f, fmin((x-x1f),0), h1[], x1e-x1f);
    fclose(f);
}


    /**
     The fluxes are written as: $$Q_0= FQ_{0} (S-\partial_xp_0) + FQ_{0B} (\tau{10}-B)$$
     $$Q_1= FQ_{10} (S-\partial_xp_0) + MFQ_{1}(RS - \partial_xp_1) + FQ_{1B}(\tau_{10}-B)$$

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
    Here we add a condition that allows us to transition from Bingham-Newtonian model to the Newtonian-Newtonian bi-layer flow.
    Yield surface can be writtern as:
    $$Y = h_0 +\frac{(\tau_{10}-B)}{\lvert(S -\partial_xp_0)\lvert}\hspace{0.7cm}\textrm{if}\hspace{0.2cm} \tau_{10}<B $$
    $$Y = h_0.\hspace{0.7cm}\textrm{if}\hspace{0.2cm} \tau_{10}>B$$
    
    $\beta$ represent derivative of the yield surface $Y$ in respect to  $h_1$ used in the calculation of the advection celerity $c_i$
    
    $$Y = h_0 + \underbrace{\frac{(RS - \partial_xp_1)}{(S -\partial_xp_0)}}_\textrm{$\beta$} h _1 - \frac{B}{(S -\partial_xp_0)}$$
 
    */
        
    foreach(){
        tau[] = (R*S - px1[])*h1[];
        beta[] = (R*S - px1[])/(fabs(S -  pxc0[])+eps);
      //  Y[] =  max(h0[] + (tau[] - B)/(fabs(S -  pxc0[])+eps)  ,0);
      // PYL 
        Y[] =  max(min(h0[],h0[] + (tau[] - B)/(fabs(S -  pxc0[])+eps)),0);
        if(tau[]>B){  
            // Y[] = h0[];
            // PYL
             n=1.;  //n ??
             beta[]=0.;
        }
       
    }
      
    /**
    The flux is taken with the mean value with the next cell, $\nu$ is a kind of viscosity (that is why we center)
    */
        
    foreach(){
          
        nu0[] = (FQ0(h0[-1,0],Y[-1,0]) + FQ0(h0[0,0],Y[0,0]))/2;
        nu0B[] = (FQ01(h0[-1,0],Y[-1,0]) + FQ01(h0[0,0],Y[0,0]))/2;
            
        nu1[] = (FQ1(h1[-1,0]) + FQ1(h1[0,0]))/2;
        nu10[] = (FQ10(h0[-1,0],h1[-1,0],Y[-1,0]) + FQ10(h0[0,0],h1[0,0],Y[0,0]))/2;
        nu1B[] = (FQB1(h1[-1,0],Y[-1,0]) + FQB1(h1[0,0],Y[0,0]))/2;

    }
        
    /**
    The flux is decomposed into its diffusive and advective components. With the help of the advection celerity, which is  expressed as $c_i = \frac{\partial Q_i}{\partial h_i}$, we introduce corrections to the expressions of each fluxe to enhance the stability of the scheme, as:
    $$Q_i^* = Q_i - c_i \frac{(h_i^j - h_i^{j-1})}{2}$$
    with
    $$ c_0 = h_0YS + h_0h_1RS.$$
    $$c_1 =(\beta h_1Y + \frac{Y^2}{2})S + (Mh_1^2 + 2h_1h_0)RS.$$
    */
        
    foreach(){
        c0[] = (h0[]*Y[])*S + n*h0[]*h1[]*R*S;
        c1[] = (Y[]*Y[]/2. + beta[]*h1[]*Y[])*S + (M*pow(h1[],2)+n*2*h1[]*h0[])*R*S;

        Q0[] = - nu0[]*(-S + px0[]) + nu0B[]*(tau[]-B);
        Q1[] = - nu10[]*(-S + px0[]) - M*nu1[]*(-R*S + px1[]) + nu1B[]*(tau[]-B);

        Q0star[] = Q0[] - c0[]*(h0[] - h0[-1])/2.;
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
#plot
     

~~~gnuplot Huppert's second problem

 set xlabel "xt^{-1/3}" 
 set ylabel "ht^{1/3}" 
 tt=300
 p[-1:1][0:] 'shape-l0.00-B0.00.txt' u ($1/($8**(1./3.))):($8>tt?($2*($8**(1./3.))):NaN) t'h' w l,\
  sqrt(x+.1) w l
~~~

~~~gnuplot collapse of a one layer Bingham heap
 set xlabel "x"
 set ylabel "h + (B/S)Log((B-(S h)/B)" 
 B=0.5
 tt=50
 p[-1:2][:1]'shape-l0.00-B0.50.txt' u 1:($8>tt?$2:NaN) t'B=0.50' w l,\
  ''u (($2) + B*log((B-$2)/B) )+.88:($8>tt?$2:NaN) t'exact' w l

~~~
 
~~~gnuplot collapse of a Newtonian-Niewtonian bi-layer

 set xlabel "xt^{-1/3}" 
 set ylabel "ht^{1/3}" 
 R=.2
 M=20
 tt=300
 p[-1:2][] 'shape-l1.00-B0.00.txt' u ($1/($8**(1./3.))):($8>tt?($2*($8**(1./3.))):NaN) t'lower' w l,\
 '' u ($1/($8**(1./3.))):($8>tt?($3*($8**(1./3.))):NaN) t'upper' w l,sqrt(x+.1) w l t'auto-sim',\
  sqrt((x+.1)/(R*M)) t'auto-sim' w l 

~~~

~~~gnuplot collapse of a Bingham-Niewtonian bi-layer

 set xlabel "xt^{-1/3}" 
 R=.2
 M=20
 tt=300
  p[-1:2][] 'shape-l1.00-B0.50.txt' u ($1/($8**(1./3.))):($8>tt?($3*($8**(1./3.))):NaN) t'upper' w l,\
  sqrt((x+.15)/(R*M)) t'auto-sim' w l


~~~

~~~gnuplot

set xlabel "x"
set ylabel "h + (B/S)Log((B-(S h)/B)" 
B=0.5
p[][]'front-l1.00-B0.50.txt' t 'B=0.5' w l,''u (($2) + B*log((B-$2)/B) ):2 t'exact' w l

~~~

*/

