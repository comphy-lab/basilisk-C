/**
# Bicubic interpolation test case


This is a simple test case for bicubic interpolation, we use a 2D-sinusoidal
function as reference and interpolate it on a 4 times finer grid.

We also do a biquadratic interpolation and correct it using the interpolation
coefficient of the center value of the 3x3 stencil.
*/

#define BGHOSTS 2
#include "../alex_functions.h"

void get_interp(Point point, coord * p_interp, double xinterp, double yinterp){
  p_interp->x = (xinterp - x)/Delta;
  p_interp->y = (yinterp - y)/Delta;
}

#define Pi 3.14159265358979323846
double ref_function(double x, double y){
  // return (x-0.2)*(x-0.4)*(x-0.45)*(x+0.1);
  return cos(6.2*Pi*x+60);
}

int main(){
  for(int j = 0; j<=2; j++){
    int MAXLEVEL = 4+j;
    init_grid (1 << MAXLEVEL); 
    origin(-L0/2.,-L0/2.);
    scalar s[];
    foreach(){
      s[] = ref_function(x,y);
    }
    boundary({s});


/**
We only plot the reference function on the coarsest grid
*/
    int output = 4;
    if(MAXLEVEL == output){
      char filename [100];
      FILE * fp1;  
      snprintf(filename, 100,  "log0");
      fp1 = fopen (filename,"w");
  
      foreach(){
        if(y == Delta/2.){
          fprintf(fp1, "%g %g %g\n", x, y, s[]);
        }
      }
      fclose(fp1);
    }
    int nbpoint = (1 << MAXLEVEL)*4;

    double sumlin = 0., sumcub = 0., sumquad = 0.;

/**
We check the resulting interpolation on a line with y fixed.
*/
    for(int k=0;k<=nbpoint-1;k++){

      double xinterp = -L0/2.+1.*k/nbpoint;
      double yinterp = -1./(1 << MAXLEVEL);
/**
We use the locate() function the get the Point and the InterpStencil() function
to use the proper stencil for interpolation.
*/
      Point point  = locate(xinterp,yinterp);
      coord p_interp;
      get_interp(point, &p_interp,xinterp,yinterp);
      int Stencil[2];
      InterpStencil(p_interp, Stencil);

      double f_temp = bicubic( point , s, Stencil, p_interp);
      sumcub += fabs(f_temp -ref_function(xinterp,yinterp));
      double coeff[4];
      double f_temp2 = mybilin( point , s, Stencil, p_interp, coeff);
      sumlin += fabs(f_temp2 -ref_function(xinterp,yinterp));
      double f_temp3 = mybiquadratic( point , s, p_interp);
/**
Here we do a simple correction where needed, the coefficient for the central
value of the interpolated field is : $(1-x^2)(1-y^2)$.
*/      
      double error3  = f_temp - ref_function(xinterp,yinterp);
      double corr3 = (1-sq(p_interp.x))*(1.-sq(p_interp.y));
      if(corr3)f_temp3 = f_temp + error3/(corr3);
      sumquad += fabs(f_temp3 -ref_function(xinterp,yinterp));

/**
Once again, the results are only plot for the coarsest grid.
*/
      if(MAXLEVEL == output)fprintf(stderr, "%g %g %g %g %g\n",
       xinterp,yinterp,f_temp,f_temp2,f_temp3);
    }
/**
We output the L1-error.
*/
    fprintf(stdout, "%d %g %g %g\n", nbpoint, sumlin/nbpoint, sumcub/nbpoint,
      sumquad/nbpoint);
  }
}
/**
~~~gnuplot Final Results
set term svg size 1000,500
plot  'log' u 1:3 w l t 'bicub', '' u 1:4 w l t 'bilin', \
  '' u 1:5 w l t 'biquad',\
  'log0' u 1:3 w p pt 7 ps 1 lc -1 t 'ref function'   
~~~

~~~gnuplot L1-Error
f(x)  = a + b*x
f2(x) = a2 + b2*x
f3(x) = a3 + b3*x
fit f(x)  'out' u (log($1)):(log($2)) via a,b
fit f2(x) 'out' u (log($1)):(log($3)) via a2,b2
fit f3(x) 'out' u (log($1)):(log($4)) via a3,b3
ftitle(a,b) = sprintf("%.2f-%4.2f x", a, -b)
set logscale xy
set xrange [32:512]
plot  'out' u 1:2 t 'bilin', 'out' u 1:3 t 'bicubic', \
      '' u 1:4 t 'biquad',\
  exp(f(log(x))) t ftitle(a,b) , \
  exp(f2(log(x))) t ftitle(a2,b2),\
  exp(f3(log(x))) t ftitle(a3,b3)
~~~

We plot the results of the error and it works fine. Our simple correction on the
biquadratic interpolation gives us a third order accurate method.
*/