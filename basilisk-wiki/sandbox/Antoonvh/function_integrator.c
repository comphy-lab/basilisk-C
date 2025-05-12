/**
# Numerical integration

We aim to integrate a function numerically, and employ an adaptive method.

This is our concern:
$$\int_{-100}^{100}e^{-x^2}\mathrm{d}x \approx \int_{-\infty}^{\infty}e^{-x^2}\mathrm{d}x = \sqrt{\pi}$$
*/

#include "grid/bitree.h"
#include "utils.h"

#define fun(x) (exp(-sq(x)))
#define DIGITS 6

scalar s[];
  
int main(){
  L0 = 200;
  X0 = -L0/2;
  init_grid(128);
  double se = 0.01;
  double an = HUGE, ao = HUGE;
  long unsigned int n = 0;
  do{
    ao = an;
    foreach(){
       s[] = fun(x);
       n++;
      }
    do{
      boundary({s});
    }while(adapt_wavelet({s}, (double[]){se}, 19).nf);
    an = statsf(s).sum;
    se /= 2.;
  }while(fabs(an-ao) > pow(10, -DIGITS)*an);
  printf("Estimate:   %.10g\nAnalytical: %.10g\nReduction factor: %g\n", 
         an, pow (pi, 0.5), (1 << depth())/(double)n);
}

/**
The result may look like:

~~~bash
Estimate:   1.772460639
Analytical: 1.772453851
Reduction factor: 27.5246
~~~

Meaning that the function worked well and profited from grid adaptivity. 

*/
    
