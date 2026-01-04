/**
# Ising toy model

 Add description 

# Analytical formula
 In the case of no neighbour coupling (J = 0), an analytical formula for the average magnetization of the system exists as a function of the magnetic moment of a spin.

 $$\langle M \rangle = \mu N \text{tanh}\left(\frac{\mu H}{k_bT}\right)$$
 */

#include "grid/multigrid.h"
#include "run.h"
#include "utils.h"
#include "output.h"

double mcstep2d(scalar f, double J, double kbt, double H, double mu, int N){
  
  boundary({f}); //sync periodic ghosts
  
  int posx = ((rand() % (N)));
  int posy = ((rand() % (N)));
  double deltaE = 0.;
  double sumf = 0.;
  
  /**
# Horrible for performance. 
Why can I not just access a scalar field like f[i, j]?
*/
  foreach(serial){
    if (point.i - 1== posx && point.j - 1 == posy){
      int spin_curr = f[];
      int spin_above = f[0, 1];
      int spin_below = f[0, -1];
      int spin_right = f[1];
      int spin_left = f[-1];
      deltaE = 2. * J * spin_curr * (spin_right + spin_left + spin_above + spin_below) + 2 * spin_curr * H * mu;
      //fprintf(stderr, "posx posy %g %g %d %d\n %d %d %d %d %d\n", x, y, posx, posy, spin_curr, spin_left, spin_right, spin_above, spin_below);
      // fprintf(stderr, "deltae %g\n", deltaE);
      if (deltaE < 0.){
        f[] *= -1.;
      }
      else{
        double weight = exp(-deltaE/kbt);
        double rnd = rand()/(double) RAND_MAX;
        // fprintf(stderr, "rand %g %g %g %g\n", rnd, weight, deltaE, kbt);
        if (weight >= rnd){
          f[] *= -1.;
        }
      }
    }
    if (f[] != nodata) sumf += f[];
  } 
  
  return sumf;
}

double J, kbt, mu, H;
double mag, avg_mag;
double Niter;
int main(){
  J = 0.;
  kbt = 1.;
  H = 1.;

  Niter = 4e5;
  N = 64;
  periodic(right);
  periodic(left);
  for (mu = -2.; mu <= 2.; mu += 0.4){
   avg_mag = mag = 0.;
    run();
  }
}

scalar ising[], ising_init[];
event init(t = 0){
  double sum_init = 0.;
  foreach(reduction(+:sum_init)){
    double rand = noise();
    ising_init[] = rand > 0 ? ceil(rand) : floor(rand);
    sum_init += ising_init[];
    
  }
  ising = ising_init;
}

event spin_flip(i++){
  mag = mcstep2d(ising, J, kbt, H, mu, N);
  avg_mag += mag;
}


event logfile(i += 1000){
  double sum_after_step = 0.;
  foreach(reduction(+:sum_after_step)) 
  sum_after_step += ising[];
  fprintf(stderr, "i, mu, sum_after_step sum %d %g %g %g\n", i, mu, sum_after_step, avg_mag);

  FILE *fpmag = fopen("mag_time.dat", "a");
  fprintf(fpmag, "%d %g %g\n", i, mag, avg_mag);
  fclose(fpmag);
}

event movie(i += 1000., i <= Niter){
  output_ppm(ising, file = "ising.mp4", linear = true);    
}

event profile(t = end){
  dump();

  FILE *fp = fopen("avg_mag.dat", "a");
  fprintf(fp, "%g %g\n", mu, mu*avg_mag/iter/sq(N));
}


/**
## Results

Ok but not very good
~~~gnuplot Average magnetization as a function of mu
set xrange[-2.1:2.1]
set yrange[-0.1:2.1]

set xlabel 'mu'
set ylabel '<M>'

f(x) = x * tanh(x)
plot 'avg_mag.dat' title 'Sim. Data', \
     f(x) title 'Analytical'
~~~

Movie with N = 512, J = 0.43, mu = 1, H = 0, kbt = 0.25

![Movie with N = 512, J = 0.43, mu = 1, H = 0, kbt = 0.25](ising.mp4)(width="774" height="468")
*/

