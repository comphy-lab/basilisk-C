/**
#Shallow Water Test Case with WENO 3rd order reconstruction - Le Metayer Test Case
*/

#include "grid/multigrid1D.h"
#include "saint-venant_modified.h"

void Compute_Analytical(double time){
 
  scalar u_analytical[], h_analytical[];
  double h_star,h_plus,h_minus,sigma,u_star;
  
  h_plus = 1.;
  h_minus = 1.8;
  h_star = 1.36898;  
  u_star = 2.*(sqrt(G*h_minus) - sqrt(G*h_star));
  sigma = (h_star*u_star)/(h_star - h_plus);

  foreach(){

    if(x < -1.*time*sqrt(G*h_minus)){
        u_analytical[] = 0.;
        h_analytical[] = h_minus;
      }
    else if( x > -1.*time*sqrt(G*h_minus) && x < time*(u_star - sqrt(G*h_star)) ){
        u_analytical[] = (2.*sqrt(G*h_minus) + 2.*x/time)/3.;
        h_analytical[] = sq(2.*sqrt(G*h_minus) - x/time)/(9.*G);
      }
    else if(x > time*(u_star - sqrt(G*h_star)) && x < time*sigma ){
        u_analytical[] = u_star;
        h_analytical[] = h_star;  
      }
    else if(x > time*sigma){
        u_analytical[] = 0.;
        h_analytical[] = h_plus;
      }
    }

   FILE *fp1 = fopen("Analytical_Solution.dat","w");
   FILE *fp3 = fopen("Error.dat","w");
   foreach(){
     fprintf(fp1,"%g %g %g \n",x,h_analytical[],u_analytical[]);
     fprintf(fp3,"%g %g %g \n",x,h[]-h_analytical[],u.x[]-u_analytical[]);
     }
   fclose(fp1);
   fclose(fp3);
}

int main()
{
  X0 = -300.;
  L0 = 600.;
  G = 10.;
  N = 256;
  run();
}

event init (i = 0)
{
  foreach()
    h[] = x < 0. ? 1.8 : 1.;
}

event output (t = 48) {
  FILE *fp2 = fopen("Numerical_Solution.dat","w");
  foreach()
    fprintf (fp2, "%g %g %g\n", x, h[], u.x[]);
  fclose(fp2);
  Compute_Analytical(t);
}

/**
## Results

~~~gnuplot Height
set output 'Heights.png'
set xlabel 'X'
set ylabel 'H(X)'
set grid
plot 'Numerical_Solution.dat' u 1:2 w p t 'Numerical Height', 'Analytical_Solution.dat' u 1:2 w p t 'Analytical Heights' 
~~~

~~~gnuplot Velocity
set output 'Velocity.png'
set xlabel 'X'
set ylabel 'U(X)'
set grid
plot 'Numerical_Solution.dat' u 1:3 w p t 'Numerical Velocity', 'Analytical_Solution.dat' u 1:3 w p t 'Analytical Velocity'
~~~

~~~gnuplot Error in Heights
set output 'Error_Heights.png'
set xlabel 'X'
set ylabel 'Error(H(X))'
set grid
plot 'Error.dat' u 1:2 w p t 'Error in Numerical Heights'
~~~

~~~gnuplot Error in Velocity
set output 'Error_Velocity.png'
set xlabel 'X'
set ylabel 'Error(U(X))'
set grid
plot 'Error.dat' u 1:3 w p t 'Error in Numerical Velocity'
~~~

*/

