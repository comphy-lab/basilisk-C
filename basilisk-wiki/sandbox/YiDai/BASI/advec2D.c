/**
# 2D advection equation
*/
#include "grid/cartesian.h"
#include "run.h"

scalar T[], dTx[], dTy[]; 
double dt;
double U1 = 3, U2 = 3;
#define EPS 1.

int main(){
    L0 = 10;
    X0 = -L0/2.;
    Y0 = -L0/2.;
    N = 1 << 8;
    DT = L0/N/U1/4;
    run();
}

// T[left] = dirichlet(3.);
// T[top] = dirichlet(0.);

event init(t = 0){
    foreach(){
        T[] = 10.*(sqrt(sq(x)+sq(y))<2);
        // T[] = 10 * sin(sqrt(sq(x) + sq(y))); 
    }
    boundary ({T});
}

event printdata (t = 0; t <= 0.3; t += 0.01) {
  static FILE * fp = fopen ("case0.dat","w");
  for (double y = -L0/2; y<L0/2; y+=L0/N){
    fprintf (fp, "%g %g\n", y, interpolate(T, 0, y));}
  fprintf (fp, "\n \n");
  fflush (fp);
}

event integration (i++) {
  double alpha = 1;    // uniform material
  double dt = DT;
  dt = dtnext (dt);
  foreach(){
    dTx[] = -U1 * (T[1,0] - T[-1,0])/(2 * Delta);
    dTy[] = -U2 * (T[0,1] - T[0,-1])/(2 * Delta);
    }
  foreach()
    T[] = 0.25* (T[1,0] + T[-1,0] + T[0,1] + T[0,-1]) + dt*alpha*(dTx[]+dTy[]);
  boundary ({T});
}

event movie(t = 0; t <= 0.3; t += 0.01) {
  output_ppm (T, file = "temp_case0.mp4", min = 0, max = 5, linear = true, map = cool_warm);
}
/**
![](advec2D/temp_case0.mp4)

*/