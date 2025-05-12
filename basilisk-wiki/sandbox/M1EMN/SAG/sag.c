/**
SAG equation

$$
\frac{\partial c}{\partial t} =\nabla^2 c 
$$
f8701e3f07

*/
#include "run.h"

scalar c[];

event init(i=0) {
  double L = 0.1;
  foreach()
    c[] = exp(-(sq(x)+sq(y))/L);
}

event integration (i++){
  scalar dc[];
  double min = 1.0;
  
  foreach(){
    dc[]= (c[1,0] + c[0,1] + c[-1,0] + c[0,-1] - 4*c[])/sq(Delta); 
    if(Delta<min)
      min = Delta;
  }
  
  double dt  = sq(min)/2.;
    
  foreach()
    c[]+= dt *dc[];
  boundary({c});
}    

event graph ( t = 10){
  output_ppm(c,file="start1.png", n = 200);
}
/**
![init](/sandbox/M1EMN/SAG/sag/start1.png)
*/

int main () {
  X0 = -0.5; 
  run();

}
