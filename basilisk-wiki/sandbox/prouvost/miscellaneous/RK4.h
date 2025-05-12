
/**
Runge-Kutta 4th order

Integration between $t$ and $t+dt$ of the state vector $Y$ using the function $f\_rk4$ : $\frac{\text{d} Y}{\text{d} t} = f\_rk4(t,Y)$.
*/

// fixme: rename the structure my_vec into state_vec ?
typedef struct { double y0, y1;}  my_vec;   // vector (y,y') for RK4 integration
// why typedef struct and not just struct ? 
// It allows to wirte   my_vec Y;   instead of  struct my_vec Y;

my_vec RK4 (double t, double dt, my_vec Y, my_vec (* f_rk4) (double t, my_vec Y) ){

  my_vec k1 = f_rk4(t, Y);
  my_vec Y2;
  Y2.y0 = Y.y0+dt/2.*k1.y0;
  Y2.y1 = Y.y1+dt/2.*k1.y1;
  my_vec k2 = f_rk4(t+dt/2., Y2);
  Y2.y0 = Y.y0+dt/2.*k2.y0;
  Y2.y1 = Y.y1+dt/2.*k2.y1;
  my_vec k3 = f_rk4(t+dt/2., Y2);
  Y2.y0 = Y.y0+dt*k3.y0;
  Y2.y1 = Y.y1+dt*k3.y1;
  my_vec k4 = f_rk4(t+dt/2., Y2);

  Y2.y0 = Y.y0 + dt/6.*(k1.y0 + 2.*k2.y0 + 2.*k3.y0 + k4.y0);
  Y2.y1 = Y.y1 + dt/6.*(k1.y1 + 2.*k2.y1 + 2.*k3.y1 + k4.y1);
  
  return Y2;
}



/**
Euler 1st order

*/

my_vec Euler (double t, double dt, my_vec Y, my_vec (* f_rk4) (double t, my_vec Y) ){

  my_vec k1 = f_rk4(t, Y);
  my_vec Y2;
  Y2.y0 = Y.y0+dt*k1.y0;
  Y2.y1 = Y.y1+dt*k1.y1;
    
  return Y2;
}















