/**
# Dynamical Update of Volume Fractions in Ghost Cells for a bubble in expansion

We recall that for a bubble in expansion, the normal (directed from the liquid 
phase towards the gas one) reads in axisymmetric coordinates:

$$
n_z = - \dfrac{\cos \theta}{|\cos \theta| + |\sin \theta|}
\quad ; \quad
n_r = - \dfrac{\sin \theta}{|\cos \theta| + |\sin \theta|}
$$
where $\theta = \arctan(r/z)$ and using the $L_1-$norm convention.


![Update of the ghost cells'volume fraction on the **left** border](img_rd_CA_bubble/fig_f_BC_matcha.png)(height=500 width=500) 

See the algorithm documentation
[here](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/f_BC_vel_tan.h) 
for updating ghost cell values.
*/

double f_BC_left(double x, double y, double Delta, 
                 double ff0, double ff1, double ff2){
  // ff0, ff1, ff2 = f[], f[0,1], f[0,-1] respectively
  if ((ff0 > 0.0) && (ff0 < 1.0)){
    double th = atan2(y,x);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff0, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-1.5,-0.5}, (coord){-0.5,0.5});
  }

  if ((ff1 > 0.0) && (ff1 < 1.0)){
    double th = atan2(y+Delta,x);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff1, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-1.5,-1.5}, (coord){-0.5,-0.5});
  }

/**
The following condition is the one schematized on the opening figure: 
*/

  if ((ff2 > 0.0) && (ff2 < 1.0)){
    double th = atan2(y-Delta,x);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff2, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-1.5,0.5}, (coord){-0.5,1.5});
  }

  return ff0;
}

/** 
## Analogous treatment for other boundaries

*/


double f_BC_bottom(double x, double y, double Delta, 
                   double ff0, double ff1, double ff2){
  // ff0, ff1, ff2 = f[], f[1,0], f[-1,0] respectively
  if ((ff0 > 0.0) && (ff0 < 1.0)){
    double th = atan2(y,x);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff0, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-0.5,-1.5}, (coord){0.5,-0.5});
  }

  if ((ff1 > 0.0) && (ff1 < 1.0)){
    double th = atan2(y,x+Delta);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff1, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-1.5,-1.5}, (coord){-0.5,-0.5});
  }

  if ((ff2 > 0.0) && (ff2 < 1.0)){
    double th = atan2(y,x-Delta);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff2, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){0.5,-1.5}, (coord){1.5,-0.5});
  }

  return ff0;
}

double f_BC_right(double x, double y, double Delta, 
                  double ff0, double ff1, double ff2){
  // ff0, ff1, ff2 = f[], f[0,1], f[0,-1] respectively
  if ((ff0 > 0.0) && (ff0 < 1.0)){
    double th = atan2(y,x);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff0, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){0.5,-0.5}, (coord){1.5,0.5});
  }

  if ((ff1 > 0.0) && (ff1 < 1.0)){
    double th = atan2(y+Delta,x);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff1, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){0.5,-1.5}, (coord){1.5,-0.5});
  }

  if ((ff2 > 0.0) && (ff2 < 1.0)){
    double th = atan2(y-Delta,x);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff2, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){0.5,0.5}, (coord){1.5,1.5});
  }

  return ff0;
}

double f_BC_top(double x, double y, double Delta, 
                double ff0, double ff1, double ff2){
  // ff0, ff1, ff2 = f[], f[1,0], f[-1,0] respectively
  if ((ff0 > 0.0) && (ff0 < 1.0)){
    double th = atan2(y,x);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff0, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-0.5,0.5}, (coord){0.5,1.5});
  }

  if ((ff1 > 0.0) && (ff1 < 1.0)){
    double th = atan2(y,x+Delta);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff1, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-1.5,0.5}, (coord){-0.5,1.5});
  }

  if ((ff2 > 0.0) && (ff2 < 1.0)){
    double th = atan2(y,x-Delta);
    coord ntemp = {
      - cos(th)/(fabs(- cos(th)) + fabs(- sin(th))), 
      - sin(th)/(fabs(- cos(th)) + fabs(- sin(th)))
    };
    double alpha_temp = plane_alpha (ff2, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){0.5,0.5}, (coord){1.5,1.5});
  }

  return ff0;
}