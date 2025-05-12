/** 
See the algorithm documentation
[here](http://basilisk.fr/sandbox/cailler/test-cases/vel_tan/f_BC_vel_tan.h) 
for updating ghost cell values.
*/

double f_BC_right(double ff0, double ff1, double ff2){
  // ff0, ff1, ff2 = f[], f[0,1], f[0,-1] respectively
  if ((ff0 > 0.0) && (ff0 < 1.0)){
    coord ntemp = {
      -sin(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0))), 
      cos(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)))
    };
    double alpha_temp = plane_alpha (ff0, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){0.5,-0.5}, (coord){1.5,0.5});
  }

  if ((ff1 > 0.0) && (ff1 < 1.0)){
    coord ntemp = {
      -sin(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0))), 
      cos(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)))
    };
    double alpha_temp = plane_alpha (ff1, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){0.5,-1.5}, (coord){1.5,-0.5});
  }

  if ((ff2 > 0.0) && (ff2 < 1.0)){
    coord ntemp = {
      -sin(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0))), 
      cos(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)))
    };
    double alpha_temp = plane_alpha (ff2, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){0.5,0.5}, (coord){1.5,1.5});
  }

  return ff0;
}

double f_BC_top(double ff0, double ff1, double ff2){
  // ff0, ff1, ff2 = f[], f[1,0], f[-1,0] respectively
  if ((ff0 > 0.0) && (ff0 < 1.0)){
    coord ntemp = {
      -sin(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0))), 
      cos(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)))
    };
    double alpha_temp = plane_alpha (ff0, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-0.5,0.5}, (coord){0.5,1.5});
  }

  if ((ff1 > 0.0) && (ff1 < 1.0)){
    coord ntemp = {
      -sin(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0))), 
      cos(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)))
    };
    double alpha_temp = plane_alpha (ff1, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){-1.5,0.5}, (coord){-0.5,1.5});
  }

  if ((ff2 > 0.0) && (ff2 < 1.0)){
    coord ntemp = {
      -sin(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0))), 
      cos(THETA_0)/(fabs(sin(THETA_0)) + fabs(cos(THETA_0)))
    };
    double alpha_temp = plane_alpha (ff2, ntemp);
    return rectangle_fraction (ntemp, alpha_temp, 
      (coord){0.5,0.5}, (coord){1.5,1.5});
  }

  return ff0;
}