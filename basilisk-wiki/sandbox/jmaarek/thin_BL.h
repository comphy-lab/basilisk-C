#define BGHOSTS 2

#ifdef interface_normal
#undef interface_normal
#endif

coord interface_normal (Point point, scalar c);
#define interface_normal(point, c) interface_normal (point, c)

#include "fractions.h"
#include "curvature.h"

attribute {
  double c_inf;
  scalar delta_b;
  scalar c_boundary; // private
}

/*to achieve second order convergence on surface area we need a second order approximation of the normal vector.
We do that by using a normal computed by height functions */

coord interface_normal(Point point, scalar c){
  coord n;
  if (!c.height.x.i || (n = height_normal (point, c, c.height)).x == nodata)
    n = mycs (point, c);
  else {
    double nn = 0.;
    foreach_dimension()
      nn += fabs(n.x);
    foreach_dimension()
      n.x /= nn;
  }

  if ((fabs(n.x) + fabs(n.y) + fabs(n.z)) == 0.0){
    foreach_dimension()
      n.x = 1.0;
    double nn = 0.;
    foreach_dimension()
      nn += fabs(n.x);
    foreach_dimension()
      n.x /= nn;
  }  

  return n;
}

/***

compared to the Bothe paper, we choose to define a vof tracer defined in the explicitly defined phase as opposed to the implicit phase
This is consistent with standard tracer methods in basilisk. At the current moment this is done by changing the sign of alpha in the plane transformation

***/
//Function converts planar representation of interface described in geometry.h to description used in Bothe paper
//signed distance origin in the edge of the disperse phase farthest away from the interface
void plane_transformation(const coord n_in, const double alpha_in, coord * n_out, double * alpha_out){
  *n_out = n_in;
  /*changing sign inverts the plane such that cut cell is now describing
  the implcit phase. As such we can leave all other functions the same and
  now we are tracking a tracer defined in the explicit phase given the volume
  and surface integrals are written for implicit phase*/

  *alpha_out = -alpha_in;

  //translation
  foreach_dimension()
    *alpha_out += 0.5*n_in.x;

  //reflection
  foreach_dimension()
    if (n_in.x < 0.0)
      *alpha_out -= n_in.x;

  foreach_dimension()
    (*n_out).x = fabs((*n_out).x);

  //sort normals
  double nmax = -HUGE;
  double nmin = HUGE;
  double nmid = 0.0;

  foreach_dimension(){
    if ((*n_out).x > nmax){
      nmax = (*n_out).x;}
    if ((*n_out).x < nmin){
      nmin = (*n_out).x;}
    if ((((*n_out).x >= (*n_out).y) && ((*n_out).x <= (*n_out).z)) || (((*n_out).x >= (*n_out).z) && ((*n_out).x <= (*n_out).y))){
      nmid = (*n_out).x;}
  }
  *n_out = (coord){nmax,nmid,nmin};

  //scale unit normal
  double nn = 0.0;
  foreach_dimension()
    nn += sq((*n_out).x);
  nn = sqrt(nn);

  foreach_dimension()
    (*n_out).x /= nn;
  *alpha_out /= nn;
}

//Volume integrals
double f_xyz(coord a, coord n, double alpha, double delta_b){
  return (a.x*n.x + a.y*n.y + a.z*n.z - alpha)/delta_b;}


double f_term(coord b, coord n, double alpha, double delta_b){
    return (sq(f_xyz(b,n,alpha,delta_b))+1.5)*f_xyz(b,n,alpha,delta_b)*erf(f_xyz(b,n,alpha,delta_b))
    + (sq(f_xyz(b,n,alpha,delta_b))+1)*exp(-sq(f_xyz(b,n,alpha,delta_b)))/sqrt(M_PI);}

double f_term2(coord b1, coord b2, coord n, double alpha, double delta_b){
    return  f_xyz(b1,n,alpha,delta_b)*erf(f_xyz(b1,n,alpha,delta_b))
           -f_xyz(b2,n,alpha,delta_b)*erf(f_xyz(b2,n,alpha,delta_b))
           + (exp(-sq(f_xyz(b1,n,alpha,delta_b))) - exp(-sq(f_xyz(b2,n,alpha,delta_b))))/sqrt(M_PI);}

double gamma_cube(coord a, coord b, coord n, double alpha, double delta_b){
  return cube(delta_b)/(6*n.x*n.y*n.z)*(
     f_term(b, n, alpha, delta_b)
     - f_term((coord){b.x,b.y,a.z}, n, alpha, delta_b)
     - f_term((coord){b.x,a.y,b.z}, n, alpha, delta_b)
     + f_term((coord){b.x,a.y,a.z}, n, alpha, delta_b)
     - f_term((coord){a.x,b.y,b.z}, n, alpha, delta_b)
     + f_term((coord){a.x,b.y,a.z}, n, alpha, delta_b)
     + f_term((coord){a.x,a.y,b.z}, n, alpha, delta_b)
     - f_term(a, n, alpha, delta_b));}


double Dgamma_cube_Ddelta_b(coord a, coord b, coord n, double alpha, double delta_b){
  return sq(delta_b)/(2*n.x*n.y*n.z)*(
      f_term2(b,  (coord){b.x,b.y,a.z}, n, alpha, delta_b)
      -f_term2((coord){b.x,a.y,b.z}, (coord){b.x,a.y,a.z}, n, alpha, delta_b)
      -f_term2((coord){a.x,b.y,b.z}, (coord){a.x,b.y,a.z}, n, alpha, delta_b)
      +f_term2((coord){a.x,a.y,b.z}, a, n, alpha, delta_b));}


/*
double gamma_wedge(coord a, coord b, coord n, double alpha, double delta_b){
  double K_x = n.x*(b.x-a.x)/delta_b;
  return cube(delta_b)/(6*n.x*n.y*n.z)*(
    3*n.z*(b.z-a.z)/delta_b*((sq(K_x)+0.5)*erf(K_x)+K_x*exp(-sq(K_x))/sqrt(M_PI))
    - f_term((coord){b.x,a.y,b.z}, n, alpha, delta_b)
    + f_term((coord){b.x,a.y,a.z}, n, alpha, delta_b)
    + f_term((coord){a.x,a.y,b.z}, n, alpha, delta_b)
    - f_term(a, n, alpha, delta_b));}


double Dgamma_wedge_Ddelta_b(coord a, coord b, coord n, double alpha, double delta_b){
  double K_x = n.x*(b.x-a.x)/delta_b;
  return sq(delta_b)/(2*n.x*n.y*n.z)*(
    n.z*(b.z-a.z)/delta_b*erf(K_x)
    -f_xyz((coord){b.x,a.y,b.z},n,alpha,delta_b)*erf(f_xyz((coord){b.x,a.y,b.z},n,alpha,delta_b))
       + f_xyz((coord){b.x,a.y,a.z},n,alpha,delta_b)*erf(f_xyz((coord){b.x,a.y,a.z},n,alpha,delta_b))
       - (exp(-sq(f_xyz((coord){b.x,a.y,b.z},n,alpha,delta_b))) + exp(-sq(f_xyz((coord){b.x,a.y,a.z},n,alpha,delta_b))))/sqrt(M_PI)
       + f_xyz((coord){a.x,a.y,b.z},n,alpha,delta_b)*erf(f_xyz((coord){a.x,a.y,b.z},n,alpha,delta_b))
         - f_xyz(a,n,alpha,delta_b)*erf(f_xyz(a,n,alpha,delta_b))
         + (exp(-sq(f_xyz((coord){a.x,a.y,b.z},n,alpha,delta_b))) - exp(-sq(f_xyz(a,n,alpha,delta_b))))/sqrt(M_PI));}*/

double gamma_tet(coord a, coord n, double alpha, double delta_b){
    double bz = (alpha - a.x*n.x - a.y*n.y)/n.z;
    return cube(delta_b)/(6*n.x*n.y*n.z)*(
      6*n.z*(bz-a.z)/(sqrt(M_PI)*sq(delta_b))*(alpha - a.x*n.x - a.y*n.y - n.z*(bz+a.z)/2)
      + f_term((coord){a.x,a.y,bz}, n, alpha, delta_b)
      - f_term(a, n, alpha, delta_b));}

double Dgamma_tet_Ddelta_b(coord a, coord n, double alpha, double delta_b){
    double bz = (alpha - a.x*n.x - a.y*n.y)/n.z;
    return sq(delta_b)/(2*n.x*n.y*n.z)*(
      2*n.z*(bz-a.z)/(sqrt(M_PI)*sq(delta_b))*(alpha - a.x*n.x - a.y*n.y - n.z*(bz+a.z)/2)
      + f_term2((coord){a.x,a.y,bz}, a, n, alpha, delta_b));}


//compute SGS model function in cell from VOF
double compute_eta_sgs (coord n_transform, double alpha_transform, double delta_b, double f, coord a, coord b){
  double eta_sgs = 0.0;
  if (f_xyz(a, n_transform, alpha_transform, delta_b) < 0.0){
    if (f_xyz(b, n_transform, alpha_transform, delta_b) < 0.0){
      return 0.0;
    }
    else{

      eta_sgs +=  gamma_cube(a, b, n_transform, alpha_transform, delta_b)
                  - gamma_tet(a, n_transform, alpha_transform, delta_b);
      foreach_dimension(){
          if (alpha_transform - b.x*n_transform.x - a.y*n_transform.y - a.z*n_transform.z > 0.0){
            coord a_1 = a;
            a_1.x = b.x;
            eta_sgs += gamma_tet(a_1, n_transform, alpha_transform, delta_b);}

          if (alpha_transform - b.x*n_transform.x - b.y*n_transform.y - a.z*n_transform.z > 0.0){
            coord a_2 = a;
            a_2.x = b.x;
            a_2.y = b.y;
            eta_sgs -= gamma_tet(a_2, n_transform, alpha_transform, delta_b);}
          }
    }
  }
  else{
    eta_sgs = gamma_cube(a, b, n_transform, alpha_transform, delta_b);
  }
  return eta_sgs/f;
}

double compute_grad_eta_sgs (coord n_transform, double alpha_transform, double delta_b, double f, coord a, coord b){
  double eta_sgs = 0.0;
  if (f_xyz(a, n_transform, alpha_transform, delta_b) < 0.0){
    if (f_xyz(b, n_transform, alpha_transform, delta_b) < 0.0){
      return 0.0;
    }
    else{

      eta_sgs +=  Dgamma_cube_Ddelta_b(a, b, n_transform, alpha_transform, delta_b)
                  - Dgamma_tet_Ddelta_b(a, n_transform, alpha_transform, delta_b);
      foreach_dimension(){
          if (alpha_transform - b.x*n_transform.x - a.y*n_transform.y - a.z*n_transform.z > 0.0){
            coord a_1 = a;
            a_1.x = b.x;
            eta_sgs += Dgamma_tet_Ddelta_b(a_1, n_transform, alpha_transform, delta_b);}

          if (alpha_transform - b.x*n_transform.x - b.y*n_transform.y - a.z*n_transform.z > 0.0){
            coord a_2 = a;
            a_2.x = b.x;
            a_2.y = b.y;
            eta_sgs -= Dgamma_tet_Ddelta_b(a_2, n_transform, alpha_transform, delta_b);}
          }
    }
  }
  else{
    eta_sgs = Dgamma_cube_Ddelta_b(a, b, n_transform, alpha_transform, delta_b);
  }
  return eta_sgs/f;
}

//fit the boundary layer thickness to the concentration in the cell (cube)
double Bisection_Newton(coord n_transform, double alpha_transform, double f, double cc, double c_boundary, double c_inf, double Deltax, double delta_0){

  double delta_min = 1e-5;
  double delta_max = 1000;

  delta_0 = clamp((delta_0 > 0.0 ? delta_0/Deltax : 1.0), delta_min*10, delta_max/10);


  coord a = {0.0,0.0,0.0};
  coord b = {1.0,1.0,1.0};

  double res_max = compute_eta_sgs(n_transform, alpha_transform, delta_max, f, a,b)  - (cc/f-c_boundary)/(c_inf-c_boundary);
  double delta_n = delta_0;
  double res = compute_eta_sgs(n_transform, alpha_transform, delta_n, f, a,b)  - (cc/f-c_boundary)/(c_inf-c_boundary);
  double delta_n1;
  int ii = 0;


  //CASE IN WHICH THE NEWTON'S ALGORITHM WILL CRASH
  if ((cc/f-c_boundary)/(c_inf-c_boundary) <= 0.0)
    return -1;

  if ((cc/f-c_boundary)/(c_inf-c_boundary) >= 1.0)
    return -1;

  if ( (cc/f-c_boundary)/(c_inf-c_boundary) - compute_eta_sgs(n_transform, alpha_transform, delta_max, f, a,b) <= 0.0 )
	return delta_max*Deltax;

  if ( (cc/f-c_boundary)/(c_inf-c_boundary) - compute_eta_sgs(n_transform, alpha_transform, delta_min, f, a,b) >= 0.0 )
        return delta_min*Deltax;

  //while(fabs((res*(c_inf-c_boundary)+c_boundary) - cc/f)/(cc/f) > 1e-9){
  while(fabs(res/((cc/f-c_boundary)/(c_inf-c_boundary))) > 1e-9){
    double grad = compute_grad_eta_sgs(n_transform, alpha_transform, delta_n, f, a,b);
    if (grad == 0.0)
      return -1;
    delta_n1 = delta_n - res/grad;
    if ((delta_n1 < delta_min) || (delta_n1 > delta_max)){
      if (res*res_max > 0.0){
        delta_n1 = (delta_min+delta_n)/2;
        delta_max = delta_n;
      }
      else{
        delta_n1 = (delta_max+delta_n)/2;
        delta_max = delta_n;
      }
    }
  delta_n = delta_n1;
  res = compute_eta_sgs(n_transform, alpha_transform, delta_n, f, a,b)  - (cc/f-c_boundary)/(c_inf-c_boundary);
  ii++;


  //in cases where the boundary layer thickness is much greater than the grid size (>100) or the vof fraction is very small the algorithm might struggle to converge
  //in these cases the numerical flux is usually sufficient and flux correction based on boundary layer modeling is unecessary

  //fprintf(ferr, "%f %d\n", delta_n, ii);

  if (ii > 30){
    return -1.0;
  }
}
  return delta_n*Deltax;
}

double Bisection_Newton_2(coord n_transform, double alpha_transform, double f, double cc, double c_boundary, double c_inf, double Deltax, double delta_0, coord b){

  double delta_min = 1e-5;
  double delta_max = 1000;

  delta_0 = clamp((delta_0 > 0.0 ? delta_0/Deltax : 1.0), delta_min*10, delta_max/10);


  coord a = {0.0,0.0,0.0};

  double res_max = compute_eta_sgs(n_transform, alpha_transform, delta_max, f, a,b) - (cc/f-c_boundary)/(c_inf-c_boundary);
  double delta_n = delta_0;
  double res = compute_eta_sgs(n_transform, alpha_transform, delta_n, f, a,b) - (cc/f-c_boundary)/(c_inf-c_boundary);
  double delta_n1;
  int ii = 0;


  //CASE IN WHICH THE NEWTON'S ALGORITHM WILL CRASH
  if ((cc/f-c_boundary)/(c_inf-c_boundary) <= 0.0)
    return -1;

  if ((cc/f-c_boundary)/(c_inf-c_boundary) >= 1.0)
    return -1;

  if ( (cc/f-c_boundary)/(c_inf-c_boundary) - compute_eta_sgs(n_transform, alpha_transform, delta_max, f, a,b) <= 0.0 )
        return delta_max*Deltax;

  if ( (cc/f-c_boundary)/(c_inf-c_boundary) - compute_eta_sgs(n_transform, alpha_transform, delta_min, f, a,b) >= 0.0 )
        return delta_min*Deltax;

  //while(fabs((res*(c_inf-c_boundary)+c_boundary) - cc/f)/(cc/f) > 1e-9){
  while(fabs(res/((cc/f-c_boundary)/(c_inf-c_boundary))) > 1e-9){
    double grad = compute_grad_eta_sgs(n_transform, alpha_transform, delta_n, f, a,b);
    if (grad == 0.0)
      return -1;
    delta_n1 = delta_n - res/grad;
    if ((delta_n1 < delta_min) || (delta_n1 > delta_max)){
      if (res*res_max > 0.0){
        delta_n1 = (delta_min+delta_n)/2;
        delta_max = delta_n;
      }
      else{
        delta_n1 = (delta_max+delta_n)/2;
        delta_max = delta_n;
      }
    }
  delta_n = delta_n1;
  res = compute_eta_sgs(n_transform, alpha_transform, delta_n, f, a,b) - (cc/f-c_boundary)/(c_inf-c_boundary);
  ii++;


  //in cases where the boundary layer thickness is much greater than the grid size (>100) or the vof fraction is very small the algorithm might struggle to converge
  //in these cases the numerical flux is usually sufficient and flux correction based on boundary layer modeling is unecessary

  //fprintf(ferr, "%f %d\n", delta_n, ii);

  if (ii > 30){
    return -1.0;
  }
}
  return delta_n*Deltax;
}



static inline void restriction_delta_b(Point point, scalar s){
  double sum = 0.;
  double count = 0.0;
  foreach_child(){
    sum += s[] > 0.0 ? s[] : 0.0;
    count += s[] > 0.0 ? 1.0 : 0.0;
  }
  s[] = count > 0.0 ? sum/count : -1;
}

//fit the boundary layer thickness to the concentration in the cell (cube)


//populates bl_tracer.delta_b with boundary layer thickness computed from the vof interface f and the concetration tracer c
//if the cell is a bulk cell or if the bisection newton solver does not converge populate with -1
void boundary_layer(scalar bl_tracer){

  scalar f = bl_tracer.c;
  scalar delta_b = bl_tracer.delta_b;
  scalar c_boundary = bl_tracer.c_boundary;
  double c_inf = bl_tracer.c_inf;

  delta_b.restriction = restriction_delta_b;
  delta_b.refine = delta_b.prolongation = refine_injection;

  //fprintf(fout, "Start bisection newton\n");

  foreach(){

    coord n;
    coord n_out;
    double alpha;
    double alpha_out;


    if ((f[] > 0.0) && (f[] < 1)){ //decides if a cell is interfacial
        delta_b[] = -1;
        n = interface_normal(point, f);

        foreach_dimension(){ // done to avoid issues that might occur if the magnitude of normal components are too small
          n.x = sign(n.x)*max(fabs(n.x), 1e-5);
        }

        //L1 normalization
        double summ = 0.0;
        foreach_dimension()
          summ += fabs(n.x);
        foreach_dimension()
          n.x /= summ;

        //fprintf(ferr, "test 1 \n");

        //plane recontruction and transformation.
        //Plane transformation is a local signed distance function defined as seen in bothe paper
        if (bl_tracer.inverse)
    		  alpha = plane_alpha (1-f[], n);
        else
          alpha = plane_alpha (f[], n);
        plane_transformation(n, alpha, &n_out, &alpha_out);

        if (bl_tracer.inverse){
          if ((1 - f[]) >= 0.2){
              delta_b[] = Bisection_Newton(n_out, alpha_out, 1 - f[], bl_tracer[], c_boundary[], c_inf, Delta, delta_b[]);
          }
          else{
              foreach_dimension() //cell linking
                if ((fabs(n.x)/sqrt(sq(n.x)+sq(n.y)+sq(n.z)) > 0.9) && (1 - f[sign(n.x)] >= 1.0))
                    delta_b[] = Bisection_Newton_2(n_out, alpha_out, 2 - f[], bl_tracer[] + bl_tracer[sign(n.x)], c_boundary[], c_inf, Delta, delta_b[], (coord){2.0,1.0,1.0});
          }
        }

        else{
          if (f[] >= 0.2){
            delta_b[] = Bisection_Newton(n_out, alpha_out, f[], bl_tracer[], c_boundary[], c_inf, Delta, delta_b[]);
          }
          else{
            foreach_dimension()
              if ((fabs(n.x)/sqrt(sq(n.x)+sq(n.y)+sq(n.z)) > 0.9) && (f[-sign(n.x)] >= 1.0)) //cell linking
                delta_b[] = Bisection_Newton_2(n_out, alpha_out, 1 + f[], bl_tracer[] + bl_tracer[-sign(n.x)], c_boundary[], c_inf, Delta, delta_b[], (coord){2.0,1.0,1.0});
        }
      }
    }
    else{
      delta_b[] = -1;
    }
 }
 //fprintf(fout, "Finish foreach\n");
 boundary({delta_b});
 restriction({delta_b});


 //fprintf(fout, "Start tag\n");

 //for interfacial cells that did not converge declare boundary layer thickness to the mean of the nearest neighbors


 //declare temporary scalar tagging values that can be used for running mean
 scalar tagged_cells[];
 tagged_cells.restriction = restriction_average;
 tagged_cells.refine = tagged_cells.prolongation = refine_injection;

 //tag cells that the fitting algorithm converged
 foreach(){
   if ((f[] > 0) && (f[] < 1) && (delta_b[] != -1.0)){
    tagged_cells[] = 1.0;}
   else{
    tagged_cells[] = 0.0;}
  }
  boundary({tagged_cells});

  //for interfacial cell that the fitting algorithm failed to converge, declare boundary layer thickness to the volume weighted mean of the nearest neighbors
  foreach(){
    if ((f[] > 0) && (f[] < 1) && (tagged_cells[] == 0.0) && ((bl_tracer[]/f[]-c_boundary[])/(c_inf-c_boundary[]) > 0.0) && ((bl_tracer[]/f[]-c_boundary[])/(c_inf-c_boundary[]) < 1.0)){
      double bl_average = 0.0;
      double count = 0.0;
      foreach_neighbor(1){
        if (tagged_cells[] > 0.0){
          bl_average += delta_b[];
          count += 1.0;
        }
      }
      if (count > 0){
        delta_b[] = bl_average/count;}
    }
  }
  boundary({delta_b});
  /*if the value of delta_b is still -1 that means that there were no nearest neighbors where the fitting method converged
  in this case we will treat the tracer as a normal vof tracer. percentage should be small. */
}


void boundary_layer2(scalar bl_tracer){

  scalar f = bl_tracer.c;
  scalar delta_b = bl_tracer.delta_b;
  scalar c_boundary = bl_tracer.c_boundary;
  double c_inf = bl_tracer.c_inf;

  delta_b.restriction = restriction_delta_b;
  delta_b.refine = delta_b.prolongation = refine_injection;

  //fprintf(fout, "Start bisection newton\n");

  foreach(){

    coord n;
    coord n_out;
    double alpha;
    double alpha_out;


    if ((f[] > 0.0) && (f[] < 1)){ //decides if a cell is interfacial
        delta_b[] = -1;
        n = interface_normal(point, f);

        foreach_dimension(){ // done to avoid issues that might occur if the magnitude of normal components are too small
          n.x = sign(n.x)*max(fabs(n.x), 1e-5);
        }

        //L1 normalization
        double summ = 0.0;
        foreach_dimension()
          summ += fabs(n.x);
        foreach_dimension()
          n.x /= summ;

        //fprintf(ferr, "test 1 \n");

        //plane recontruction and transformation.
        //Plane transformation is a local signed distance function defined as seen in bothe paper
        if (bl_tracer.inverse)
    		  alpha = plane_alpha (1 - f[], n);
        else
          alpha = plane_alpha (f[], n);
        plane_transformation(n, alpha, &n_out, &alpha_out);

        if (bl_tracer.inverse){
              foreach_dimension(){ //cell linking
                if ((fabs(n.x) > fabs(n.y)) &&  (fabs(n.x) > fabs(n.z))){
                  if ((1 - f[sign(n.x)]) >= 1.0){
                    delta_b[] = Bisection_Newton_2(n_out, alpha_out, 2 - f[], bl_tracer[] + bl_tracer[sign(n.x)], c_boundary[], c_inf, Delta, delta_b[], (coord){2.0,1.0,1.0});}
                  else{
                      if ((1 - f[]) >= 0.2){
                        delta_b[] = Bisection_Newton(n_out, alpha_out, 1 - f[], bl_tracer[], c_boundary[], c_inf, Delta, delta_b[]);}
                    }
                }
              }
        }

        else{
          foreach_dimension(){ //cell linking
            if ((fabs(n.x) >= fabs(n.y)) &&  (fabs(n.x) >= fabs(n.z))){
              if ((f[-sign(n.x)] >= 1.0)){
                delta_b[] = Bisection_Newton_2(n_out, alpha_out, f[] + 1, bl_tracer[] + bl_tracer[-sign(n.x)], c_boundary[], c_inf, Delta, delta_b[], (coord){2.0,1.0,1.0});}
              else{
                  if (f[] >= 0.2){
                    delta_b[] = Bisection_Newton(n_out, alpha_out, f[], bl_tracer[], c_boundary[], c_inf, Delta, delta_b[]);}
                }
            }
          }
      }
    }
    else{
      delta_b[] = -1;
    }
 }
 //fprintf(fout, "Finish foreach\n");
 boundary({delta_b});


 //fprintf(fout, "Start tag\n");

 //for interfacial cells that did not converge declare boundary layer thickness to the mean of the nearest neighbors


 //declare temporary scalar tagging values that can be used for running mean
 scalar tagged_cells[];
 tagged_cells.restriction = no_restriction;
 tagged_cells.refine = tagged_cells.prolongation = refine_injection;

 //tag cells that the fitting algorithm converged
 foreach(){
   if ((f[] > 0) && (f[] < 1) && (delta_b[] != -1.0)){
    tagged_cells[] = 1.0;}
   else{
    tagged_cells[] = 0.0;}
  }
  boundary({tagged_cells});

  //for interfacial cell that the fitting algorithm failed to converge, declare boundary layer thickness to the volume weighted mean of the nearest neighbors
  foreach(){
    if ((f[] > 0) && (f[] < 1) && (tagged_cells[] == 0.0) && ((bl_tracer[]/f[]-c_boundary[])/(c_inf-c_boundary[]) > 0.0) && ((bl_tracer[]/f[]-c_boundary[])/(c_inf-c_boundary[]) < 1.0)){
      double bl_average = 0.0;
      double count = 0.0;
      foreach_neighbor(1){
        if (tagged_cells[] > 0.0){
          bl_average += delta_b[];
          count += 1.0;
        }
      }
      if (count > 0){
        delta_b[] = bl_average/count;}
    }
  }
  boundary({delta_b});
}



/*for transporting and refining boundary layer tracer we use a stricter gradient and refining criteria
where both the cell and an adjacent cell must be full to use the numerical gradient otherwise it is zero.
This ensures that concentrations will always remain physical*/

foreach_dimension()
static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
{
  static const double cmin = 1.0;
  double cl = c[-1], cc = c[], cr = c[1];
  if (t.inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && t.gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
	if (t.gradient)
	  return t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
	else
	  return (t[1]/cr - t[-1]/cl)/(2.*Delta);
      }
      else
	return (t[1]/cr - t[]/cc)/Delta;
    }
    else if (cl >= cmin)
      return (t[]/cc - t[-1]/cl)/Delta;
  }
  return 0.;
}

static void vof_concentration_refine(Point point, scalar s)
{
  scalar f = s.c;
  if (cm[] == 0. || (!s.inverse && f[] <= 0.) || (s.inverse && f[] >= 1.))
    foreach_child()
      s[] = 0.;
  else {
    coord g;
    foreach_dimension()
      g.x = Delta*vof_concentration_gradient_x (point, f, s);
    double sc = s.inverse ? s[]/(1. - f[]) : s[]/f[], cmc = 4.*cm[];
    foreach_child() {
      s[] = sc;
      foreach_dimension()
	     s[] += child.x*g.x*cm[-child.x]/cmc;
      s[] *= s.inverse ? 1. - f[] : f[];
    }
  }
}



static void bl_tracer_refine(Point point, scalar s){
  scalar f = s.c;
  scalar g = s.delta_b;
  scalar h = s.c_boundary;

  coord n_out;
  double alpha_out;

  //redistribute the vof tracer in the child cells based on the model function and the boundary layer thickness
  if ((g[] > 0.0) && (f[] > 0) && (f[] < 1)){
    coord n = interface_normal(point, f);

    foreach_dimension(){ // done to avoid issues that might occur if the magnitude of normal compoents are too small
      n.x = sign(n.x)*max(fabs(n.x), 1e-5);
    }

    //L1 normalization
    double summ = 0.0;
    foreach_dimension()
      summ += fabs(n.x);
    foreach_dimension()
      n.x /= summ;

    double alpha;

    if (s.inverse){
      alpha = plane_alpha (1 - f[], n);}
    else{
      alpha = plane_alpha (f[], n);}

   //sort normals
   double n_temp[3] = {n.x, n.y, n.z};
   int idx[3] = {0,1,2};

    for (int ii=0; ii<3; ii++)
    {
      for (int jj=ii+1; jj<3; jj++)
      {
        if (fabs(n_temp[idx[ii]]) < fabs(n_temp[idx[jj]]))
        {
          swap (double, idx[ii], idx[jj]);
        }
      }
    }
  //fprintf(fp2, "%g %g %g %d %d %d\n", n.x, n.y, n.z, idx[0], idx[1], idx[2]);


    plane_transformation(n, alpha, &n_out, &alpha_out);
    //coord a = {0.0,0.0,0.0};
    //coord b = {1.0,1.0,1.0};
    //double temp = (compute_eta_sgs(n_out, alpha_out, g[]/Delta, f[], a,b)*(s.c_inf-1.0)+1.0)*f[];
    double cc_total = 0.0;
    for (int ii = 0; ii <= 1; ii++){
      for (int jj = 0; jj <= 1; jj++){
        for (int kk = 0; kk <= 1; kk++){
          double a_temp[3] = {0.5*(-sign(n.x)*(((double)ii)-0.5) + 0.5), 0.5*(-sign(n.y)*(((double)jj)-0.5) + 0.5), 0.5*(-sign(n.z)*(((double)kk)-0.5) + 0.5)};
          coord a = {a_temp[idx[0]], a_temp[idx[1]], a_temp[idx[2]]};
          coord b = {a.x+0.5, a.y + 0.5, a.z + 0.5};
          //fprintf(fp2, "%g %g %g %g %g %g %d %d %d %g %g %g\n", n.x, n.y, n.z, n_out.x, n_out.y, n_out.z, ii, jj, kk, a.x, a.y, a.z);
          double child_volume_fraction = rectangle_fraction(n, alpha, (coord){0.5*((double)ii)-0.5,0.5*((double)jj)-0.5,0.5*((double)kk)-0.5}, (coord){0.5*((double)ii),0.5*((double)jj),0.5*((double)kk)});
          if (child_volume_fraction > 0.0){

            //fprintf(fp2, "%g %g %g %g\n", child_volume_fraction, compute_eta_sgs(n_out, alpha_out, g[]/Delta, 0.125, a,b), f[], compute_eta_sgs(n_out, alpha_out, g[]/Delta, 1, (coord){0.0,0.0,0.0},(coord){1.0,1.0,1.0}));

            //compute_eta_sgs(n_out, alpha_out, g[]/Delta, 1/8, a,b);

            fine(s,ii,jj,kk) = clamp(compute_eta_sgs(n_out, alpha_out, g[]/Delta, child_volume_fraction/8, a,b),0.0,1.0);
            fine(s,ii,jj,kk) = (fine(s,ii,jj,kk)*(s.c_inf-h[])+h[])*child_volume_fraction/8;
            cc_total += fine(s,ii,jj,kk);
          }
          else{
            fine(s,ii,jj,kk) = 0.0;
          }
        }
      }
    }

    //normalize the tracer to conserve mass and then convert tracer magnitude to concentration
    for (int ii = 0; ii <= 1; ii++)
      for (int jj = 0; jj <= 1; jj++)
        for (int kk = 0; kk <= 1; kk++){
          //fprintf(fp2, "%g %g %g %g\n", cc_total, s[], fine(s,ii,jj,kk), f[]);
          fine(s,ii,jj,kk) = cc_total > 0.0 ? fine(s,ii,jj,kk)*s[]/cc_total*8 : 0.0;
        }
    //fprintf(fp2, "\n\n");
  }
  else{ //vof_concentration_refine
    vof_concentration_refine(point, s);
  }
}

//Functions to compute average diffusion flux on an interface

double zeta_rect(coord a, coord b, coord n, double alpha, double delta_b){
    return delta_b/(n.y*n.z)*(
      f_xyz((coord){0.0,b.y,b.z},n,alpha,delta_b)*erf(f_xyz((coord){0.0,b.y,b.z},n,alpha,delta_b))
    - f_xyz((coord){0.0,b.y,a.z},n,alpha,delta_b)*erf(f_xyz((coord){0.0,b.y,a.z},n,alpha,delta_b))
    - f_xyz((coord){0.0,a.y,b.z},n,alpha,delta_b)*erf(f_xyz((coord){0.0,a.y,b.z},n,alpha,delta_b))
    + f_xyz((coord){0.0,a.y,a.z},n,alpha,delta_b)*erf(f_xyz((coord){0.0,a.y,a.z},n,alpha,delta_b))
    + 1/sqrt(M_PI)*(exp(-sq(f_xyz((coord){0.0,b.y,b.z},n,alpha,delta_b)))
                  - exp(-sq(f_xyz((coord){0.0,b.y,a.z},n,alpha,delta_b)))
                  - exp(-sq(f_xyz((coord){0.0,a.y,b.z},n,alpha,delta_b)))
                  + exp(-sq(f_xyz((coord){0.0,a.y,a.z},n,alpha,delta_b)))));}


double zeta_tri(coord b, coord n, double alpha, double delta_b){
    double az = (alpha-n.y*b.y)/n.z;
    return delta_b/(n.y*n.z)*(
      f_xyz((coord){0.0,b.y,b.z},n,alpha,delta_b)*erf(f_xyz((coord){0.0,b.y,b.z},n,alpha,delta_b))
    - f_xyz((coord){0.0,b.y,az},n,alpha,delta_b)*erf(f_xyz((coord){0.0,b.y,az},n,alpha,delta_b))
    + 1/sqrt(M_PI)*(exp(-sq(f_xyz((coord){0.0,b.y,b.z},n,alpha,delta_b)))
                  - exp(-sq(f_xyz((coord){0.0,b.y,az},n,alpha,delta_b)))));}

//function computes average flux of model function on chosen interface given input face, cut cell data, and boundary_layer thickness

double compute_zeta_sgs (coord n_transform, double alpha_transform, double delta_b, coord plane){ //plane should look like {-1,-1,a} where positive value components is plane of interest and value is z = a plane
  double n_temp[2] = {0.0,0.0};
  int ii = 0;
  double alpha_2d = 0.0;
  foreach_dimension(){
    if (plane.x != -1){
      alpha_2d = alpha_transform - n_transform.x*plane.x;
    }
    else{
      n_temp[ii] = n_transform.x;
      ii++;
    }
  }
  coord n_transform_2d = {0.0, n_temp[0], n_temp[1]};

  if (f_xyz((coord){0.0,0.0,0.0}, n_transform_2d, alpha_2d, delta_b) > 0.0){
    //return -1;
    return zeta_rect((coord){0.0,0.0,0.0}, (coord){1.0,1.0,1.0}, n_transform_2d, alpha_2d, delta_b);
  }
  else{
    if (f_xyz((coord){1.0,1.0,1.0}, n_transform_2d, alpha_2d, delta_b) < 0.0){
      //return -1;
      return 0.0;
    }
    else{
      double zeta_sgs = 0.0;
      zeta_sgs += zeta_rect((coord){0.0,0.0,0.0}, (coord){1.0,1.0,1.0}, n_transform_2d, alpha_2d, delta_b)
              - zeta_tri((coord){0.0,0.0,0.0}, n_transform_2d, alpha_2d, delta_b);
      if (alpha_2d - n_transform_2d.z > 0.0){
        zeta_sgs += zeta_tri((coord){0.0,0.0,1.0}, n_transform_2d, alpha_2d, delta_b);
        }
      if (alpha_2d - n_transform_2d.y > 0.0){
        zeta_sgs += zeta_tri((coord){0.0,1.0,0.0}, n_transform_2d, alpha_2d, delta_b);
        }
      return zeta_sgs;
    }
  }
}

double f_term3(coord b, coord n, double alpha, double delta_b){
    return (sq(f_xyz(b,n,alpha,delta_b))+0.5)*erf(f_xyz(b,n,alpha,delta_b));}

double lambda_rect(coord a, coord b, coord n, double alpha, double delta_b){
    return sq(delta_b)/(2*n.y*n.z)*(
      f_term3((coord){0.0,b.y,b.z}, n, alpha, delta_b) - f_term3((coord){0.0,b.y,a.z},n, alpha, delta_b)
    - f_term3((coord){0.0,a.y,b.z}, n, alpha, delta_b) + f_term3((coord){0.0,a.y,a.z}, n, alpha, delta_b)
    + 1/sqrt(M_PI)*(f_xyz((coord){0.0,b.y,b.z},n,alpha,delta_b)*exp(-sq(f_xyz((coord){0.0,b.y,b.z},n,alpha,delta_b)))
      -f_xyz((coord){0.0,b.y,a.z},n,alpha,delta_b)*exp(-sq(f_xyz((coord){0.0,b.y,a.z},n,alpha,delta_b)))
      -f_xyz((coord){0.0,a.y,b.z},n,alpha,delta_b)*exp(-sq(f_xyz((coord){0.0,a.y,b.z},n,alpha,delta_b)))
      +f_xyz((coord){0.0,a.y,a.z},n,alpha,delta_b)*exp(-sq(f_xyz((coord){0.0,a.y,a.z},n,alpha,delta_b)))));}

double lambda_tri(coord b, coord n, double alpha, double delta_b){
    double az = (alpha-n.y*b.y)/n.z;
    return sq(delta_b)/(2*n.y*n.z)*(
      f_term3((coord){0.0,b.y,b.z}, n, alpha, delta_b) - f_term3((coord){0.0,b.y,az}, n, alpha, delta_b)
    + 1/sqrt(M_PI)*(f_xyz((coord){0.0,b.y,b.z},n,alpha,delta_b)*exp(-sq(f_xyz((coord){0.0,b.y,b.z},n,alpha,delta_b)))
      -f_xyz((coord){0.0,b.y,az},n,alpha,delta_b)*exp(-sq(f_xyz((coord){0.0,b.y,az},n,alpha,delta_b))))
    - 2*n.z*(b.z-az)/(sqrt(M_PI)*delta_b));}

double compute_lambda_sgs (coord n_transform, double alpha_transform, double delta_b, coord plane){ //plane should look like {-1,-1,a} where positive value components is plane of interest and value is z = a plane
  double n_temp[2] = {0.0,0.0};
  int ii = 0;
  double alpha_2d = 0.0;
  foreach_dimension(){
    if (plane.x != -1){
      alpha_2d = alpha_transform - n_transform.x*plane.x;
    }
    else{
      n_temp[ii] = n_transform.x;
      ii++;
    }
  }
  coord n_transform_2d = {0.0, n_temp[0], n_temp[1]};

  if (f_xyz((coord){0.0,0.0,0.0}, n_transform_2d, alpha_2d, delta_b) > 0.0){
    //return -1;
    return lambda_rect((coord){0.0,0.0,0.0}, (coord){1.0,1.0,1.0}, n_transform_2d, alpha_2d, delta_b);
  }
  else{
    if (f_xyz((coord){1.0,1.0,1.0}, n_transform_2d, alpha_2d, delta_b) < 0.0){
      //return -1;
      return 0.0;
    }
    else{
      double lambda_sgs = 0.0;
      lambda_sgs += lambda_rect((coord){0.0,0.0,0.0}, (coord){1.0,1.0,1.0}, n_transform_2d, alpha_2d, delta_b)
              - lambda_tri((coord){0.0,0.0,0.0}, n_transform_2d, alpha_2d, delta_b);
      if (alpha_2d - n_transform_2d.z > 0.0){
        lambda_sgs += lambda_tri((coord){0.0,0.0,1.0}, n_transform_2d, alpha_2d, delta_b);
        }
      if (alpha_2d - n_transform_2d.y > 0.0){
        lambda_sgs += lambda_tri((coord){0.0,1.0,0.0}, n_transform_2d, alpha_2d, delta_b);
        }
      return lambda_sgs;
    }
  }
}
