/**
  We want to detect thin sheets and ligaments in 2 phase flows.
  First, we compute the quadratic form $f (x, x) = x_i x_j T_{ij}$,
  where $T_{ij}$ are the quadratic moments, by
  integrating over a shell of arbitrary thickness and radius.
  Finally, we compute the signature s of the quadratic form.  */

int find_moments_level(scalar f, double length, double target_delta){
  
  int level = depth();
  bool found = false;
  
  while (level >= 0 && !found){
    double num = 1<<level;
    if (length/num > target_delta) found = true;    
    level--;
  } 
    
  return level;
}

/**
  To compute the eigenvalues used in the signature method in 3D
  we use the GNU scientific library.
*/
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

void eigen ( double data[], double tau[], double dim){


  gsl_matrix_view m
    = gsl_matrix_view_array (data, dim, dim);

  gsl_vector_complex *eval = gsl_vector_complex_alloc (dim);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (dim, dim);

  gsl_eigen_nonsymmv_workspace * w =
    gsl_eigen_nonsymmv_alloc (dim);

  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec,
                           GSL_EIGEN_SORT_ABS_DESC);


    for (int i = 0; i < dim; i++)
      {
        gsl_complex eval_i
           = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i
           = gsl_matrix_complex_column (evec, i);

//         printf ("eigenvalue = %g + %gi\n", GSL_REAL(eval_i), GSL_IMAG(eval_i));
        tau[i]=GSL_REAL(eval_i);
      }

  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);

}

/**
This function is used to compute the signature at a given level of refinement `lev`. It can be found using the function `find_moments_level`.
*/
void compute_signature_neigh_level(scalar f, scalar phii, scalar s, int lev){
  
  scalar int_xx[], int_xy[], int_yy[];
#if dimension == 3 
  scalar int_xz[], int_yz[], int_zz[];
#endif  
  vector Tij[], sign[];
  
  foreach_level(lev){
    int_xx[] = 0;
    int_xy[] = 0;
    int_yy[] = 0;
#if dimension == 3    
    int_xz[] = 0;
    int_yz[] = 0;
    int_zz[] = 0;
#endif    
  }
  
  boundary ({phii});
  
  /** We integrate over the cells of the 5x5 stencil and subtract
   the values of those in the 3x3 stencil. 
   We can say the thickness of the shell is equal to the grid-size
   and the radius is twice the grid-size*/
  
  double mom_xx, mom_yy, mom_xy;
#if dimension == 3 
  double mom_xz, mom_yz, mom_zz;
#endif  
  
  foreach_level(lev){
    
    mom_xx = mom_xy = mom_yy = 0.;
#if dimension == 3 
    mom_xz = mom_yz = mom_zz = 0.;
#endif 
    
//     if (phii[] > -0.99){
      double xp = x;                               
      double yp = y;
#if dimension == 3 
      double zp = z;
#endif   
      
      foreach_neighbor() {
        mom_xx += sq(x - xp)*phii[];
        mom_yy += sq(y - yp)*phii[];
        mom_xy += (x - xp)*(y - yp)*phii[];
#if dimension == 3        
        mom_xz += (x - xp)*(z - zp)*phii[];
        mom_yz += (y - yp)*(z - zp)*phii[];
        mom_zz += sq(z - zp)*phii[];
#endif
      }
      
      foreach_neighbor(1) {
        mom_xx -= sq(x - xp)*phii[];
        mom_yy -= sq(y - yp)*phii[];
        mom_xy -= (x - xp)*(y - yp)*phii[];
#if dimension == 3        
        mom_xz -= (x - xp)*(z - zp)*phii[];
        mom_yz -= (y - yp)*(z - zp)*phii[];
        mom_zz -= sq(z - zp)*phii[];
#endif
      }
      
      int_xx[] = mom_xx;
      int_yy[] = mom_yy;
      int_xy[] = mom_xy;
#if dimension == 3  
      int_xz[] = mom_xz;
      int_yz[] = mom_yz;
      int_zz[] = mom_zz; 
#endif    
  }
  
  foreach_level(lev)
    foreach_dimension()
      Tij.x[] = 0;
    
    /** We now compute the eigenvalues of the quadratic form.
      For a 2x2 symmetric matrix the eigenvalues can be found analytically. 
      In the 3D case we call the function `eigen` that uses the GSL library.
      */
#if dimension == 2     
  double Txx, Tyy, Txy;
  double tau[2]={0.};
  
  foreach_level(lev){
    Txx = int_xx[]/sq(Delta);
    Tyy = int_yy[]/sq(Delta);
    Txy = int_xy[]/sq(Delta);
    double  data[] = { Txx, Txy, 
                       Txy, Tyy };
    
    eigen(data, tau, 2);
  
    Tij.x[] = tau[0];
    Tij.y[] = tau[1];

  }
#else            //dimension == 3  
  double Txx, Tyy, Tzz, Txy, Txz, Tyz=0.;
  double tau[3]={0.};
  
  foreach_level(lev){
    Txx = int_xx[]/sq(Delta);
    Tyy = int_yy[]/sq(Delta);
    Txy = int_xy[]/sq(Delta);
    Tzz = int_zz[]/sq(Delta); 
    Txz = int_xz[]/sq(Delta);
    Tyz = int_yz[]/sq(Delta);
    double  data[] = { Txx, Txy, Txz,
                       Txy, Tyy, Tyz,
                       Txz, Tyz, Tzz };
    eigen(data, tau, 3);
  
    Tij.x[] = tau[0];
    Tij.y[] = tau[1];
    Tij.z[] = tau[2];
  }
#endif
  
  /** The signature of the quadratic form are the signs of the eigenvalues.  */
  
  foreach_level(lev){
    foreach_dimension(){
        sign.x[] = fabs(Tij.x[])/(Tij.x[] + 1.e-10);  
        if (fabs(Tij.x[]) < 10.) sign.x[] = 0;
    }
  }
  
  /** We identify thin ligaments by looking at the signature. In 2D we can have the signatures: 
   
   * (+,+) bulk of phase
   * (+,-) thin ligament
   * (+,0) or (-,0) interface  
   * (-,-) outside of phase          */
#if dimension ==2  
  foreach_level(lev){
    if (sign.x[]*sign.y[] > 0.999 && sign.x[] > 0.999) s[] = 2;   //bulk of phase
    if (sign.x[]*sign.y[] > 0.999 && sign.x[] < -0.999) s[] = -1; //outside of phase
    if (fabs(sign.x[]*sign.y[]) < 0.01) s[] = 0;                 //interface
    if (sign.x[]*sign.y[] < -0.999) s[] = 1;                     //ligament
  }
#else //dimension == 3
  foreach_level(lev){
    s[] = 1; //thin
    if (sign.x[] == 0. || sign.y[] == 0. || sign.z[] == 0.) s[] = 0;   //interface
    if (sign.x[] > 0.999 && sign.y[] > 0.999 && sign.z[] > 0.999) s[] = 2;   //bulk of phase
    if (sign.x[] < -0.999 && sign.y[] < -0.999 && sign.z[] < -0.999) s[] = -1; //outside of phase
                  
  }
#endif
}

/**
 This function is used to perforate thin sheets. The scalar `s` must be filled using the `compute_signature_neigh_level` function.
*/
void change_topology (scalar f, scalar s, scalar M, int lev, const int max_change, bool large){
  
  double f_avg[max_change];  // average f in the neighbor
  int num = 0.; int i_change = 0;
  srand(time(NULL));   // Initialization of random seed.
  int r;      
  
  Cache to_perf = {0};
#ifndef SQUARES
  int i;
  scalar hole[];
  double xc[max_change], yc[max_change], zc[max_change];
#endif  
  foreach_level(lev){  
    M[] = 0.;
    r = rand() % 50;
    
    if (s[] > 0.99 && s[] < 1.01 && r==0 && i_change < max_change){
      f_avg[i_change] = 0.;
      num = 0;
#ifndef SQUARES      
      xc[i_change] = x; yc[i_change] = y; zc[i_change] = z;
#endif
#ifdef SQUARES      
      if (large == false){
        foreach_neighbor(1){       //compute average f 
          f_avg[i_change] += f[];
          cache_append (&to_perf, point, i_change);
          num += 1;
        }
      }
      else{
        foreach_neighbor(){       //compute average f 
          f_avg[i_change] += f[];
          cache_append (&to_perf, point, i_change);
          num += 1;
        }
      }
      
      f_avg[i_change]/=num;
      cache_shrink (&to_perf);
#endif            
      i_change++;
    }
  }
  
#ifndef SQUARES
  printf("Holes are spheres \n");
  vertex scalar phi[];
  double Ddelta = 0.001/64.;
  double R = 1.2*2.5*Ddelta;
  foreach_vertex() {
    phi[] = HUGE;
    for (i = 0; i < max_change; i++){
      phi[] = intersection (phi[], (sq(x - xc[i])/sq(1)+sq(y -yc[i])/sq(1)+sq(z -zc[i])/sq(1) - sq(R))); 
//       phi[] = intersection (phi[], (/*sq(x - xc[i]) +*/ (sq(y - yc[i]) + sq(z - zc[i]) - sq(R)))); 
    }
  }
  boundary ({phi});
  fractions (phi, hole);
  foreach(){
    if (f[] > 1e-4 && hole[] < 1 - 1.e-4){
      f[] = hole[];
    }
  }
  boundary({f});
#endif
  
#ifdef SQUARES  
  printf("Holes are squares \n");
  foreach_cache(to_perf){
    s[] = -2.;
    M[] = f[] - (1 - f_avg[_flags]); 
    if (f_avg[_flags] < 0.3) {
      s[] = -3.;
      f[] = 1.;
    }
    else
      f[] = 0.;    
  }
#endif  
  
  for (int ilev = lev; ilev < depth(); ilev++)  
    foreach_level(ilev){
      s.prolongation = M.prolongation = f.prolongation = refine_injection;
      if(is_refined(cell) && s[] == -2){
        s.prolongation (point, s);
        M.prolongation (point, M);
        f.prolongation (point, f);
      }
    }
  free(to_perf.p);
  boundary(all);
}



