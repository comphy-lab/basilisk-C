/**
## Operator split explicit diffusive scheme for computing interfacial diffusion on a VOF interface.

To develop an approximately 2nd order numerical interfacial diffusion scheme we use principles seen in embed.h
to compute flux at the interface and between interfacial/interfacial and interfacial/bulk cells*/

#define BGHOSTS 2

/*determines maximum time step for diffusion solver dtmax = stability_condition/FLUX_LIMIT.
Increasing this number will increase the overall transfer rate which is limited to mantain stability
might be beneficial to use a larger number in the case of high schmidt numbers on coarse meshes,
where the concentration gradient at the interface is much larger than the numerical gradient.
This will also affect the flux limiting for small cells (both interfacial transfer and bulk transfer).
Using a higher number will increase the number of subloops of the diffusion solver linearly.
*/
#define FLUX_LIMIT 100

//diffusion stability probably needs to be modified for axi symmetric case but not sure
double FTCS_diffusion_stability(face vector D, double dt_global){
  //double tmin = dt_global;
  double tmin = HUGE;
  foreach_face(reduction(min:tmin)){
    tmin = min(tmin,0.5*cm[]*sq(Delta)/(D.x[] + SEPS));
  }
  return tmin;
}


/*function computes aperture from the VOF reconstruction, we can compute the aperture from both sides by inputing 0 or -1 for the side argument
where 0 corresponds to the right side and -1 corresponds to the left side
*/

foreach_dimension()
  static inline double face_fraction_x(Point point, scalar f, vector nn, scalar alpha, int side, int i, int j, int k){

    //to ensure compatibility with the stencil if i == -2 side = 0
    if (i == -2)
      side = 0;

    if ((f[side+i,j,k] <= 0.0) || (f[side+i,j,k] >= 1.0)){
      return f[side+i,j,k];}

    else{

      /*
      in cases where the cell of interest is refined the normal and the alpha value
      stored in the stencil will be wrong. In this case we change the local stencil to have the correct normal and alpha
      and build the face fraction accordingly

      */
#if TREE
      if (is_refined(neighbor(side+i,j,k))){

        #if dimension == 2
        coord idx = {(double)(side + 2*i), (double)(2*j), (double)(2*k)};
        int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS;
        point.level++;
        point.i = _i + (int)(idx.x); point.j = _j + (int)(idx.y);
        POINT_VARIABLES;

        coord n_temp;
        double output = 0.0;

        for (int qq = 0; qq <=1; qq++){

          n_temp.x = 0.0;
          n_temp.y = nn.y[0,qq,0];
          n_temp.z = nn.z[0,qq,0];

          double alpha_temp= alpha[0,qq,0] + (0.5 + (float)side)*nn.x[0,qq,0];
          output += (((f[0,qq,0] >= 1.0) || (f[0,qq,0] <= 0.0)) ? f[0,qq,0] : rectangle_fraction(n_temp, alpha_temp, (coord){-0.5, -0.5,-0.5}, (coord){0.5,0.5,0.5}))/2;

        }

          point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
          point.level--;

        #else //dimension == 3

        coord idx = {(double)(side + 2*i), (double)(2*j), (double)(2*k)};
        int _i = 2*point.i - GHOSTS, _j = 2*point.j - GHOSTS, _k = 2*point.k - GHOSTS;
        point.level++;
        point.i = _i + (int)(idx.x); point.j = _j + (int)(idx.y); point.k = _k + (int)(idx.z);
        POINT_VARIABLES;

        coord n_temp;
        double output = 0.0;

        for (int qq = 0; qq <=1; qq++){

          for (int ww = 0; ww <=1; ww++){

            n_temp.x = 0.0;
            n_temp.y = nn.y[0,qq,ww];
            n_temp.z = nn.z[0,qq,ww];

            double alpha_temp = alpha[0,qq, ww] + (0.5 + (float)side)*nn.x[0,qq,ww];
            output += (((f[0,qq,ww] >= 1.0) || (f[0,qq,ww] <= 0.0)) ? f[0,qq,ww] : rectangle_fraction(n_temp, alpha_temp, (coord){-0.5, -0.5,-0.5}, (coord){0.5,0.5,0.5}))/4;
          }

        }

        point.i = (_i + GHOSTS)/2; point.j = (_j + GHOSTS)/2;
        point.level--;

        #endif
        return output;
    }
#endif

      double alpha_temp;
      coord n_temp;

      n_temp.x = 0.0;
      n_temp.y = nn.y[side+i,j,k];
      n_temp.z = nn.z[side+i,j,k];

      alpha_temp = alpha[side + i, j,k] + (0.5 + (float)side)*nn.x[side + i, j,k];

      return rectangle_fraction(n_temp, alpha_temp, (coord){-0.5, -0.5,-0.5}, (coord){0.5,0.5,0.5});
  }
}


/*approximation of the face condition seen in embed.h we use both one sided face fractions to provide a strict condition given
that a PLIC VOF reconstruction is discontinuous at cell faces
*/


//scalar cs[];
//face vector fs[];

#if dimension == 2
foreach_dimension()
bool face_condition_x(Point point, scalar f, vector n, scalar alpha, int i, int j, int k){
    return (face_fraction_x(point, f, n, alpha, 0, i, j, k) > 0.5 && \
            face_fraction_y(point, f, n, alpha, 0, j + (j < 0), i, k) && \
            face_fraction_y(point, f, n, alpha, 0, j +  (j < 0), i-1, k) && \
            face_fraction_x(point, f, n, alpha, -1, i, j, k) > 0.5 && \
            face_fraction_y(point, f, n, alpha, -1, j + (j < 0), i, k) && \
            face_fraction_y(point, f, n, alpha, -1, j +  (j < 0), i-1, k) && \
            f[i,j,k] && f[i-1,j,k]);
}
/*
#define face_condition(fs, cs)						\
  (fs.x[i,j] > 0.5 && fs.y[i,j + (j < 0)] && fs.y[i-1,j + (j < 0)] &&	\
   cs[i,j] && cs[i-1,j])
*/

foreach_dimension()
static inline double face_gradient_2_x(Point point, scalar c, scalar f, vector n, scalar alpha, int i, double face_area){
 int j = sign(face_fraction_x(point, f, n, alpha, 0, i, 1, 0) + face_fraction_x(point, f, n, alpha, -1, i, 1, 0) - face_fraction_x(point, f, n, alpha, 0, i, -1, 0) - face_fraction_x(point, f, n, alpha, -1, i, -1, 0));
 if (face_condition_x(point, f, n, alpha, i, j, 0))
   return ((1. + face_area)*(c[]/f[] - c[-1]/f[-1]) + (1. - face_area)*(c[0,j]/f[0,j] - c[-1,j]/f[-1,j]))/(2.*Delta);
 return (c[]/f[] - c[-1]/f[-1])/Delta;
}


#else //dimension == 3
foreach_dimension()
bool face_condition_x(Point point, scalar f, vector n, scalar alpha, int i, int j, int k){
  return (face_fraction_x(point, f, n, alpha, 0, i, j, k) > 0.5 && \
          (face_fraction_x(point, f, n, alpha, 0, i, j, 0) > 0.5 || \
          face_fraction_x(point, f, n, alpha, 0, i, 0, k) > 0.5) && \
          face_fraction_y(point, f, n, alpha, 0, j + (j < 0), 0, i) && \
          face_fraction_y(point, f, n, alpha, 0, j +  (j < 0), 0, i-1) && \
          face_fraction_y(point, f, n, alpha, 0, j + (j < 0), k, i) && \
          face_fraction_y(point, f, n, alpha, 0, j +  (j < 0), k, i-1) && \
          face_fraction_z(point, f, n, alpha, 0, k + (k < 0), i, 0) && \
          face_fraction_z(point, f, n, alpha, 0, k +  (k < 0), i-1, 0) && \
          face_fraction_z(point, f, n, alpha, 0, k + (k < 0), i, j) && \
          face_fraction_z(point, f, n, alpha, 0, k +  (k < 0), i-1, j) && \
          face_fraction_x(point, f, n, alpha, -1, i, j, k) > 0.5 && \
         (face_fraction_x(point, f, n, alpha, -1, i, j, 0) > 0.5 || \
          face_fraction_x(point, f, n, alpha, -1, i, 0, k) > 0.5) && \
          face_fraction_y(point, f, n, alpha, -1, j + (j < 0), 0, i) && \
          face_fraction_y(point, f, n, alpha, -1, j +  (j < 0), 0, i-1) && \
          face_fraction_y(point, f, n, alpha, -1, j + (j < 0), k, i) && \
          face_fraction_y(point, f, n, alpha, -1, j +  (j < 0), k, i-1) && \
          face_fraction_z(point, f, n, alpha, -1, k + (k < 0), i, 0) && \
          face_fraction_z(point, f, n, alpha, -1, k +  (k < 0), i-1, 0) && \
          face_fraction_z(point, f, n, alpha, -1, k + (k < 0), i, j) && \
          face_fraction_z(point, f, n, alpha, -1, k +  (k < 0), i-1, j) && \
          f[i-1,j,0] && f[i-1,0,k] && f[i-1,j,k] &&				\
          f[i,j,0] && f[i,0,k] && f[i,j,k]);
}
/*
#define face_condition(fs, cs)						\
  (fs.x[i,j,k] > 0.5 && (fs.x[i,j,0] > 0.5 || fs.x[i,0,k] > 0.5) &&	\
   fs.y[i,j + (j < 0),0] && fs.y[i-1,j + (j < 0),0] &&			\
   fs.y[i,j + (j < 0),k] && fs.y[i-1,j + (j < 0),k] &&			\
   fs.z[i,0,k + (k < 0)] && fs.z[i-1,0,k + (k < 0)] &&			\
   fs.z[i,j,k + (k < 0)] && fs.z[i-1,j,k + (k < 0)] &&			\
   cs[i-1,j,0] && cs[i-1,0,k] && cs[i-1,j,k] &&				\
   cs[i,j,0] && cs[i,0,k] && cs[i,j,k])
*/

/*function computes face barycentre from the VOF reconstruction, we compute the
barycentre by taking the average of the one sided barycentres. if one of the one sided
faces is cut and the other is empty we simply use barycentre of the cut side
*/

foreach_dimension()
static inline coord face_barycentre_x(Point point, scalar f, vector nn, scalar alpha, int i)
{
  coord p = {0.0,0.0,0.0}, p1;
  double num_used = 0.0;

  coord n_temp;

  for (int side = -1; side <= 0; side++){
    double face_area = face_fraction_x(point, f, nn, alpha, side, i, 0, 0);
    if (face_area >= 1.0){
      //p.x += 0.0; p.y += 0.0; p.z += 0.0;
      num_used ++;}
    else{
      if (face_area > 0.0){
        n_temp = (coord){nn.y[i+side,0,0], nn.z[i+side,0,0], 0};
        double nsum = 0.0;
        foreach_dimension()
          nsum += fabs(n_temp.x);

        if (nsum != 0.0){
        foreach_dimension()
          n_temp.x /= nsum;
        }

        double alpha = line_alpha (face_area, n_temp);
        line_center (n_temp, alpha, face_area, &p1);
        p.y += ((double *)&p1)[0], p.z += ((double *)&p1)[1], p.x += 0.;
        num_used++;
      }
    }
  }
  if (num_used != 0.0){
  foreach_dimension()
    p.x/=num_used;}
  return p;
}

/*
foreach_dimension()
static inline coord embed_face_barycentre_z(Point point, int i)
{
  // Young's normal calculation
  coord n1 = {0};
  double nn = 0.;
  scalar f = fs.z;
  foreach_dimension(2) {
    n1.x = (f[-1,-1,i] + 2.*f[-1,0,i] + f[-1,1,i] -
	    f[+1,-1,i] - 2.*f[+1,0,i] - f[+1,1,i]);
    nn += fabs(n1.x);
  }
  if (!nn)
    return (coord){0.,0.,0.};
  foreach_dimension(2)
    n1.x /= nn;
  // Position `p` of the face barycentre
  coord n, p1, p;
  ((double *)&n)[0] = n1.x, ((double *)&n)[1] = n1.y;
  double alpha = line_alpha (f[0,0,i], n);
  line_center (n, alpha, f[0,0,i], &p1);
  p.x = ((double *)&p1)[0], p.y = ((double *)&p1)[1], p.z = 0.;
  return p;
}*/

foreach_dimension()
static inline double face_gradient_2_x(Point point, scalar c, scalar f, vector n, scalar alpha, int i, double face_area){
  coord p = face_barycentre_x(point, f, n, alpha, i);
  // Bilinear interpolation of the gradient (see Fig. 1 of Schwartz et al., 2006)
  int j = sign(p.y), k = sign(p.z);
  if (face_condition_x(point, f, n, alpha, i, j, k)) {
    p.y = fabs(p.y), p.z = fabs(p.z);
    return (((c[i,0,0]/f[i,0,0] - c[i-1,0,0]/f[i-1,0,0])*(1. - p.y) +
	     (c[i,j,0]/f[i,j,0] - c[i-1,j,0]/f[i-1,j,0])*p.y)*(1. - p.z) +
	    ((c[i,0,k]/f[i,0,k] - c[i-1,0,k]/f[i-1,0,k])*(1. - p.y) +
	     (c[i,j,k]/f[i,j,k] - c[i-1,j,k]/f[i-1,j,k])*p.y)*p.z)/Delta;
  }
  return (c[i]/f[i] - c[i-1]/f[i-1])/Delta;
}

#endif

//one-dimensional diffusion finite volume scheme
foreach_dimension()
static void diffusion_sweep_x(scalar c, scalar f, face vector D, double dt_diffusion, vector n, scalar alpha){

  scalar cflux[];

  foreach_face(x){
    //test.x[] = 0.0;
    if ((f[] > 0.0) && (f[-1] > 0.0)){
      if ((f[] < 1.0) || (f[-1] < 1.0)){

        double face_area = (face_fraction_x(point, f, n, alpha, 0, 0, 0, 0) + face_fraction_x(point, f, n, alpha, -1, 0, 0, 0))/2;
        face_area = clamp(face_area, 0.0, 1.0);

        int i = 0;
        double grad = face_gradient_2_x(point, c, f, n, alpha, i, face_area);
        double grad_degenerate = (c[]/f[]-c[-1]/f[-1])/Delta;
        grad = (sign(grad) == sign(grad_degenerate)) ? grad : grad_degenerate;

        coord n_temp, p;
        double gamma_1 = 0.0, gamma_2 = 0.0;

        if (f[] < 1){
          n_temp.x = n.x[0,0,0];
          n_temp.y = n.y[0,0,0];
          n_temp.z = n.z[0,0,0];
          plane_center (n_temp, alpha[], f[], &p);
          gamma_1 = p.x;
        }

        if (f[-1] < 1){
          n_temp.x = n.x[-1,0,0];
          n_temp.y = n.y[-1,0,0];
          n_temp.z = n.z[-1,0,0];
          plane_center (n_temp, alpha[-1], f[-1], &p);
          gamma_2 = p.x;
        }
        double gamma = (gamma_1 - gamma_2 + 1);
        ///gamma = 1.0;

        cflux[] = grad*D.x[]*face_area*gamma;
        double limit = (c[]/f[]-c[-1]/f[-1])/dt_diffusion*f[]*f[-1]/(f[] + f[-1])*Delta;
        cflux[] = min(fabs(cflux[]), fabs(limit))*sign(cflux[]);}

      else{
        cflux[] = (c[]/f[]-c[-1]/f[-1])/Delta*D.x[];}
    }
    else{
      //test.x[] = 0.0;
      cflux[] = 0.0;}
  }


  #if TREE


    //IF NEEDS TO BE MODIFIED I NEED TO CHANGE ALL FACES AT ONCE NOT ONE AT A TIME

    //in faces that where the grid resolution changes, ensure stability for the coarse grid as well
    foreach_face(x){
      if ((is_prolongation(neighbor(0)) || is_prolongation(neighbor(-1))) && (coarse(f,0) > 0.0) && (coarse(f,-1) > 0.0)){

        double limit = (coarse(c,0)/coarse(f,0)-coarse(c,-1)/coarse(f,-1))/dt_diffusion*coarse(f,0)*coarse(f,-1)/(coarse(f,0) + coarse(f,-1))*Delta*2;

        #if dimension == 2

          coord temp_idx = {(float)(point.i + GHOSTS), (float)(point.j + GHOSTS), 0.0};
          int yy = -(((int)temp_idx.y)%2);

          double coarse_flux = 0.0;
          double wrong_sign_sum = 0.0;

          for (int ii = yy; ii <= yy+1; ii++){
                coarse_flux += cflux[0,ii]/2;
                wrong_sign_sum += (sign(cflux[0,ii])*sign(limit) == -1 ? cflux[0,ii]/2 : 0.0);}


          //if opposite sign scale fine fluxes such that coarse flux will be zero
          if (sign(coarse_flux)*sign(limit) == -1 && coarse_flux*limit !=  0.0){
            for (int ii = yy; ii <= yy+1; ii++){
              if (wrong_sign_sum == 0.0){
                cflux[0,ii] = 0.0;}
              else{
                cflux[0,ii] = (sign(cflux[0,ii]) != sign(limit) ? cflux[0,ii]*max(fabs((wrong_sign_sum - coarse_flux)/wrong_sign_sum) - 1e-16, 0.0) : cflux[0,ii]);
                }
          }
         }
         //if magnitude is too large scale fluxes such that coarse flux is equal to limit
         if ((sign(coarse_flux)*sign(limit) == 1 && fabs(coarse_flux) > fabs(limit)) || fabs(limit) == 0.0){
            for (int ii = yy; ii <= yy+1; ii++){
              if (coarse_flux-wrong_sign_sum == 0.0){
                cflux[0,ii] = 0.0;}
              else{
                cflux[0,ii] = (sign(cflux[0,ii]) == sign(limit) ? cflux[0,ii]*max(fabs((limit-wrong_sign_sum)/(coarse_flux-wrong_sign_sum))-1e-16, 0.0) : cflux[0,ii]);
              }
          }
        }

      #else //dimension == 3

        coord temp_idx = {(float)(point.i + GHOSTS), (float)(point.j + GHOSTS), (float)(point.k + GHOSTS)};
        int yy = -(((int)temp_idx.y)%2);
        int zz = -(((int)temp_idx.z)%2);

        double coarse_flux = 0.0;
        double wrong_sign_sum = 0.0;


        for (int ii = yy; ii <= yy+1; ii++)
          for (int jj = zz; jj <= zz+1; jj++){
              coarse_flux += cflux[0,ii,jj]/4;
              wrong_sign_sum += (sign(cflux[0,ii,jj]) != sign(limit) ? cflux[0,ii,jj]/4 : 0.0);}

        //if wrong sign scale fine fluxes such that coarse flux will be zero
        double temp_temp = 0.0;

        if (sign(coarse_flux)*sign(limit) == -1 && coarse_flux*limit != 0.0){
          for (int ii = yy; ii <= yy+1; ii++)
            for (int jj = zz; jj <= zz+1; jj++){
              if (wrong_sign_sum == 0.0){
                cflux[0,ii,jj] = 0.0;}
              else{
                cflux[0,ii,jj] = (sign(cflux[0,ii,jj]) != sign(limit) ? cflux[0,ii,jj]*max(fabs((wrong_sign_sum - coarse_flux)/wrong_sign_sum) - 1e-16, 0.0) : cflux[0,ii,jj]);

              temp_temp +=  cflux[0,ii,jj];}
          }
        }
        if((sign(coarse_flux)*sign(limit) == 1 && fabs(coarse_flux) > fabs(limit)) || fabs(limit) == 0.0){
          for (int ii = yy; ii <= yy+1; ii++)
            for (int jj = zz; jj <= zz+1; jj++){
              if (coarse_flux-wrong_sign_sum == 0.0){
                cflux[0,ii,jj] = 0.0;}
              else{
                cflux[0,ii,jj] = (sign(cflux[0,ii,jj]) == sign(limit) ? cflux[0,ii,jj]*max(fabs((limit-wrong_sign_sum)/(coarse_flux-wrong_sign_sum))-1e-16, 0.0) : cflux[0,ii,jj]);
              temp_temp +=  cflux[0,ii,jj]/4;}
          }
        }
    #endif
    }
  }
  #endif


#if TREE
  for (int l = depth() - 1; l >= 0; l--)
    foreach_halo (prolongation, l) {
#if dimension == 1
      if (is_refined (neighbor(-1)))
	       cflux[] = fine(cflux);
      if (is_refined (neighbor(1)))
	       cflux[1] = fine(cflux,2);
#elif dimension == 2
      if (is_refined (neighbor(-1))){
	       cflux[] = (fine(cflux,0,0) + fine(cflux,0,1))/2;}
      if (is_refined (neighbor(1))){
	       cflux[1] = (fine(cflux,2,0) + fine(cflux,2,1))/2;}
#else // dimension == 3
      if (is_refined (neighbor(-1)))
	       cflux[] = (fine(cflux,0,0,0) + fine(cflux,0,1,0) +
		       fine(cflux,0,0,1) + fine(cflux,0,1,1))/4.;
      if (is_refined (neighbor(1)))
	       cflux[1] = (fine(cflux,2,0,0) + fine(cflux,2,1,0) +
		       fine(cflux,2,0,1) + fine(cflux,2,1,1))/4.;
#endif
    }
#endif

  foreach() {
    if (f[] > 0.0){
      c[] += dt_diffusion*(-cflux[] + cflux[1])/Delta;
    }
    else{
      c[] = 0.0;}
  }

  boundary ({c});
}

/**
## Multi-dimensional diffusion */

void operator_split_diffusion(scalar c, scalar f, face vector D, int i, double dt_diffusion, vector n, scalar alpha)
{

  void (* diffusion_sweep[dimension]) (scalar, scalar, vector, double, vector, scalar);
  int d = 0;
  foreach_dimension()
    diffusion_sweep[d++] = diffusion_sweep_x;
  for (d = 0; d < dimension; d++)
    diffusion_sweep[(i + d) % dimension] (c, f, D, dt_diffusion, n, alpha);
}

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradient_x (Point point, scalar s, scalar f,
  					   vector n, scalar alpha, coord p, coord p2, double bc, double bc_inf){

    coord n_temp;
    foreach_dimension()
      n_temp.x = - n.x[];

    normalize (&n_temp);

    double d[2], v[2] = {nodata,nodata};
    bool defined = true;

    foreach_dimension()
      if (defined && !face_fraction_x(point, f, n, alpha, 0, (n_temp.x > 0.), 0, 0) && !face_fraction_x(point, f, n, alpha, -1, (n_temp.x > 0.), 0, 0))
        defined = false;

    if (defined)
      for (int l = 0; l <= 1; l++) {
        int i = (l + 1)*sign(n_temp.x);
        d[l] = (i - p.x)/n_temp.x;
        double y1 = p.y + d[l]*n_temp.y;
        int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
        y1 -= (double)j;

  #if dimension == 2
        if (face_fraction_x(point, f, n, alpha, 0, i + (i < 0), j, 0) && face_fraction_y(point, f, n, alpha, 0, j, i, 0) && face_fraction_y(point, f, n, alpha, 0, j+1, i, 0) &&
        face_fraction_x(point, f, n, alpha, -1, i + (i < 0), j, 0) && face_fraction_y(point, f, n, alpha, 0, j, i, 0) && face_fraction_y(point, f, n, alpha, -1, j+1, i, 0) &&
  	  f[i,j-1] && f[i,j] && f[i,j+1]){
  	v[l] = clamp(quadratic (y1, (s[i,j-1]/f[i,j-1]), (s[i,j]/f[i,j]), (s[i,j+1]/f[i,j+1])), min(bc, bc_inf), max(bc, bc_inf));
    }
  #else // dimension == 3
        double z = p.z + d[l]*n_temp.z;
        int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
        z -= (double)k;
        bool defined = face_fraction_x(point, f, n, alpha, 0, i + (i < 0), j, k) && face_fraction_x(point, f, n, alpha, -1, i + (i < 0), j, k);
        for (int m = -1; m <= 1 && defined; m++){

          if (!face_fraction_y(point, f, n, alpha, 0, j, k+m, i) ||
              !face_fraction_y(point, f, n, alpha, 0, j+1, k+m, i) ||
            	!face_fraction_z(point, f, n, alpha, 0, k, i, j+m) ||
              !face_fraction_z(point, f, n, alpha, 0, k+1, i, j+m) ||
              !face_fraction_y(point, f, n, alpha, -1, j, k+m, i) ||
              !face_fraction_y(point, f, n, alpha, -1, j+1, k+m, i) ||
              !face_fraction_z(point, f, n, alpha, -1, k, i, j+m) ||
              !face_fraction_z(point, f, n, alpha, -1, k+1, i, j+m) ||
            	!f[i,j+m,k-1] || !f[i,j+m,k] || !f[i,j+m,k+1])
                defined = false;
      }

        if (defined){
  	// bi-quadratic interpolation
        	v[l] =
        	  clamp(quadratic (z,
        		     clamp(quadratic (y1,
        				(s[i,j-1,k-1]/f[i,j-1,k-1]), (s[i,j,k-1]/f[i,j,k-1]), (s[i,j+1,k-1]/f[i,j+1,k-1])), min(bc, bc_inf), max(bc, bc_inf)),
        		     clamp(quadratic (y1,
        				(s[i,j-1,k]/f[i,j-1,k]),   (s[i,j,k]/f[i,j,k]),   (s[i,j+1,k]/f[i,j+1,k])), min(bc, bc_inf), max(bc, bc_inf)),
        		     clamp(quadratic (y1,
        				(s[i,j-1,k+1]/f[i,j-1,k+1]), (s[i,j,k+1]/f[i,j,k+1]), (s[i,j+1,k+1]/f[i,j+1,k+1])), min(bc, bc_inf), max(bc, bc_inf))),min(bc, bc_inf), max(bc, bc_inf));
          }
  #endif // dimension == 3
        else{
  	break;}
      }
    if (v[0] == nodata) {

      /**
      This is a degenerate case, we use the boundary value and the
      cell-center value to define the gradient. */

      d[0] = max(fabs(n.x[]*p2.x + n.y[]*p2.y + n.z[]*p2.z - alpha[])/sqrt(sq(n.x[])+sq(n.y[])+sq(n.z[])), 1e-16);
      return (bc - s[]/f[])/(d[0]*Delta);
    }

    /**
    For non-degenerate cases, the gradient is obtained using either
    second- or third-order estimates. We impose a condition that the gradient must be
    the same sign and greater in magnitude than the degenerate case. this helps mantain the
    stability of the algorithm and is in accord with the expected physical gradient at the interface*/

  /*  if (v[1] != nodata){ // third-order gradient

     //we place a constaint on v[1] to stabilize the gradient
      if (bc > bc_inf)
        v[1] = v[1] < v[0] ? v[1] : v[0];
      else
        v[1] = v[1] > v[0] ? v[1] : v[0];

      double grad_second_order = (bc - v[0])/(d[0]*Delta);
      double grad_degenerate = (bc - s[]/f[])/(max(fabs(n.x[]*p2.x + n.y[]*p2.y + n.z[]*p2.z - alpha[])/sqrt(sq(n.x[])+sq(n.y[])+sq(n.z[])), 1e-16)*Delta);
      //first we confirm that the second order gradient is the same sign as the first order degenerate case. This is required for the stability of the algorithm
      grad_second_order = (sign(grad_second_order) == sign(grad_degenerate)) ? grad_second_order : grad_degenerate;
      double output = (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
      // we place a constraint the that 3nd order gradient must be the same sign and greater than or equal to the magnitude of the second order gradient. This follows the physics of physio-absorption
      return (fabs(output) > fabs(grad_second_order) && sign(output) == sign(grad_second_order)) ? output : grad_second_order;}*/
    //second-order gradient
    double grad_degenerate = (bc - s[]/f[])/(max(fabs(n.x[]*p2.x + n.y[]*p2.y + n.z[]*p2.z - alpha[])/sqrt(sq(n.x[])+sq(n.y[])+sq(n.z[])), 1e-16)*Delta);
    double output = (bc - v[0])/(d[0]*Delta);
    return (sign(output) == sign(grad_degenerate)) ? output : grad_degenerate;
  }



double dirichlet_gradient2(Point point, scalar s, scalar f,
             coord n, double alpha, coord p2, double bc, double bc_inf){
  double d = max(fabs(n.x*p2.x + n.y*p2.y + n.z*p2.z - alpha)/sqrt(sq(n.x)+sq(n.y)+sq(n.z)), 1e-12);
  return (bc - s[]/f[])/(d*Delta);
}





double dirichlet_gradient (Point point, scalar s, scalar f,
			   vector n, scalar alpha, coord p, coord p2, double bc, double bc_inf)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x[]) >= fabs(n.y[]))
      return dirichlet_gradient_x (point, s, f, n, alpha, p, p2, bc, bc_inf);
#else // dimension == 3
  if (fabs(n.x[]) >= fabs(n.y[])) {
    if (fabs(n.x[]) >= fabs(n.z[]))
      return dirichlet_gradient_x (point, s, f, n, alpha, p, p2, bc, bc_inf);
  }
  else if (fabs(n.y[]) >= fabs(n.z[]))
    return dirichlet_gradient_y (point, s, f, n, alpha, p, p2, bc, bc_inf);
  return dirichlet_gradient_z (point, s, f, n, alpha, p, p2, bc, bc_inf);
#endif // dimension == 3
  return nodata;
}

//

void interfacial_transfer(scalar c, scalar f, face vector D, double dt_diffusion, double bc, vector n, scalar alpha, double bc_inf) {

  foreach()
    if ((f[] > 0) && (f[] < 1.0)){
      coord n_temp = {n.x[], n.y[], n.z[]}, p, p2;
      double area = plane_area_center(n_temp, alpha[], &p);
      #if dimension == 2
        line_center(n_temp, alpha[], f[], &p2);
      #else
        plane_center (n_temp, alpha[], f[], &p2);
      #endif
      double grad = dirichlet_gradient(point, c, f, n, alpha, p, p2, bc, bc_inf);
      double flux = D.x[]*area*grad*dt_diffusion/Delta;
      flux = min(fabs(flux), fabs(bc-c[]/f[])*f[])*sign(flux); //flux_limiter to mantain stability in small cells
      c[] += flux;
    }
    boundary ({c});
}

/*compute diffusion in subloop for global timestep where loop time step is based on stability condition of diffusion solver

  inside loop first compute species transfer across the interface and then compute diffusion in bulk

*/
void explicit_transfer_and_diffusion (scalar c, scalar f, face vector Diff, double dt_global, double bc, double bc_inf, int i) {

  vector n[];
  scalar alpha[];

  /*if height is declared use height functions to compute the normal. this is needed to have second order
  convergence of interfacial area. mixed young centered scheme produces 0th order convergence of surface
  area see http://basilisk.fr/sandbox/Antoonvh/interface.c
  */

  if (f.height.x.i){
    heights (f, f.height);
  }

  //modified vof reconstruction algorithm that computes normals with height function

  reconstruction2(f, n, alpha);

  //compute bulk diffusion timestep from stabilty condition and number of subloop iterations

  double dt_diffusion = FTCS_diffusion_stability(Diff, dt_global)/FLUX_LIMIT;

  double num_it = ceil(dt_global/dt_diffusion);
  dt_diffusion = dt_global/num_it;

  //double min_c;
  //double max_c;

  for (int ii = 0; ii<(int)num_it; ii++){
    interfacial_transfer(c, f, Diff, dt_diffusion, bc, n, alpha, bc_inf);
    /*min_c = HUGE;
    max_c = -HUGE;
    foreach(reduction(min:min_c), reduction(max:max_c))
      if (f[] > 0.0){
        min_c = min(min_c, c[]/f[]);
        max_c = max(max_c, c[]/f[]);
      }
    fprintf(stderr, "%g transfer %g %g\n", t, min_c, max_c);*/

    operator_split_diffusion(c,f,Diff,i+ii,dt_diffusion, n, alpha);

  /*  min_c = HUGE;
    max_c = -HUGE;
    foreach(reduction(min:min_c), reduction(max:max_c))
      if (f[] > 0.0){
        min_c = min(min_c, c[]/f[]);
        max_c = max(max_c, c[]/f[]);
      }
    fprintf(stderr, "%g bulk %g %g\n", t, min_c, max_c);*/
  }
}