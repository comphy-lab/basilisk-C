/**
# Surface derivative calculation of surface tension coefficient Seric. et. al. 2018
*/
attribute {
    vector chi;
  }
scalar sigmaf[];
vector iforce[];
bool addmarangoni = 1;
#include "iforce.h"
#include "curvature.h"
event stability (i++)
{

  /**
  We first compute the minimum and maximum values of $\alpha/f_m =
  1/\rho$, as well as $\Delta_{min}$. */

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin))
    if (fm.x[] > 0.) {
      if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
      if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
      if (Delta < dmin) dmin = Delta;
    }
  double rhom = (1./amin + 1./amax)/2.;

  /**
  The maximum timestep is set using the sum of surface tension
  coefficients. */

  //~ double sigma = 0.;
  for (scalar c in interfaces)
    //~ sigma += c.sigma;
  
    foreach(serial){
      if (sigmaf[] > 0) {
        double dt = sqrt (rhom*cube(dmin)/(pi*sigmaf[]));
        if (dt < dtmax)
          dtmax = dt;
      }
      if (addmarangoni == 1){
      foreach_dimension()
       if (fabs(iforce.x[]) > 1.e-6) {
        rhom = fmin(rho1,rho2);
        double dt = 2.*dmin/(fabs(u.x[]) + sqrt(sq(u.x[]) + 4.*fabs(iforce.x[])/rhom ));  //sqrt (rhom*sq(dmin)/(iforce.x[]));
        if (dt < dtmax)
          dtmax = dt;
      }
      }
  
   }
}
static inline bool interf (double c)
{
  if (c > 1.e-6 && c < 1. - 1.e-6)
    return true;
  else
    return false;
}

#if dimension == 2
foreach_dimension()
static double marangoni_y (Point point, vector h, coord * dsigt)
{
   int ori = orientation(h.y[]);
  for (int i = -1; i <= 1; i++)
    if (h.y[i] == nodata || orientation(h.y[i]) != ori){
	      return nodata;
      }
  double hx = (h.y[1] - h.y[-1])/2.;
  double hxx = (h.y[1] + h.y[-1] - 2.*h.y[])/Delta;
  
  
   double sigmaleft = 0., fsumleft = 0., sigmaright = 0., fsumright = 0., sigmacenter = 0., fsumcenter = 0.;
      //I get mean value of sigma on 3 cells in a column (y-direction)
  for (int i = -1; i <= 1; i++){
        //  printf("\n%d\t%g", i, f[-1,i]); fflush(stdout);

//PALAS: CHECK HOW THIS WORKS IN X-ORIENTATION
    if (interf(f[-1,i])){
      sigmaleft += f[-1,i]*sigmaf[-1,i];
      fsumleft += f[-1,i];
    }

    if (interf(f[0,i])){
      sigmacenter += f[0,i]*sigmaf[0,i];
      fsumcenter += f[0,i];
    }
      
    if (interf(f[1,i])){
      sigmaright += f[1,i]*sigmaf[1,i];
      fsumright += f[1,i];
    }
   
  }

 

double ds = 0., dsigmadx = 0.;
  if (fsumleft >= 1.e-6 && fsumright >= 1.e-6){
    sigmaleft = sigmaleft/fsumleft; //mean of multiple cells in column
    sigmaright = sigmaright/fsumright; //mean of multiple cells in column
     ds = 2.*Delta*pow(1. + (hx*hx), 0.5); // length of h 
     dsigmadx = (sigmaright - sigmaleft)/ds;
  }
  else if (fsumleft < 1.e-6 && fsumright > 1.e-6 && fsumcenter > 1.e-6){
    sigmacenter = sigmacenter/fsumcenter; //mean of multiple cells in column
    sigmaright = sigmaright/fsumright; //mean of multiple cells in column
     ds = Delta*pow(1. + (hx*hx), 0.5); // length of h 
     dsigmadx = (sigmaright - sigmacenter)/ds;
  }
  else if (fsumleft > 1.e-6 && fsumright < 1.e-6 && fsumcenter > 1.e-6){
     sigmaleft = sigmaleft/fsumleft; //mean of multiple cells in column
    sigmacenter = sigmacenter/fsumcenter; //mean of multiple cells in column
     ds = Delta*pow(1. + (hx*hx), 0.5); // length of h 
     dsigmadx = (sigmacenter - sigmaleft)/ds;
  }
  else {
      dsigmadx = 0.;
      // Fix thses errors to find bugs
      // printf("\nError in calculating dsigmadx\n");
      // dump("error_file");
      // exit(0);

  }
   


  
    //printf("\n ds = %g",ds);
    //exit(0);
    coord signtangent;
    signtangent.x = 1;
    signtangent.y = 1;

  
 
  coord n;

  n.x = dsigt->x; n.y = dsigt->y; // get normal values passed using dsigt

    //1st quad normal, tangent dir is in increasing x direction (tripathi et. al.)
    if ( (n.x > 0. && n.y > 0.) ){
      //tangent is 4th quad
      signtangent.x = 1.; 
      signtangent.y = -1.; 
    }
      //2nd quad
    if ( (n.x < 0. && n.y > 0.) ){
      //tangent is in 1st quad
      signtangent.x = 1.; 
      signtangent.y = 1.; 
    }
      //3rd quad
    if ( (n.x < 0. && n.y < 0.) ){
      //tangent is in 4th quad
      signtangent.x = 1.; 
      signtangent.y = -1.; 
    }
      //4th quad
    if ( (n.x > 0. && n.y < 0.) ){
      //tangent is in 1st quad
      signtangent.x = 1.; 
      signtangent.y = 1.; 
    }


      dsigt->x = signtangent.x * dsigmadx / pow(1. + (hx*hx), 0.5 ); //unity component of tangent vector
      dsigt->y = signtangent.y * dsigmadx * fabs(hx)/pow(1. + (hx*hx), 0.5);
    
  
  
  return hxx/pow(1. + sq(hx), 3/2.);
}

#else // dimension == 3
foreach_dimension()
static double marangoni_z (Point point, vector h, coord * dsigt)
{
  bool hf_failed_for_curvature = 0;

  int ori = orientation(h.z[]);
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (h.z[i,j] == nodata || orientation(h.z[i,j]) != ori)
	      hf_failed_for_curvature = 1;

  double hx = (h.z[1] - h.z[-1])/2.;
  double hy = (h.z[0,1] - h.z[0,-1])/2.;

   double fl1,fl2,fl3;
   double fr1,fr2,fr3;
  double fc1, fc2, fc3;

   double sl1,sl0,slm1;
   double sr1,sr0,srm1;
  double hl = 0., hc = 0., hr = 0., hyl = 0., hyc =0. , hyr = 0.;


  fl1 = f[-1,0,-1]; fl2 = f[-1,0,0]; fl3 = f[-1,0,1];
  fc1 = f[0,0,-1]; fc2 = f[0,0,0]; fc3 = f[0,0,1];
  fr1 = f[1,0,-1]; fr2 = f[1,0,0]; fr3 = f[1,0,1];
  //We check whether the four heights on the x-y stencil are available to find the 
  //tangent of the interface in both directions
  if (h.z[1] == nodata || h.z[-1] == nodata || orientation(h.z[1]) != ori || orientation(h.z[-1]) != ori){
   
     for (int i = -2; i <= 2; i++){
      hl += f[-1,0,i];
      hr += f[1,0,i];
     }
     hx = (hr - hl)/2.;

  }
    
  if (h.z[0,-1] == nodata || h.z[0,1] == nodata || orientation(h.z[0,1]) != ori || orientation(h.z[0,-1]) != ori)
   {
     for (int i = -2; i <= 2; i++){
      hyl += f[0,-1,i];
      hyr += f[0,1,i];
     }
     hy = (hyr - hyl)/2.;
   }

  /**
  We "filter" the curvature using a weighted sum of the three
  second-derivatives in the $x$ and $y$ directions. This is necessary
  to avoid a numerical mode when the curvature is used to compute
  surface tension. */
  
  double filter = 0.2;
  double hxx = (filter*(h.z[1,1] + h.z[-1,1] - 2.*h.z[0,1]) +
		(h.z[1] + h.z[-1] - 2.*h.z[]) +
		filter*(h.z[1,-1] + h.z[-1,-1] - 2.*h.z[0,-1]))/
    ((1. + 2.*filter)*Delta);
  double hyy = (filter*(h.z[1,1] + h.z[1,-1] - 2.*h.z[1]) +
		(h.z[0,1] + h.z[0,-1] - 2.*h.z[]) +
		filter*(h.z[-1,1] + h.z[-1,-1] - 2.*h.z[-1]))/
    ((1. + 2.*filter)*Delta);
  double hxy = (h.z[1,1] + h.z[-1,-1] - h.z[1,-1] - h.z[-1,1])/(4.*Delta);

    
   double sigmaleft = 0., fsumleft = 0., sigmaright = 0., fsumright = 0., sigmacenter = 0., fsumcenter = 0.;
      //I get mean value of sigma on 3 cells in a column (y-direction)
  for (int i = -2; i <= 2; i++){
        //  printf("\n%d\t%g", i, f[-1,i]); fflush(stdout);

//Calculate derivative in X-Z plane
    if (interf(f[-1,0,i])){
      sigmaleft += f[-1,0,i]*sigmaf[-1,0,i];
      fsumleft += f[-1,0,i];
    }

    if (interf(f[0,0,i])){
      sigmacenter += f[0,0,i]*sigmaf[0,0,i];
      fsumcenter += f[0,0,i];
    }
      
    if (interf(f[1,0,i])){
      sigmaright += f[1,0,i]*sigmaf[1,0,i];
      fsumright += f[1,0,i];
    }
   
  }

  

//  sl1 = sigmaf[-1,-1]; sl0 = sigmaf[-1,0]; slm1 = sigmaf[-1,-1];
//   sr1 = sigmaf[1,-1]; sr0 = sigmaf[1,0]; srm1 = sigmaf[1,-1];

double dsx = 0., dsigmadx = 0.;
  if (fsumleft >= 1.e-6 && fsumright >= 1.e-6){

    sigmaleft = sigmaleft/fsumleft; //mean of multiple cells in column
    sigmaright = sigmaright/fsumright; //mean of multiple cells in column

     dsx = 2.*Delta*pow(1. + (hx*hx), 0.5); // length of h 
     dsigmadx = (sigmaright - sigmaleft)/dsx;
  }
  else if (fsumleft < 1.e-6 && fsumright > 1.e-6 && fsumcenter > 1.e-6){
    sigmacenter = sigmacenter/fsumcenter; //mean of multiple cells in column
    sigmaright = sigmaright/fsumright; //mean of multiple cells in column

     dsx = Delta*pow(1. + (hx*hx), 0.5); // length of h 
     dsigmadx = (sigmaright - sigmacenter)/dsx;
  }
  else if (fsumleft > 1.e-6 && fsumright < 1.e-6 && fsumcenter > 1.e-6){
     sigmaleft = sigmaleft/fsumleft; //mean of multiple cells in column
    sigmacenter = sigmacenter/fsumcenter; //mean of multiple cells in column

     dsx = Delta*pow(1. + (hx*hx), 0.5); // length of h 
     dsigmadx = (sigmacenter - sigmaleft)/dsx;
  }
  else {
    //   dsigmadx = 0.;
    // printf("\nError in calculating dsigmadx\n");
    // exit(0);

  }

sigmaleft = 0.; fsumleft = 0.; sigmaright = 0.; fsumright = 0.; sigmacenter = 0.; fsumcenter = 0.;
  //Calculate derivative in Y-Z plane
    for (int i = -2; i <= 2; i++){

      if (interf(f[0,-1,i])){
        sigmaleft += f[0,-1,i]*sigmaf[0,-1,i];
        fsumleft += f[0,-1,i];
      }

      if (interf(f[0,0,i])){
        sigmacenter += f[0,0,i]*sigmaf[0,0,i];
        fsumcenter += f[0,0,i];
      }
        
      if (interf(f[0,1,i])){
        sigmaright += f[0,1,i]*sigmaf[0,1,i];
        fsumright += f[0,1,i];
      }
   
    }

//    double fl1,fl0,flm1;
//    double fr1,fr0,frm1;

//    double sl1,sl0,slm1;
//    double sr1,sr0,srm1;
//   double hl,hr;
//   hl = h.y[-1];
//   hr = h.y[1];

//   fl1 = f[-1,-1]; fl0 = f[-1,0]; flm1 = f[-1,-1];
//   fr1 = f[1,-1]; fr0 = f[1,0]; frm1 = f[1,-1];

//  sl1 = sigmaf[-1,-1]; sl0 = sigmaf[-1,0]; slm1 = sigmaf[-1,-1];
//   sr1 = sigmaf[1,-1]; sr0 = sigmaf[1,0]; srm1 = sigmaf[1,-1];

double dsy = 0., dsigmady = 0.;
  if (fsumleft >= 1.e-6 && fsumright >= 1.e-6){
    sigmaleft = sigmaleft/fsumleft; //mean of multiple cells in column
    sigmaright = sigmaright/fsumright; //mean of multiple cells in column

     dsy = 2.*Delta*pow(1. + (hy*hy), 0.5); // length of h 
     dsigmady = (sigmaright - sigmaleft)/dsy;
  }
  else if (fsumleft < 1.e-6 && fsumright > 1.e-6 && fsumcenter > 1.e-6){
    sigmacenter = sigmacenter/fsumcenter; //mean of multiple cells in column
    sigmaright = sigmaright/fsumright; //mean of multiple cells in column

     dsy = Delta*pow(1. + (hy*hy), 0.5); // length of h 
     dsigmady = (sigmaright - sigmacenter)/dsy;
  }
  else if (fsumleft > 1.e-6 && fsumright < 1.e-6 && fsumcenter > 1.e-6){
     sigmaleft = sigmaleft/fsumleft; //mean of multiple cells in column
    sigmacenter = sigmacenter/fsumcenter; //mean of multiple cells in column

     dsy = Delta*pow(1. + (hy*hy), 0.5); // length of h 
     dsigmady = (sigmacenter - sigmaleft)/dsy;
  }
  else {
    //   dsigmady = 0.;
    // printf("\nError in calculating dsigmady\n");
    // exit(0);

  }
   
  

    //  if (fabs(fsumleft - f[-1,0]) < 1e-3){ // if there is only one cell in the column
    //   sigmaleft = sigmaf[-1,0]; //keep the ST coeff same as in the cell, discarding the mean value above
    // }
    //  if (fabs(fsumright - f[1,0]) < 1e-3){ // if there is only one cell in the column
    //   sigmaright = sigmaf[1,0]; //keep the ST coeff same as in the cell, discarding the mean value above
    // }
        //~ printf("\n%g\t%g", sigmaleft,sigmaright); fflush(stdout);

  
    //printf("\n ds = %g",ds);
    //exit(0);
    coord signtangent1;
    signtangent1.x = 1;
    signtangent1.y = 1;
    signtangent1.z = 1;

    coord signtangent2;
    signtangent2.x = 1;
    signtangent2.y = 1;
    signtangent2.z = 1;
 
  coord n;

  n.x = dsigt->x; n.y = dsigt->y; n.z = dsigt->z;  // get normal values stored using dsigt

    // normal xz which is in x-z plane 1st quad
    if ( (n.x > 0. && n.z > 0.)){
      //tangent is in 4th quad
      signtangent1.x = 1.; 
      signtangent1.z = -1.; 
    }
      // normal xz which is in x-z plane 2nd quad
    if ( (n.x < 0. && n.z > 0.)){
      //tangent is in 1st quad
      signtangent1.x = 1.; 
      signtangent1.z = 1.; 
    }
    // normal xz which is in x-z plane 3rd quad
    if ( (n.x < 0. && n.z < 0.) ){
      //tangent is in 4th quad
      signtangent1.x = 1.; 
      signtangent1.z = -1.; 
    }
      // normal xz which is in x-z plane 4th quad
    if ( (n.x < 0. && n.z > 0.)){
      // tangent is in 1st quad
      signtangent1.x = 1.; 
      signtangent1.z = 1.; 
    }

    // normal yz which is in y-z plane 1st quad
    if ( (n.y > 0. && n.z > 0.) ){
      //tangent is in 4th
      signtangent2.y = 1.; 
      signtangent2.z = -1.; 
    }
        // normal yz which is in y-z plane 2nd quad
    if ( (n.y < 0. && n.z > 0.) ){
      //1st
      signtangent2.y = 1.; 
      signtangent2.z = 1.; 
    }
       // normal yz which is in y-z plane 3rd quad
    if ( (n.y < 0. && n.z < 0.) ){
      //4th quad
      signtangent2.y = 1.; 
      signtangent2.z = -1.; 
    }
      // normal yz which is in y-z plane 4th quad
    if ( (n.y > 0. && n.z < 0.) ){
      signtangent2.y = 1.; 
      signtangent2.z = 1.; 
    }

  


      dsigt->x = signtangent1.x * dsigmadx / pow(1. + (hx*hx), 0.5 ); //unity component of tangent vector
      dsigt->y = signtangent2.y * dsigmady / pow(1. + (hy*hy), 0.5 ); 
      dsigt->z = signtangent1.z * dsigmadx * fabs(hx)/pow(1. + (hx*hx), 0.5) + signtangent2.z * dsigmady * fabs(hy)/pow(1. + (hy*hy), 0.5);

  
  //  dsigt->x = dsigmadx ;
  //   dsigt->y = dsigmadx ;
    // if ((dsigt->x) > 1.){
    //     double dsix = dsigt->x;
    //       exit(0);
    // }
    // if ((dsigt->y) > 10.){
    //     double dsiy = dsigt->y;
    //      exit(0);
    // }
   //~ dsigt->x = dsigmadx / sqrt(1 + sq(hx)); //unity component of tangent vector
    // dsigt->x = 0.;
    // dsigt->y = 0.;

 if (hf_failed_for_curvature == 0)
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
    else
    return nodata;
}
#endif

// foreach_dimension()
// static coord normal2_z (Point point, vector h)
// {
//   scalar hz = h.z;
//   if (hz[] == nodata)
//     return (coord){nodata, nodata, nodata};
//   int ori = orientation(hz[]);
//   double a = ori ? -1. : 1.;
//   coord n;
//   n.z = a;
//   foreach_dimension(2) {
//     if (allocated(-1) && hz[-1] != nodata && orientation(hz[-1]) == ori) {
//       if (allocated(1) && hz[1] != nodata && orientation(hz[1]) == ori)
// 	n.x = a*(hz[-1] - hz[1])/2.;
//       else
// 	n.x = a*(hz[-1] - hz[]);
//     }
//     else if (allocated(1) && hz[1] != nodata && orientation(hz[1]) == ori)
//       n.x = a*(hz[] - hz[1]);
//     else
//       n.x = nodata;
//   }
//   return n;
// }

// foreach_dimension()
// static coord normal_z (Point point, vector h) {
//   coord n = normal2_z (point, h);
//   double nn = fabs(n.x) + fabs(n.y) + fabs(n.z);
//   if (nn < nodata) {
//     foreach_dimension()
//       n.x /= nn;
//     return n;
//   }
//   return (coord){nodata, nodata, nodata};
// }

/**
We now need to choose one of the $x$, $y$ or $z$ height functions to
compute the curvature. This is done by the function below which
returns the HF curvature given a volume fraction field *c* and a
height function field *h*. */

static double height_marangoni (Point point, scalar c, vector h, coord * dsigt, coord * dsigt_hnc, int i)
{

  /**
  We first define pairs of normal coordinates *n* (computed by simple
  differencing of *c*) and corresponding HF curvature function *kappa*
  (defined above). */

  typedef struct {
    double n;
    double (* kappa) (Point, vector, coord *);
  } NormKappa;
  struct { NormKappa x, y, z; } n;  
  foreach_dimension()
    n.x.n = c[1] - c[-1], n.x.kappa = marangoni_x;
  double (* kappaf) (Point, vector, coord *) = NULL; NOT_UNUSED (kappaf);
  
  //Palas: Start working here set signs on dsigt
  //Palas: Get the relative sign of tangents and store it in dsigt
  foreach_dimension()
    dsigt->x = c[1] - c[-1];
  /**
  We sort these pairs in decreasing order of $|n|$. */
  
  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormKappa, n.x, n.y);
#if dimension == 3
  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormKappa, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormKappa, n.y, n.z);
#endif

  /**
  We try each curvature function in turn. */

    double kappa = nodata;
    int count_dimension = 1;
    foreach_dimension()
    {
      if (kappa == nodata)
      {
        kappa = n.x.kappa(point, h, dsigt); // pass dsigt with normal coordinates, it helps in calculating marangoni force
        // save the dsigmadx for the highest normal component
        if (count_dimension == 1)
        {
          dsigt_hnc->x = dsigt->x;
          dsigt_hnc->y = dsigt->y;
#if dimension == 3
          dsigt_hnc->z = dsigt->z;
#endif
        }
        if (kappa != nodata)
        {
          kappaf = n.x.kappa;
          if (n.x.n < 0.)
          {
            kappa = -kappa;
          }
        }
      }
      count_dimension++;
    }

    if (kappa != nodata)
    {

      /**
       We limit the maximum curvature to $1/\Delta$. */

      if (fabs(kappa) > 1. / Delta)
        kappa = sign(kappa) / Delta;

        /**
         We add the axisymmetric curvature if necessary. */
      
#if AXI
    double nr, r = y, hx;
    if (kappaf == marangoni_x) {
      hx = (height(h.x[0,1]) - height(h.x[0,-1]))/2.;
      nr = hx*(orientation(h.x[]) ? 1 : -1);
    }
    else {
      r += height(h.y[])*Delta;
      hx = (height(h.y[1,0]) - height(h.y[-1,0]))/2.;
      nr = orientation(h.y[]) ? -1 : 1;
    }
    /* limit the minimum radius to half the grid size */
    kappa += nr/max (sqrt(1. + sq(hx))*r, Delta/2.);
#endif
  }
  
  return kappa;
}

/**
## General curvature computation

We first need to define "interfacial cells" i.e. cells which contain
an interface. A simple test would just be that the volume fraction is
neither zero nor one. As usual things are more complicated because of
round-off errors. They can cause the interface to be exactly aligned
with cell boundaries, so that cells on either side of this interface
have fractions exactly equal to zero or one. The function below takes
this into account. */




/**
The function below computes the mean curvature *kappa* of the
interface defined by the volume fraction *c*. It uses a combination of
the methods above: statistics on the number of curvatures computed
which each method is returned in a *cstats* data structure. 

If *sigma* is different from zero the curvature is multiplied by *sigma*.

If *add* is *true*, the curvature (optionally multiplied by *sigma*)
is added to field *kappa*. */



struct Marangoni {
  scalar c, kappa;
	vector dsigtv;
  //double sigma; // I need variable sigma declared above
  bool add;
  int i;
};

trace
// cstats marangoni (struct Marangoni p)
cstats marangoni (scalar c, scalar kappa, vector dsigtv, bool add, int i )
{
  // scalar c = p.c, kappa = p.kappa;
	// vector dsigtv = p.dsigtv;
  //~ double sigma = p.sigma ? p.sigma : 1.;
  // scalar sigmafL = p.sigmafL ;
  int sh = 0, sf = 0, sa = 0, sc = 0;
  vector ch = c.height, h = automatic (ch);
  if (!ch.x.i)
    heights (c, h);
    // marangoni_heights (c, sgc);

  /**
  On trees we set the prolongation and restriction functions for
  the curvature. */
 //What needs to be done for variable sigma similar to kappa 
  #if TREE
    kappa.refine = kappa.prolongation = curvature_prolongation;
    kappa.restriction = curvature_restriction;

    dsigtv.x.refine = dsigtv.x.prolongation = curvature_prolongation;
    dsigtv.x.restriction = curvature_restriction;

    dsigtv.y.refine = dsigtv.y.prolongation = curvature_prolongation;
    dsigtv.y.restriction = curvature_restriction;
    #if dimension==3 
      dsigtv.z.refine = dsigtv.z.prolongation = curvature_prolongation;
      dsigtv.z.restriction = curvature_restriction;
    #endif
  #endif

  /**
  We first compute a temporary curvature *k*: a "clone" of
  $\kappa$. */
  
  scalar k[];
  scalar_clone (k, kappa);
  coord  dsigt;
  coord  dsigt_hnc; //for highest normal component if height function curvature is failing

  


  //~ fprintf (fp, "%g %g\n", t*11.1366559937, max);
  //~ fflush (fp);
  foreach(reduction(+:sh) reduction(+:sf)) {

    /**
    If we are not in an interfacial cell, we set $\kappa$ to *nodata*. */


    if (!interfacial (point, c)){
      k[] = nodata;
	    foreach_dimension(){
	    dsigtv.x[] = nodata;
      // iforce.x[] = 0.;
      } 
    }

    /**
    Otherwise we try the standard HF curvature calculation first, and
    the "mixed heights" HF curvature second. */ 
    
    else if ((k[] = height_marangoni (point, c, h, &dsigt, &dsigt_hnc, p.i)) != nodata){
	    foreach_dimension(){
	      dsigtv.x[] = dsigt.x;
        //  iforce.x[] = dsigt.x;
	    }
      // diff[] = k[];	
		  //~ fprintf( fp, "%g %g\n",x, dsigtv.x[]);
//~ exit(0);	    
      sh++;
    }
    else if ((k[] = height_curvature_fit (point, c, h)) != nodata){
      
      //if height_function fails
    foreach_dimension()
      dsigtv.x[] = dsigt_hnc.x; //use values from the dominant component of the normal
      // printf("\n %e %e %e", dsigt_hnc.x, dsigt_hnc.y, dsigt_hnc.z); fflush(stdout);
       sf++;

       
    }
     
  }
  
  //~ fclose(fp);
	//~ exit(0);
        //~ boundary({dsigtv});


foreach (reduction(+:sa) reduction(+:sc)) {
    
    /**
    We then construct the final curvature field using either the
    computed temporary curvature... */

    double kf;
    if (k[] < nodata)
      kf = k[];
    else if (interfacial (point, c)) {

      /**
      ...or the average of the curvatures in the $3^{d}$ neighborhood
      of interfacial cells. */
      
      double sk = 0., a = 0.;
      foreach_neighbor(1)
	if (k[] < nodata)
	  sk += k[], a++;
      if (a > 0.)
	kf = sk/a, sa++;
      else

	/**
	Empty neighborhood: we try centroids as a last resort. */

	kf = centroids_curvature_fit (point, c), sc++;
    //if height_function failed but get the dsigmadx
    double junk = height_marangoni (point, c, h, &dsigt, &dsigt_hnc, p.i);
    foreach_dimension()
      dsigtv.x[] = dsigt_hnc.x;
      // printf("\n %e %e %e", dsigt_hnc.x, dsigt_hnc.y, dsigt_hnc.z); fflush(stdout);
    }
    else
      kf = nodata;

    /**
    We add or set *kappa*. */
    
    if (kf == nodata)
      kappa[] = nodata;
    else if (add)
      kappa[] += sigmaf[]*kf;
    else
      kappa[] = sigmaf[]*kf;      
  }
	//~ printf("\n %d %d %d %d", sh, sf, sa, sc);
  //~ exit(0);
  return (cstats){sh, sf, sa, sc};
}

event acceleration (i++)
{
 


  /**
  We check for all VOF interfaces for which $\sigma$ is non-zero. */
 
  for (scalar f in interfaces){
    // if (f.sigma) {
      
      /**
      If $\phi$ is already allocated, we add $\sigma\kappa$, otherwise
      we allocate a new field and set it to $\sigma\kappa$. */

      scalar phi = f.phi;
      vector chi = f.chi;
      if (phi.i){
	marangoni (f, phi, chi, true, i);
      }
	else {
	phi = new scalar;
	chi = new vector;
	marangoni (f, phi, chi, false, i); // FIX the arguements
		
	f.phi = phi;
	f.chi = chi;
      }


  }

     /**
    We check for all VOF interfaces for which $\phi$ is allocated. The
    corresponding volume fraction fields will be stored in *list*. */

    scalar * list = NULL;
    for (scalar f in interfaces)
      if (f.phi.i) {
        list = list_add (list, f);

        /**
        To avoid undeterminations due to round-off errors, we remove
        values of the volume fraction larger than one or smaller than
        zero. */

        foreach()
    f[] = clamp (f[], 0., 1.);
      }

    /**
    On trees we need to make sure that the volume fraction gradient
    is computed exactly like the pressure gradient. This is necessary to
    ensure well-balancing of the pressure gradient and interfacial force
    term. To do so, we apply the same prolongation to the volume
    fraction field as applied to the pressure field. */
    
  #if TREE
    for (scalar f in list) {
      f.prolongation = p.prolongation;
      f.dirty = true; // boundary conditions need to be updated
    }
  #endif

 
    list = list_concat (NULL, {f});
    face vector ia = a;


      foreach_face(){

             vector chi = f.chi;
      
       double chif =
      (fabs(chi.x[]) < nodata && fabs(chi.x[-1]) < nodata) ?
      (chi.x[] + chi.x[-1])/2. :
      fabs(chi.x[]) < nodata ? chi.x[] :
      fabs(chi.x[-1]) < nodata ? chi.x[-1] :
      0.;

 

      double area = 0;
       if (interf(f[])) {
      coord n = interface_normal (point, f), p;      
      double alpha = plane_alpha (f[], n);
      // area of the bubble interface
      #if AXI
        area = y*pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);   
      #else
        area = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);   
      #endif

    }

       double deltam = area/dv() ; //+ ( f[0,1] - f[] )/Delta ;

      if (addmarangoni == 1){
        iforce.x[] = alpha.x[]/fm.x[]*chif*deltam; //save this quantity for time step restriction
        ia.x[] += iforce.x[]; //  *fabs(deltas);
        
      }

    }




    /**
    On trees, we need to restore the prolongation values for the
    volume fraction field. */
    
  #if TREE
    for (scalar f in list) {
      f.prolongation = fraction_refine;
      f.dirty = true; // boundary conditions need to be updated
    }
  #endif
    
for (scalar f in list) {
    vector chi = f.chi;
    delete ((scalar *) {chi});

  }
    free (list);


}
