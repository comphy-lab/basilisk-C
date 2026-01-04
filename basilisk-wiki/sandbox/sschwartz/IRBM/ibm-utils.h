/**
This file contains multiple helper/utility functions for the immersed boundary method. */

extern scalar ibm;
extern face vector ibmf;
extern coord vc;

/**
Default spreading function is Peskin's 4-point one, named phi_func below. */

int spreading_func = 4;

#define distance(a,b) (sqrt(sq(a) + sq(b)))
#define distance3D(a,b,c) (sqrt(sq(a) + sq(b) + sq(c)))

void bilinear_interpolation (Point point, vector uv, coord pc, coord * uc)
{
    coord uci;
    double xc = pc.x;
    double yc = pc.y;

    double xx = x - xc;
    double yy = y - yc;

    double x2 = x - sign(xx)*Delta;
    double y2 = y - sign(yy)*Delta;

    foreach_dimension() {
        uci.x = (x2 - xc)*(y2 - yc) * uv.x[] +
                (xc - x) * (y2 - yc) * uv.x[-sign(xx)] +
        	    (x2 - xc) * (yc - y) * uv.x[0,-sign(yy)] +
	            (xc - x) * (yc - y) * uv.x[-sign(xx),-sign(yy)];
        uci.x /=  ((x2 - x)*(y2 - y));
    }

    *uc = uci;
}

double scalar_bilinear_interpolation (Point point, scalar p, coord pc) 
{
    double pci;
    double xc = pc.x;
    double yc = pc.y;

    double xx = x - xc;
    double yy = y - yc;

    double x2 = x - sign(xx)*Delta;
    double y2 = y - sign(yy)*Delta;

    pci = (x2 - xc)*(y2 - yc) * p[] +
          (xc - x) * (y2 - yc) * p[-sign(xx)] +
          (x2 - xc) * (yc - y) * p[0,-sign(yy)] +
	      (xc - x) * (yc - y) * p[-sign(xx),-sign(yy)];
    pci /=  ((x2 - x)*(y2 - y));

  return pci;
}

/**
The following few functions are different types of spreading functions for smearing
the body force to nearby cells. The one used in the paper is phi_func. Note the number
in the function's name denotes the diameter (in cells) of the circle/sphere of 
influence. */

/**
Basic 4-Point spreading funciton from Peskin.*/
double phi_func (double x, double h) 
{
    double r = x / h;
    double phi;

    if (fabs(r) <= 1)
      phi = (1./8.) * (3 - (2 * fabs(r)) + sqrt (1 + (4 * fabs(r)) - (4 * sq(r))));
    else if (fabs(r) > 1 && fabs(r) <= 2)
      phi = (1./8.) * (5 - (2 * fabs(r)) - sqrt (-7 + (12 * fabs(r)) - ( 4 * sq(r))));
    else
      phi = 0;
    return phi;
}

double phi_five_smooth (double x, double h)
{
    double r = x / h;
    double phi;

    if (fabs(r) <= 0.5)
        phi = 3./8. + M_PI/32. - sq(r)/4.;
    else if (fabs(r) > 0.5 && fabs(r) <= 1.5)
        phi = 1./4. + (1. - fabs(r)) / (8.) * sqrt (-2. + 8. * fabs(r) - 4*sq(r)) - 
              1./8. * asin(sqrt(2.) * (fabs(r) - 1.));
    else if (fabs(r) > 1.5 && fabs(r) <= 2.5)
        phi = 17./16. - M_PI/64. - (3*fabs(r))/4. + sq(r)/8. + (fabs(r) - 2)/16. * 
              sqrt (-14. + 16*fabs(r) - 4*sq(r)) + 1./16. * asin(sqrt(2.) * (fabs(r) - 2));
    else
        phi = 0.;

    return phi;
}

double phi_three (double x, double h)
{
    double r = x / h;
    double phi;

    if (fabs(r) <= 0.5)
        phi = 1./3. * (1 + sqrt(1 - 3 * sq(r)));
    else if (fabs(r) <= 1.5 && fabs(r) > 0.5)
        phi = 1./6. * (5 - 3 * fabs(r) - sqrt (1 - 3 * sq(1 - fabs(r))));
    else
        phi = 0.;

    return phi;
}

double delta_func (coord sPoint, coord markerPoint, double dv, double Delta)
{
    double dphi = 1.;
    foreach_dimension(){
        if (spreading_func == 4)
            dphi = dphi * phi_func(sPoint.x - markerPoint.x,Delta);
        else if (spreading_func == 5)
    	    dphi = dphi * phi_five_smooth(sPoint.x - markerPoint.x,Delta);
        else if (spreading_func == 3)
        	dphi = dphi * phi_three(sPoint.x - markerPoint.x,Delta);
        else
            dphi = dphi * phi_func(sPoint.x - markerPoint.x,Delta);
    }
    return dphi/dv; 
}

bool empty_neighbor (Point point, coord * pc, scalar ibm)
{
    coord pc_temp;
    double temp_ibm = ibm[];
    double xc = x, yc = y;
    double max_d = 1e6;
    bool neighbor = false;

    foreach_neighbor(1)
        if (ibm[] == 1 && temp_ibm == 0 && (distance(x - xc, y - yc) < max_d)) {
            pc_temp.x = (xc + x) / 2.;
            pc_temp.y = (yc + y) / 2.;
            max_d = distance(x - xc, y - yc);
            neighbor = true;
            *pc = pc_temp;
        }
   return neighbor;
}


/**
Returns the marker point (the centroid of the interface fragment) for a cell
on the solid interface. Note the point is returned in global coordinates. */

double marker_point (Point point, scalar ibm, face vector ibmf, coord* markerPoint, coord* nu = NULL)
{
    coord cellCenter = {x, y, z};
    coord n = facet_normal (point, ibm, ibmf);

    double alpha = plane_alpha (ibm[], n);
    double area = plane_area_center (n, alpha, markerPoint);
    foreach_dimension()
        n.x *= -1;
    normalize(&n);
    if (nu) *nu = n;

    foreach_dimension()
        markerPoint->x = cellCenter.x + markerPoint->x*Delta;
    return area;
}

/**
The next few functions are for various types of extrapolation using the approach
outlined in the paper. This could probably be cleaned up in a more concise way. */

coord extrapolate_gradient (Point point, scalar s, coord markerCoord, coord n, vector v)
{
    double weight[5][5] = {0};
    double weightSum = 0.;
    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            if (s[i,j] == 0) {

                coord cellCenter = {x + Delta*i,y + Delta*j}, d; // how to get x and y from cell centers?
                foreach_dimension()
                    d.x = markerCoord.x - cellCenter.x;

                double distanceMag = distance (d.x, d.y);
                double normalProjection = (n.x * d.x) + (n.y * d.y);

                weight[i][j] = sq(distanceMag) * fabs(normalProjection);

                weightSum += weight[i][j];
            }
            else
                weight[i][j] = 0.;
        }
    }

    coord dudn = {0};

    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            foreach_dimension()
                dudn.x += (weight[i][j]/weightSum) * v.x[i,j];
        }
    }

    return dudn;
}

double extrapolate_scalar (Point point, scalar s, coord interpolatePoint, coord n, scalar p)
{
#if dimension == 2
    double weight[5][5] = {0};
    double weightSum = 0.;
    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            if (s[i,j] > 0.5) {

                coord cellCenter = {x + Delta*i,y + Delta*j}, d;
                foreach_dimension()
                    d.x = interpolatePoint.x - cellCenter.x;

                double distanceMag = distance (d.x, d.y);
                double normalProjection = (n.x * d.x) + (n.y * d.y);

                weight[i][j] = sq(distanceMag) * fabs(normalProjection);

                weightSum += weight[i][j];
            }
            else
                weight[i][j] = 0.;
        }
    }

    double interpolatedScalar = 0;

    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            interpolatedScalar += (weight[i][j]/(weightSum + SEPS)) * p[i,j];
        }
    }
#else
    double weight[5][5][5] = {0};
    double weightSum = 0.;
    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            for (int k = -2; k <= 2; k++) {
                if (s[i,j,k] > 0.5) {

                    coord cellCenter = {x + Delta*i, y + Delta*j, z + Delta*k}, d;
                    foreach_dimension()
                        d.x = interpolatePoint.x - cellCenter.x;

                    double distanceMag = distance3D (d.x, d.y, d.z);
                    double normalProjection = (n.x * d.x) + (n.y * d.y) + (n.z * d.z);

                    weight[i+2][j+2][k+2] = sq(distanceMag) * fabs(normalProjection);

                    weightSum += weight[i+2][j+2][k+2];
                }
                else
                   weight[i+2][j+2][k+2] = 0.;
            }
        }
    }

    double interpolatedScalar = 0;

    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            for (int k = -2; k <= 2; k++) {
                interpolatedScalar += (weight[i+2][j+2][k+2]/(weightSum + SEPS)) * p[i,j,k];
            }
        }
    }
#endif
    return interpolatedScalar;
}

double extrapolate_scalarv2 (Point point, scalar s, coord markerCoord, coord n, scalar p)
{
    double weight[5][5] = {0};
    double weightSum = 0.;
    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            if (s[i,j] == 0) {

                coord cellCenter = {x + Delta*i,y + Delta*j}, d;
                foreach_dimension()
                    d.x = markerCoord.x - cellCenter.x;

                double distanceMag = distance (d.x, d.y);
                double normalProjection = (n.x * d.x) + (n.y * d.y);

                weight[i][j] = sq(distanceMag) * fabs(normalProjection);

                weightSum += weight[i][j];
            }
            else
                weight[i][j] = 0.;
        }
    }

    double pressure = 0;

    for (int i = -2; i <= 2; i++) {
        for (int j = -2; j <= 2; j++) {
            pressure += (weight[i][j]/weightSum) * p[i,j];
        }
    }

    return pressure;
}

/**
One way to measure the accuracy of the method is calculating how much flux there is
across the solid interface. Another way is to get the interfacial velocity via 
interpolation, which is the next function. */

typedef struct SurfaceFlux
{
    double favg, fmin, fmax;   // flux variables
    double navg, nmin, nmax;   // normal variables
    double tavg, tmin, tmax;   // tangent variables
} SurfaceFlux;


scalar unf[], utf[], uflux[];

SurfaceFlux get_surfaceflux (vector u)
{
    double fmin = 1e30, fmax = -1e30, favg = 0;
    double nmin = 1e30, nmax = -1e30, navg = 0;
    double tmin = 1e30, tmax = -1e30, tavg = 0;
    int nump = 0;
    foreach(reduction(min:fmin) reduction(max:fmax) reduction(+:favg)
            reduction(min:nmin) reduction(max:nmax) reduction(+:navg)
            reduction(min:tmin) reduction(max:tmax) reduction(+:tavg) reduction(+:nump)) {

        if (ibm[] > 0 && ibm[] < 1) {
            coord markerPoint, n;
            double area = marker_point (point, ibm, ibmf, &markerPoint, &n);
            double uint = interpolate_linear (point, u.x, markerPoint.x, markerPoint.y, 
                                                                         markerPoint.z);
            double vint = interpolate_linear (point, u.y, markerPoint.x, markerPoint.y, 
                                                                         markerPoint.z);
            // project velocity to interface normal direction
            double un = uint*n.x + vint*n.y;
            double flux = un * area * Delta;

            navg += un, favg += flux;

            if (flux < fmin) fmin = flux;
            if (flux > fmax) fmax = flux;
            if (un < nmin) nmin = un;
            if (un > nmax) nmax = un;
            unf[] = un, uflux[] = flux;

            // project velocity to interface tangent direction
            coord tn = {-n.y, n.x};
            double ut = uint*tn.x + vint*tn.y;
            tavg += ut;
            if (ut < tmin) tmin = ut;
            if (ut > tmax) tmax = ut;
            utf[] = ut;

            nump++;
        }
        else
            unf[] = 0, utf[] = 0, uflux[] = 0;
    }
    SurfaceFlux tempf = {favg/nump, fmin, fmax, navg/nump, nmin, nmax, tavg/nump, tmin, tmax};
    return tempf;
}

/**
We want to get statistic on solid interfacial/surface points to see how well the 
no-slip boundary condition is being enforced. This is done by interpolating the
velocity filed on the solid interface and then using it for various statisitcs,
e.g., min, max, L2. */

typedef struct SurfacePoints
{
    // basic velocity statistics
    double uavg, umin, umax;
    double vavg, vmin, vmax;
    double wavg, wmin, wmax;

    // more advanced velocity statistics
    double ul2, vl2, wl2;

    int total; // # of marker points

}SurfacePoints;


scalar uerror[], verror[];
#if dimension == 3
scalar werror[];
#endif

SurfacePoints get_surfacepoints (vector u)
{
    double uavg = 0, umin = 1e30, umax = -1e30;
    double vavg = 0, vmin = 1e30, vmax = -1e30;
    double wavg = 0, wmin = 1e30, wmax = -1e30;
    double ul2 = 0, vl2 = 0, wl2 = 0;

    int nump = 0;

    foreach(reduction(min:umin) reduction(max:umax) reduction(+:uavg) reduction(+:ul2)
            reduction(min:vmin) reduction(max:vmax) reduction(+:vavg) reduction(+:vl2)
            reduction(min:wmin) reduction(max:wmax) reduction(+:wavg) reduction(+:wl2)
            reduction(+:nump)) {
        if (ibm[] > 0 && ibm[] < 1) {
            coord markerPoint, n;
            double area = marker_point (point, ibm, ibmf, &markerPoint, &n);
            double uint = interpolate_linear (point, u.x, markerPoint.x, markerPoint.y, markerPoint.z);
            double vint = interpolate_linear (point, u.y, markerPoint.x, markerPoint.y, markerPoint.z);

            if (uint > umax) umax = uint;
            if (uint < umin) umin = uint;
            if (vint > vmax) vmax = vint;
            if (vint < vmin) vmin = vint;

            uavg += uint, vavg += vint;
            ul2 += sq(uint - vc.x), vl2 += sq(vint - vc.y);
            uerror[] = uint - vc.x;
            verror[] = vint - vc.y;

#if dimension == 3
            double wint = interpolate_linear (point, u.z, markerPoint.x, markerPoint.y, markerPoint.z);
            if (wint > wmax) wmax = wint;
            if (wint < wmin) wmin = wint;

            wavg += wint, wl2 += sq(wint - vc.z);
            werror[] = wint - vc.z;
#endif
            nump++;
        }
    }

    SurfacePoints temp = {uavg/nump, umin, umax, vavg/nump, vmin, vmax, 
                          wavg/nump, wmin, wmax, ul2, vl2, wl2, nump};

    return temp;
}




/**
The next few functions are taken from embed to calculate interfacial force, although
they aren't used. Instead, we rely on ibm_force to get lift and drag. */

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradient_x (Point point, scalar s, scalar ibm,
					   coord n, coord p, double bc, double * coef)
{
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !ibmf.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (ibmf.x[i + (i < 0),j] && ibmf.y[i,j] && ibmf.y[i,j+1] &&
	  ibm[i,j-1] && ibm[i,j] && ibm[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = ibmf.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!ibmf.y[i,j,k+m] || !ibmf.y[i,j+1,k+m] ||
	    !ibmf.z[i,j+m,k] || !ibmf.z[i,j+m,k+1] ||
	    !ibm[i,j+m,k-1] || !ibm[i,j+m,k] || !ibm[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
  if (v[0] == nodata) {

    /**
    This is a degenerate case, we use the boundary value and the
    cell-center value to define the gradient. */
	
    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
  *coef = 0.;
  if (v[1] != nodata) // third-order gradient
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double dirichlet_gradient (Point point, scalar s, scalar ibm,
			   coord n, coord p, double bc, double * coef)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient_x (point, s, ibm, n, p, bc, coef);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient_x (point, s, ibm, n, p, bc, coef);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient_y (point, s, ibm, n, p, bc, coef);
  return dirichlet_gradient_z (point, s, ibm, n, p, bc, coef);
#endif // dimension == 3
  return nodata;
}

static inline
coord ibm_gradient (Point point, vector u, coord p, coord n)
{
    coord cellCenter = {x,y,z}, midPoint, dudn;
    foreach_dimension() {
        midPoint.x = cellCenter.x + p.x*Delta;
     }
#if dimension == 2
    midPoint.z = 0;
#endif
    double px = midPoint.x, py = midPoint.y, pz = midPoint.z;

    // calculate velocity at interface, may not be exactly equal to vc
    double uint_x = interpolate_linear (point, u.x, px, py, pz);
    double uint_y = interpolate_linear (point, u.y, px, py, pz);
#if dimension == 3
    double uint_z = interpolate_linear (point, u.z, px, py, pz);
#endif
    foreach_dimension() {
        double vb = uint_x;
        double val;
        dudn.x = dirichlet_gradient (point, u.x, ibm, n, p, vb, &val);
        if (dudn.x == nodata)
          dudn.x = 0.;
    }
    return dudn;
}

double ibm_vorticity (Point point, vector u, coord p, coord n)
{
    coord dudn = ibm_gradient (point, u, p, n);

    return -(dudn.y*n.x - dudn.x*n.y);
}


/**
ibm_force_int calculates the pressure force, Fp, and viscous force, Fmu, acting on the
immersed boundary by interpolation. This is different from ibm_force, which only sums
the body force to get the aerodynamic forces.
*/

void ibm_force_int (scalar p, vector u, face vector mu, coord * Fp, coord * Fmu)
{
    coord Fps = {0}, Fmus = {0};
    foreach (reduction(+:Fps) reduction(+:Fmus), nowarning) {

        // if cell contains boundary intercept
        if (ibm[] > 0. && ibm[] < 1.) {
            coord midPoint, b;
            coord n = facet_normal (point, ibm, ibmf);
            double alpha = plane_alpha (ibm[], n);
            double area = plane_area_center (n, alpha, &b) * pow(Delta, dimension-1);

            coord cellCenter = {x,y,z};
            foreach_dimension() {
                midPoint.x = cellCenter.x + b.x*Delta;
                n.x *= -1;
            }
            normalize(&n);

            // calculate pressure force
            double boundaryPressure = extrapolate_scalar (point, ibm, midPoint, n, p);
            double Fn = area * boundaryPressure;

            foreach_dimension()
                Fps.x -= Fn * n.x;

            // calculate shear force
            if (constant(mu.x) != 0.) {
            	double mua = 0., fa = 0.;

            	foreach_dimension() {
                    mua += mu.x[] + mu.x[1];
                    fa  += fm.x[] + fm.x[1];
	            }

                mua /= (fa + SEPS);
                coord velocityGrad = ibm_gradient (point, u, b, n);

#if dimension == 2
                foreach_dimension()
                    Fmus.x -= area * mua* (velocityGrad.x * (sq(n.x) + 1.) + 
                                           velocityGrad.y * -n.x * -n.y);
#else
                foreach_dimension()
                    Fmus.x -= area * mua * (velocityGrad.x * (sq(n.x) + 1.) + 
                                            velocityGrad.y * -n.x * -n.y +
                                            velocityGrad.z * -n.x * -n.z);
#endif
            }
        }
    }
    *Fp = Fps;
    *Fmu = Fmus;
}


double ibm_skin_friction (vector u, face vector mu, scalar cf)
{
    double cftotal = 0;
    foreach (reduction(+:cftotal)) {
        if (ibm[] > 0 && ibm[] < 1) {
            coord b;
            coord n = facet_normal (point, ibm, ibmf);
            double alpha = plane_alpha (ibm[], n);
            plane_area_center (n, alpha, &b) * pow(Delta, dimension-1);
            coord cellCenter = {x,y,z};
            foreach_dimension()
                n.x *= -1;
            normalize(&n);

            // calculate shear force
            double mua = 0., fa = 0.;

            foreach_dimension() {
                mua += mu.x[] + mu.x[1];
                fa  += fm.x[] + fm.x[1];
	        }

            mua /= (fa + SEPS);
            coord dudn = ibm_gradient (point, u, b, n);
            coord tau = {0};
           
#if dimension == 2
            coord tt = {-n.y, n.x}; // tangent vector
            coord dudt = ibm_gradient (point, u, b, tt);
            foreach_dimension()
                tau.x -= mua* (dudn.x * (sq(n.x) + 1.) + 
                               dudn.y * -n.x * -n.y);
            cf[] = distance(tau.x, tau.y);
#else
            foreach_dimension()
                tau.x -= mua * (dudn.x * (sq(n.x) + 1.) + 
                                dudn.y * -n.x * -n.y +
                                dudn.z * -n.x * -n.z);
            cf[] = distance3D(tau.x, tau.y, tau.z);
#endif
            cftotal += cf[];
        }
        else
            cf[] = 0;
    }
    return cftotal;
}

double ibm_pressure (Point point, scalar pressure, coord normal, coord markerPoint)
{
    return extrapolate_scalar (point, ibm, markerPoint, normal, pressure);
}

/*
double ibm_vorticityv2 (Point point, vector u, coord p, coord n) // needs improvement
{
    coord dudn = ibm_gradientv2 (point, u, p, n);

    return dudn.y*n.x - dudn.x*n.y;
}
*/
