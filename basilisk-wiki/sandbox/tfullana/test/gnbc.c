#include "navier-stokes/centered.h"
#include "contact.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/perfs.h"
#include "view.h"

#define GNBC 1
#define SLIP 1

double rho_1 = 1. [0], rho_2 = 0.2 [0];
double grav = -1.25 [0];
double sig = 1. [0];
double l_c = 1. [0];

double Ca = 0.04 [0];
double eps = 0.05 [0];
double angle_eq = pi/2. [0];

#define V_s sqrt(Ca)
#define mu_1 (V_s)
#define mu_2 (V_s)

#define L (10*l_c)
#define h_0 (3.5*l_c)

#define Re (rho_1*V_s*l_c/mu_1)

#define TCHAR (l_c/V_s)
#define TEND (20.*TCHAR)

int LEVEL = 11;
double f_err = 1e-3 [0], u_err = 1e-4 [0];

vector h[];

#define robin_local(b, c) ((dirichlet ((c)*Delta/(2*(b) + Delta))) + ((neumann (0))*((2*(b) - Delta)/(2*(b) + Delta) + 1.)))

double bell(double y, double eps){
	return ((1. - pow(tanh(y/eps), 2))/eps);
}

double Young_stress(double y, double eps, double theta, double theta_eq){
	return bell(y, eps)*(cos(theta_eq) - cos(theta));
}

double mybiquadratic(Point point, scalar s, coord p, int offset){
  double dx = p.x;
  double dy = p.y;
  double X[3] = {dx*dx,dx,1.};
  double Y[3] = {dy*dy,dy,1.};

double Bm1[9] = {0.5, -1, 0.5,
                -0.5,  0, 0.5,
                  0,  1,   0};
double Bm1T[9]= {0.5,-0.5, 0,
                  -1,  0,  1,
                 0.5,  0.5, 0};
double XX[3],YY[3];
  for (int kk=0;kk<=2;kk++){
    XX[kk] = Bm1[kk]*X[0] + Bm1[3+kk]*X[1] +
      Bm1[6+kk]*X[2];
    YY[kk] = Bm1T[kk*3]*Y[0] + Bm1T[kk*3+1]*Y[1] +
      Bm1T[kk*3+2]*Y[2];
  }


  for (int kk = 0;kk<=2;kk++){
    X[kk] = s[-1,-1+kk,offset]*XX[0] +
            s[ 0,-1+kk,offset]*XX[1] +
            s[+1,-1+kk,offset]*XX[2];
  }
  return X[0]*YY[0]+X[1]*YY[1]+X[2]*YY[2];
}

double linearInterpolation(double x1, double f_x1, double x2, double f_x2,
  double x)
{
  double result = (x - x1)/(x2-x1)*f_x2  + (x2-x)/(x2-x1)*f_x1;
  return result;
}

double quadraticInterpolation(double x0, double f_x0,
                              double x1, double f_x1,
                              double x2, double f_x2,
                              double x)
{
    double L0 = ((x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2));
    double L1 = ((x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2));
    double L2 = ((x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1));
    double result = f_x0 * L0 + f_x1 * L1 + f_x2 * L2;
    return result;
}

#if GNBC
scalar gnbc[];
double dynamic_angle = pi/2.;
h.t[left] = contact_angle(dynamic_angle);
#if SLIP
u.t[left] = robin_local(eps, gnbc[]);
#else
u.t[left] = dirichlet(gnbc[]);
#endif
#else
h.t[left] = contact_angle(angle_eq);
#if SLIP
u.t[left] = robin_local(eps, V_s);
#else
u.t[left] = dirichlet(V_s);
#endif
#endif

u.n[bottom]  = dirichlet(0.);
u.t[bottom] = dirichlet(0.);

u.n[right]  = neumann(0.);
u.t[right] = neumann(0.);

u.n[top]  = neumann(0.);
u.t[top] = neumann(0.);

u.n[left]  = dirichlet(0.);

vertex scalar psi[];

int main (int argc, char *argv[]) {
  if (argc > 1)
      Ca = atof(argv[1]);
  if (argc > 2)
      LEVEL = atoi(argv[2]);
  if (argc > 3)
      eps = atof(argv[3]);
  if (argc > 4)
      angle_eq = (pi/180.)*atof(argv[4]);

  origin(0., 0.);
  init_grid (128);
  DT = HUGE [0];
  f.sigma = sig;
  f.height = h;  
  rho1 = rho_1; rho2 = rho_2;
  mu1 = mu_1; mu2 = mu_2;
  size (L); 
  run();
}

event init (t = 0) {
if (!restore (file = "restart")) {
	fraction (f, -(y - h_0));
	double delta = L/pow(2,LEVEL);
	fprintf(ferr, "Re = %g Ca = %g GNBC = %d SLIP = %d LEVEL = %d\n", Re, Ca, GNBC, SLIP, LEVEL);
	fprintf(ferr, "Delta = %g eps = %g ep/Delta = %g, angle_eq = %lf\n", delta, eps, eps/delta, 180.*angle_eq/pi);
	fprintf(ferr, "--------------------------------------------------\n");
} else {
	double delta = L/pow(2,LEVEL);
	fprintf(ferr, "Re = %g Ca = %g GNBC = %d SLIP = %d LEVEL = %d\n", Re, Ca, GNBC, SLIP, LEVEL);
	fprintf(ferr, "Delta = %g eps = %g eps/Delta = %g, angle_eq = %lf\n", delta, eps, eps/delta, 180.*angle_eq/pi);
	fprintf(ferr, "--------------------------------------------------\n");
}
}

event acceleration(i++){
	face vector av = a;
	foreach_face(y)
    av.y[] += grav;
}

event adapt(i++){
	adapt_wavelet ({f,u}, (double[]){f_err,u_err,u_err,u_err},  LEVEL);
}

event extract (i++){
	char name[400];
	sprintf (name, "GNBC%d_SLIP%d_Ca%g_LEVEL%d_eps%g_angle%g", GNBC, SLIP, Ca, LEVEL, eps, angle_eq*180./pi);
	static FILE * fp = fopen(name, "w");
	double angle_extracted = 0., angle_above = 0., angle_extrapolated = 0., position = 0., maxYS = 1e-300;
	double curv = 0., Ca_loc = 0., Ca_biq0 = 0., Ca_biq = 0.;
	scalar kappa[];

	curvature (f,kappa);

	foreach(reduction(max:angle_extrapolated) reduction(max:angle_above)){
		if ((x <= (2.*Delta)) & (x > (Delta)) & (h.y[] != nodata) & (f[]< (1. - 1e-6)) & (f[]> 1e-6)){
				coord n = {0};
				n = interface_normal (point,f);
				angle_above = atan2(n.y, n.x);
				double tmp = pow(((h.y[1,0] - h.y[])/1.0), 2);
				angle_extrapolated = angle_above + 1.5*Delta*kappa[]*sqrt(1 + tmp)/sin(angle_above);
		}
	}
	scalar uy_tmp[];
	foreach()
		uy_tmp[] = u.y[];
	foreach_boundary(left, reduction(max:position) reduction(max:angle_extracted) reduction(max:curv) reduction(max:Ca_loc) reduction(max:Ca_biq0) reduction(max:Ca_biq)){
		if ((h.y[] != nodata) & (f[]< (1. - 1e-6)) & (f[]> (1e-6))){
			double loc_position = Delta*height(h.y[]);
			position = y + loc_position;
			coord n = {0};
			n = interface_normal (point,f);
			angle_extracted = atan2(n.y, n.x);
			curv = kappa[];
			// coord p_interp = {0., position};
			// Ca_loc = (mu_1/sig)*(((u.y[] - u.y[0,-1])/Delta)*(Delta*height(h.y[]) - Delta/2.) + u.y[]); 
			Ca_loc = mu_1*u.y[]/sig; // (u.y[0,1] + 2.*u.y[] + u.y[0,-1])/4.
			double _uy0 = 0.5*(u.y[0,-1] + u.y[-1,-1]);
			double _uy1 = 0.5*(u.y[0,0] + u.y[-1,0]);
			// double uy1 = linearInterpolation(-Delta/2., u.y[0,-1], Delta/2., u.y[0,0], loc_position);
			// double uy1 = linearInterpolation(-Delta/2., _uy0, Delta/2., _uy1, loc_position);
			double uy0 = quadraticInterpolation(-5.*Delta/2., u.y[0,-3], -3.*Delta/2., u.y[0,-2], -Delta/2., u.y[0,-1], loc_position);
			double uy1 = quadraticInterpolation(-3.*Delta/2., u.y[0,-2], -Delta/2., u.y[0,-1], Delta/2., u.y[0,0], loc_position);
			// double uy0 = linearInterpolation(-Delta/2., u.y[-1,-1], Delta/2., u.y[-1,0], loc_position);
			// Ca_biq = (mu_1/sig)*mybiquadratic(point, uy_tmp, p_interp, 0);
			Ca_biq0 = (mu_1/sig)*uy0;
			Ca_biq = (mu_1/sig)*uy1;
		}
	}

	double _x2 = 0., _y2 = 0., _u2 = 0., _x1 = 0., _y1 = 0., _u1 = 0., _x0 = 0., _y0 = 0., _u0 = 0.;
	foreach (reduction(max:_x2) reduction(max:_y2) reduction(max:_u2)){
		if ((h.y[] != nodata) & (f[]< (1. - 1e-6)) & (f[]> (1e-6)) & (x >= 2.1*Delta) & (x <= 2.9*Delta)){
			double tmp_pos = y + Delta*height(h.y[]);
			_x2 = x;
			_y2 = tmp_pos;
			_u2 = u.y[];
		}
	}
	foreach (reduction(max:_x1) reduction(max:_y1) reduction(max:_u1)){
		if ((h.y[] != nodata) & (f[]< (1. - 1e-6)) & (f[]> (1e-6)) & (x >= 1.1*Delta) & (x <= 1.9*Delta)){
			double tmp_pos = y + Delta*height(h.y[]);
			_x1 = x;
			_y1 = tmp_pos;
			_u1 = u.y[];
		}
	}
	foreach_boundary (left, reduction(max:_x0) reduction(max:_y0) reduction(max:_u0)){
		if ((h.y[] != nodata) & (f[]< (1. - 1e-6)) & (f[]> (1e-6))){
			double tmp_pos = y + Delta*height(h.y[]);
			_x0 = x;
			_y0 = tmp_pos;
			_u0 = u.y[];
		}
	}
	double _Delta = L/pow(2, LEVEL);
	double Ca_int = (mu_1/sig)*quadraticInterpolation(5*_Delta/2., _u2, 3.*_Delta/2., _u1, _Delta/2., _u0, 0.);
	
	#if GNBC
	foreach_boundary(left, reduction(max:maxYS)){
		double relative_position = y - position;
		gnbc[] = V_s + (eps)*(sig/mu_1)*Young_stress(relative_position, eps, angle_extracted, angle_eq);
		if (fabs(maxYS) < fabs(gnbc[]))
			maxYS = gnbc[];
	}
	dynamic_angle = angle_extrapolated;
	#endif

	fflush(fp);
	fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", t/TCHAR, (position-h_0)/l_c, 180.*angle_above/pi, 180.*angle_extrapolated/pi, 180.*angle_extracted/pi, Ca_loc, curv, maxYS, Ca_biq0, Ca_biq, Ca_int);
}

event compute_stream (i++){
  scalar omega[];
  scalar _psi[];

  vorticity(u, omega);

  _psi[bottom] = dirichlet (0.);
  _psi[top]    = neumann (0.);
  _psi[left]  = dirichlet (0.);
  _psi[right]   = neumann (0.);

  poisson (_psi, omega);
  boundary ({_psi});

  foreach_vertex()
    psi[] = (_psi[0,-1] + _psi[-1,-1] + _psi[] + _psi[-1])/4.;
}

event snapshot (t = 0; t+=0.1*TCHAR; t <= TEND){
  char name[400];
	sprintf (name, "s-GNBC%d_SLIP%d_Ca%g_LEVEL%d_eps%g_angle%g-%g", GNBC, SLIP, Ca, LEVEL, eps, angle_eq*180./pi, t/TCHAR);
	dump (name);
}

event final (t = TEND){
	char name[400], name2[400], name3[400], name4[400];
	sprintf (name, "GNBC%d_SLIP%d_pressure_Ca%g_LEVEL%d_eps%g_angle%g_%d", GNBC, SLIP, Ca, LEVEL, eps, angle_eq*180./pi, tid());
	sprintf (name2, "GNBC%d_SLIP%d_curvature_Ca%g_LEVEL%d_eps%g_angle%g_%d", GNBC, SLIP, Ca, LEVEL, eps, angle_eq*180./pi, tid());
	sprintf (name3, "GNBC%d_SLIP%d_facets_Ca%g_LEVEL%d_eps%g_angle%g_%d", GNBC, SLIP, Ca, LEVEL, eps, angle_eq*180./pi, tid());
	sprintf (name4, "GNBC%d_SLIP%d_dudy_Ca%g_LEVEL%d_eps%g_angle%g_%d", GNBC, SLIP, Ca, LEVEL, eps, angle_eq*180./pi, tid());
	FILE * fp = fopen(name, "w");
	FILE * fp2 = fopen(name2, "w");
	FILE * fp3 = fopen(name3, "w");
	FILE * fp4 = fopen(name4, "w");
	scalar kappa[];
	double position = 0.;

	curvature (f,kappa);

	foreach_boundary(left, reduction(max:position)){
		if ((h.y[] != nodata) & (f[]< (1. - 1e-6)) & (f[]> (1e-6))){
			position = y + Delta*height(h.y[]);
		}
		double dudy = (u.y[1,0] - u.y[])/Delta;
		fprintf(fp4,"%lf %lf\n", y, dudy);
	}

	foreach(){
		fprintf(fp,"%lf %lf %lf\n", x, y, p[]);
		if ((h.y[] != nodata) & (f[]< (1. - 1e-6)) & (f[]> (1e-6))){
			double ycl = y + Delta*height(h.y[]);
			double r = sqrt(pow(x,2) + pow((ycl - position), 2));
			coord n = {0};
			n = interface_normal (point,f);
			double theta = atan2(n.y, n.x);
			fprintf(fp2,"%lf %lf %lf\n", r, kappa[], 180.*theta/pi);
		}
	}

	output_facets (f, fp3);

	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
}