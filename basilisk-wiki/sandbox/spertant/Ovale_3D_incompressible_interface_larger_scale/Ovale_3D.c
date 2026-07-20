/**  

# Capillary relaxation of a soap film with two thicknesses separated by an oval shape (with smooth transition with hyperbolic tangent) with incompressible free surface.

This test case includes incompressible free surface and surface viscosity $\eta_s$. Using the lubrication approximation and especially small slopes assumption, the interface incompressibility conditions reads 
$$
\nabla \cdot \bm{u}_s = 0, 
$$
where $\nabla$ is the 2D operator (consistent with multilayer solver) and $\bm{u}_s$ is the horizontal velocity at the free surface.
The Marangoni constraint at the free surface can be written
$$
-  \left. \rho \, \nu \, \frac{\partial \bm{u}}{\partial z} \right|_{z=\eta} + \eta_s \Delta  \bm{u}_s + \nabla \sigma = 0.
$$
Using the lubrication theory, the bulk viscosity term can be transformed to obtain:
$$
\sigma_0 \, \eta \nabla \Delta \eta + \eta_s \Delta  \bm{u}_s + \nabla \sigma = 0,
$$
Where $\sigma_0$ is the surface tension of the surface at rest.
The free surface velocity $\bm{u}_s$ is thus the solution of a 2D Stokes flow defined by the equations 
$$
\begin{aligned}
\eta_s \Delta  \bm{u}_s + \nabla \sigma & = - \sigma_0 \, \eta \nabla \Delta \eta, \\
\nabla \cdot \bm{u}_s & = 0, 
\end{aligned}
$$
where the surface tension plays the role of a pressure and the right-hand side has no 3D equivalent, but can be estimated knowing the free-surface height. 
The free surface velocity $\bm{u}_s$ being determined, it can be used as a boundary condition for the layered bulk velocity solved by the multilayered solver.

Taking the divergence of the Marangoni equation and using the fact that Laplacian operator and Divergence operator can be swapped and using the incompressibility condition, we obtain a Poisson equation for the surface tension :
$$
\Delta \sigma = - \nabla \cdot \left( \sigma_0 \, \eta \nabla \Delta \eta \right),
$$
which can be solved with Dirichlet conditions $\sigma = \sigma_0$ on the boundaries if the film is at rest on the boundaries.

Using again the Marangoni equation, we have a vectorial Poisson equation for the free surface velocity :
$$
\Delta  \bm{u}_s = - \frac{1}{\eta_s} \left( \nabla \sigma + \sigma_0 \, \eta \nabla \Delta \eta \right),
$$
which is splitted into two scalar Poisson equations to solve:
$$
\begin{aligned}
\Delta  \bm{u}_{s,x} & = - \frac{1}{\eta_s} \left( \frac{\partial \sigma}{\partial x} + \sigma_0 \, \eta \frac{\partial \left( \Delta \eta \right) } {\partial x} \right), \\
\Delta  \bm{u}_{s,y} & = - \frac{1}{\eta_s} \left( \frac{\partial \sigma}{\partial y} + \sigma_0 \, \eta \frac{\partial \left( \Delta \eta \right) } {\partial y} \right). \\
\end{aligned}
$$
We thus have a system of equations verifying $\Delta \left( \nabla \cdot \bm{u}_s \right) = 0$. 
If we ensure $\nabla \cdot \bm{u}_s = 0$ on the domain boundaries (rectangle sides here), the unique solution of the problem is $\nabla \cdot \bm{u}_s = 0$ everywhere (To be proved with maximum principle ?).

The boundary condition imposed for $\bm{u}_s$ has to satisfy the condition on each boundary
$$
\frac{\partial {\bm{u}_s,n}}{\partial n} + \frac{\partial {\bm{u}_s,t}}{\partial t} = 0,
$$
where $n$ and $t$ are the normal and tangential coordinates associated with each boundary.
Arbitrarily, we impose the set of conditions
$$
\begin{aligned}
\frac{\partial {\bm{u}_s,n}}{\partial n} &  = 0 \\ ,
\bm{u}_s,t & = 0,
\end{aligned}
$$
which respects the divergence-free condition on the boundaries.

*/


/**
## Include and parameters
This test uses the multilayer solver, with surface tension and viscosity and without gravity. */
#include "grid/multigrid.h"
#include "../Solver_modified/hydro-tension.h"
#include "layered/remap.h"
#include "poisson.h"
#include "view.h"

double h_0 = 0.3e-6; // 1e-6;  // film at rest height

#define nu0 1.5e-6 // 1e-6 // mu/rho
#define eta_s0 1e-11 // 1e-8 // 1e-10 // 1e-11 // eta_s/rho 
#define sigma0 3.5e-5 // 2.0e-5 // gamma/rho
#define L 2e-2 //4e-2 //2e-2 //
#define T_END 1.0   
#define C_F_L 0.5
#define d 1.5e-4 //width of the initial thickness profile transition
#define c 20. // 3.5  //ratio of oval center-center distance with semi width of oval
#define a 5e-4 //semi width of the oval
#define delta_h 20.0 //0.3 //initial thickness difference
#define domain_ratio 2.  // ratio Ly/Lx of domain box
double xc = 0; //oval center horizontal position

/**
## Main function
*/
int main()
{
  L0 = L;
  dimensions(nx = 1., ny = domain_ratio); // Il faut que les dimensions soient superieures a 1 pour que les valeurs soient prises en compte. L'echelle se fait dans origin.
  origin(-L0/2, -domain_ratio*L0/2);
  G = 0.;
  nu = nu0;
  eta_s = eta_s0;
  N = 512; // 1024 // 512;
  nl = 1; 

  symmetric_bathymetry = true;
  navier_surface = true;

  run();
 
}

/**
## Initialisation 
*/
scalar Hi[], div_us[], div_us_nodal[], div_lapl_us[], div_rhs[], rhs_poisson[], error_poisson_sigma[], sigma_poisson[], error_poisson_us_x[], error_poisson_us_y[];
vector rhs_poisson_us[];
double maxerror;
event init (i = 0)
{
  // Il faut modifier les CFL ici sinon la modif sur CFL (et pas CFL_H) n'est pas prise en compte	
  CFL_H = C_F_L; 
  CFL = C_F_L;

  // Il faut mettre le foreach_dimension() pour bien faire le tour des boundaries
  foreach_dimension(){

    // LES CONDITIONS LIMITES SUR U DOIVENT ETRE LES MEME QUE CELLES SUR US	  
    u.n[right] = neumann(0);
    u.n[left] = neumann(0);
    
    u.t[right] = dirichlet(0);
    u.t[left] = dirichlet(0);
    
    sigma[right] = dirichlet (sigma0);
    sigma[left] = dirichlet (sigma0);

    sigma_poisson[right] = dirichlet (sigma0);
    sigma_poisson[left] = dirichlet (sigma0);
  
    us_poisson.n[right] = neumann (0);
    us_poisson.n[left] = neumann (0);
    
    us_poisson.t[right] = dirichlet (0);
    us_poisson.t[left] = dirichlet (0);

    us.n[right] = neumann (0);
    us.n[left] = neumann (0);

    us.t[right] = dirichlet (0);
    us.t[left] = dirichlet (0);
  }

  foreach(){ 
    double H = h_0;
    double r = fabs(x-xc)-a;
    if (y>c*a)
      r = sqrt(sq(x-xc) + sq(y-c*a))-a;	    
    if (y<-c*a)
      r = sqrt(sq(x-xc) + sq(y+c*a))-a;
    if (r < 3*d)    
      H = h_0 + delta_h*h_0*(1.+tanh(-r/d))/2.;
    Hi[] = H;
    sigma_poisson[] = sigma0;
    sigma[] = sigma0;
    foreach_dimension()
      us_poisson.x[] = 0;
    foreach_layer() {
      h[] = H/nl;
    }
  }
  boundary({sigma,sigma_poisson});
}


/**
## Outputs 
*/


event logfile (i+=10; i = 0; t <= T_END)
//event logfile (i++, i<2)
{
  char name[250];
  sprintf (name, "height_center_gamma_with_new_eta-h0-%g-L0-%g-nl-%d-N-%d-sigma0-%g-nu-%g-CFL-%g-d-%g-dh-%g-a-%g-c-%g-etas-%g-IMPOSE_US_POISSON_NON_SQUARE", h_0, L0, nl, N, sigma0, nu, CFL, d, delta_h, a, c, eta_s);
  static FILE * fph = fopen (name, "w");
  if (i > 0)
    fprintf(fph, "%g %g %g %g %d\n", t, statsf(eta).sum, statsf(eta).max, dt, i);
}

// ON DOIT D'ABORD PASSER DANS VISCOUS TERM ICI, PUIS DANS HYDRO-TENSION, PUIS DANS DIFFUSION (CE QUI EST FAIT)

event viscous_term (i++) 
{
  scalar lapl_eta[];
  face vector eta_grad_lapl_eta[];
  face vector grad_sigma[];

  foreach() 
    lapl_eta[] = (eta[1,0]-2*eta[0,0]+eta[-1,0])/sq(Delta) + (eta[0,1]-2*eta[0,0]+eta[0,-1])/sq(Delta); 	    
  foreach_face(){
    eta_grad_lapl_eta.x[] = (eta[-1,0]+eta[0,0])/2.*(lapl_eta[0,0]-lapl_eta[-1,0])/Delta;
  }	  

  foreach(){
    rhs_poisson[] = 0;
    foreach_dimension(){
      rhs_poisson[] += -sigma0*(eta_grad_lapl_eta.x[1,0]-eta_grad_lapl_eta.x[0,0])/Delta;
    }
  }	  
  
  TOLERANCE = 1e-9; //1e-5;
  poisson (sigma_poisson, rhs_poisson);
  
  foreach()
    sigma[] = sigma_poisson[];	  
   
  foreach_face(){
    grad_sigma.x[] = (sigma[0,0]-sigma[-1,0])/Delta;
  }

  foreach() 
    foreach_dimension() 
      rhs_poisson_us.x[] = -1./eta_s*(grad_sigma.x[]+sigma0*eta_grad_lapl_eta.x[]);

  TOLERANCE = 1e-5;

  poisson (us_poisson.x, rhs_poisson_us.x);
  
  poisson (us_poisson.y, rhs_poisson_us.y);

  foreach(){
    div_us[] = 0;
    div_lapl_us[] = 0;
    foreach_dimension() {
      div_us[] += (us_poisson.x[1,0]-us_poisson.x[0,0])/Delta;
      us.x[] = (us_poisson.x[0,0]+us_poisson.x[1,0])/2.0;
    }
    div_lapl_us[] = (us_poisson.x[2,0]-5*us_poisson.x[1,0]-us_poisson.x[-1,0]-us_poisson.x[0,1]-us_poisson.x[0,-1]+us_poisson.x[1,1]+us_poisson.x[1,-1]+5*us_poisson.x[0,0])/pow(Delta,3) + (us_poisson.y[0,2]-5*us_poisson.y[0,1]-us_poisson.y[0,-1]-us_poisson.y[1,0]-us_poisson.y[-1,0]+us_poisson.y[1,1]+us_poisson.y[-1,1]+5*us_poisson.y[0,0])/pow(Delta,3);
    div_rhs[] = -1./eta_s*((sigma[1,0]+sigma[0,1]-4*sigma[0,0]+sigma[-1,0]+sigma[0,-1])/sq(Delta)-rhs_poisson[]); 
    error_poisson_sigma[] = (sigma[1,0]+sigma[0,1]-4*sigma[0,0]+sigma[-1,0]+sigma[0,-1])/sq(Delta)-rhs_poisson[];
    error_poisson_us_x[] =  (us_poisson.x[1,0]+us_poisson.x[0,1]-4*us_poisson.x[0,0]+us_poisson.x[-1,0]+us_poisson.x[0,-1])/sq(Delta)-rhs_poisson_us.x[];
    error_poisson_us_y[] =  (us_poisson.y[1,0]+us_poisson.y[0,1]-4*us_poisson.y[0,0]+us_poisson.y[-1,0]+us_poisson.y[0,-1])/sq(Delta)-rhs_poisson_us.y[];
  }

  foreach(){
    div_us_nodal[] = 0;
    foreach_dimension()
      div_us_nodal[] += (us.x[1,0]-us.x[-1,0])/(2*Delta);
  }
}

event save_us_max (i+=10)
//event save_us_max (i++)
{
  double div_us_avg = 0, div_us_max = 0, x_div_us_max = 0, y_div_us_max = 0, area_thick_film = 0;
  double y_top=0, y_bottom=0, x_left=0, x_right=0, L_Y=0, L_X=0;
  int count = 0;  
  foreach() {
    div_us_avg += div_us_nodal[];
    count += 1;
    if (fabs(div_us_nodal[]) > div_us_max){
      div_us_max = div_us_nodal[];
      x_div_us_max = x;
      y_div_us_max = y;
    }
    if (eta[]>(2+delta_h)*h_0/2.){
      area_thick_film += sq(Delta);
    }
    if (fabs(x)<xc+L0/N && x>xc){
      if ((eta[0,0]>(2+delta_h)*h_0/2.) && (eta[0,1]<(2+delta_h)*h_0/2.))
        y_top = y;	      
      if ((eta[0,0]<(2+delta_h)*h_0/2.) && (eta[0,1]>(2+delta_h)*h_0/2.))
        y_bottom = y;	      
    }
    if (fabs(y)<L0/N && y>0){
      if ((eta[0,0]>(2+delta_h)*h_0/2.) && (eta[1,0]<(2+delta_h)*h_0/2.))
        x_right = x;	      
      if ((eta[0,0]<(2+delta_h)*h_0/2.) && (eta[1,0]>(2+delta_h)*h_0/2.))
        x_left = x;	      
    }
    L_X = x_right - x_left;
    L_Y = y_top - y_bottom;
  } 
  div_us_avg = div_us_avg/count;
  char name3[250];
  sprintf (name3, "us-max_gamma_with_new_eta-h0-%g-L0-%g-nl-%d-N-%d-sigma0-%g-nu-%g-CFL-%g-d-%g-dh-%g-a-%g-c-%g-etas-%g-IMPOSE_US_POISSON_NON_SQUARE", h_0, L0, nl, N, sigma0, nu, CFL, d, delta_h, a, c, eta_s);
  static FILE * file = fopen (name3, "w");
  fprintf(file, "%d %g %g %g %g %g %g %g %g %g %g\n", i, t, statsf(us.x).max, statsf(us.y).max, statsf(div_us_nodal).max, div_us_avg, x_div_us_max, y_div_us_max, area_thick_film, L_X, L_Y);
}

event div_us_x_profile (i+=100)
//event div_us_x_profile (i++)
{
  char name2[250];
  sprintf (name2, "div_us_x_profile_gamma_with_new_eta-h0-%g-L0-%g-nl-%d-N-%d-sigma0-%g-nu-%g-CFL-%g-d-%g-dh-%g-a-%g-c-%g-etas-%g-IMPOSE_US_POISSON_NON_SQUARE", h_0, L0, nl, N, sigma0, nu, CFL, d, delta_h, a, c, eta_s);
  static FILE * fpc = fopen (name2, "w");
  foreach(){
    if (fabs(y)<L0/N && y>0){
      fprintf (fpc, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", i, t, x, eta[], us_poisson.x[], us_poisson.y[], div_lapl_us[], div_rhs[], error_poisson_us_x[], error_poisson_us_y[], rhs_poisson_us.x[], rhs_poisson_us.y[], error_poisson_sigma[], rhs_poisson[], div_us[], (us_poisson.x[1,0] - us_poisson.x[0,0])/(Delta), (us_poisson.y[0,1] - us_poisson.y[0,0])/(Delta), div_us_nodal[], (us.x[1,0] - us.x[-1,0])/(2*Delta), (us.y[0,1] - us.y[0,-1])/(2*Delta), sigma[]);
    }  
  }

  fprintf(fpc, "\n\n");
}

event div_us_y_profile (i+=100)
//event div_us_y_profile (i++)
{
  char name5[250];
  sprintf (name5, "div_us_y_profile_gamma_with_new_eta-h0-%g-L0-%g-nl-%d-N-%d-sigma0-%g-nu-%g-CFL-%g-d-%g-dh-%g-a-%g-c-%g-etas-%g-IMPOSE_US_POISSON_NON_SQUARE", h_0, L0, nl, N, sigma0, nu, CFL, d, delta_h, a, c, eta_s);
  static FILE * fpz = fopen (name5, "w");
  foreach(){
    if (fabs(x)<xc+L0/N && x>xc)	  
      fprintf (fpz, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", i, t, y, eta[], us_poisson.x[], us_poisson.y[], div_lapl_us[], div_rhs[], error_poisson_us_x[], error_poisson_us_y[], rhs_poisson_us.x[], rhs_poisson_us.y[], error_poisson_sigma[], rhs_poisson[], div_us[], (us_poisson.x[1,0] - us_poisson.x[0,0])/(Delta), (us_poisson.y[0,1] - us_poisson.y[0,0])/(Delta), div_us_nodal[], (us.x[1,0] - us.x[-1,0])/(2*Delta), (us.y[0,1] - us.y[0,-1])/(2*Delta), sigma[]);
  }
  fprintf(fpz, "\n\n");
}

//event bnd_profile (i++)
event bnd_profile (i+=100)
{
  char name10[250];
  sprintf (name10, "bnd_profile_gamma_with_new_eta-h0-%g-L0-%g-nl-%d-N-%d-sigma0-%g-nu-%g-CFL-%g-d-%g-dh-%g-a-%g-c-%g-etas-%g-IMPOSE_US_POISSON_NON_SQUARE", h_0, L0, nl, N, sigma0, nu, CFL, d, delta_h, a, c, eta_s);
  static FILE * file10 = fopen (name10, "w");

  foreach_boundary(bottom){
    fprintf (file10, "%d %g %g %g %g %g %g %g %g %g\n", i, t, y, x, us.x[], us.y[], (us.y[0,1]-us.y[0,0])/Delta, (us.y[0,0]-us.y[0,-1])/Delta, div_us_nodal[], div_us[]);
  }
  fprintf(file10, "\n");
  
  foreach_boundary(top){
    fprintf (file10, "%d %g %g %g %g %g %g %g %g %g\n", i, t, y, x, us.x[], us.y[], (us.y[0,0]-us.y[0,-1])/Delta, (us.y[0,1]-us.y[0,0])/Delta, div_us_nodal[], div_us[]);
  }
  fprintf(file10, "\n");

  foreach_boundary(left){
    fprintf (file10, "%d %g %g %g %g %g %g %g %g %g\n", i, t, y, x, us.x[], us.y[], (us.x[1,0]-us.x[0,0])/Delta, (us.x[0,0]-us.x[-1,0])/Delta, div_us_nodal[], div_us[]);
  }
  fprintf(file10, "\n");

  foreach_boundary(right){
    fprintf (file10, "%d %g %g %g %g %g %g %g %g %g\n", i, t, y, x, us.x[], us.y[], (us.x[0,0]-us.x[-1,0])/Delta, (us.x[1,0]-us.x[0,0])/Delta, div_us_nodal[], div_us[]);
  }

  fprintf(file10, "\n\n");
}

void plot_profile_x (double t, FILE * fp)
{
  fprintf (fp,
	   "set title 't = %e'\n"
	   "plot [%g:%g][%g:%g]'-' u 1:2 w l lc rgb 'red' t ''\n", t, -L0/2, L0/2, 0.0, (3.0+delta_h)*h_0);
  foreach(){
    if (fabs(y)<L0/N && y>0)	  
      fprintf (fp, "%g %g\n", x, eta[]);
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
}


void plot_profile_y (double t, FILE * fy)
{
  fprintf (fy,
	   "set title 't = %e'\n"
	   "plot [%g:%g][%g:%g]'-' u 1:2 w l lc rgb 'red' t ''\n", t, -domain_ratio*L0/2, domain_ratio*L0/2, h_0/2., (3.+delta_h)*h_0);
  foreach(){
    if (fabs(x)<xc+L0/N && x>xc)	  
      fprintf (fy, "%g %g\n", y, eta[]);
  }
  fprintf (fy, "e\n\n");
  fflush (fy);
}

event animate_profile_x (i += 10)
{ 
static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
if (i == 0)
   fprintf (fp, "set term x11\n");
plot_profile_x (t, fp);
}
event animate_profile_y (i += 10)
{ 
static FILE * fy = popen ("gnuplot 2> /dev/null", "w");
if (i == 0)
   fprintf (fy, "set term x11\n");
plot_profile_y (t, fy);
}

event images (i+=10)
{
  char name6[250];
  sprintf (name6, "thickness_field_gamma_with_new_eta-h0-%g-L0-%g-nl-%d-N-%d-sigma0-%g-nu-%g-CFL-%g-d-%g-dh-%g-a-%g-c-%g-etas-%g-IMPOSE_US_POISSON_NON_SQUARE.mp4", h_0, L0, nl, N, sigma0, nu, CFL, d, delta_h, a, c, eta_s);
  output_ppm(eta, file=name6, min = h_0/2., max = (3.+delta_h)*h_0);
}	

event movie (i+=10; t <= T_END)
{
  view (fov = 22, quat = {0.475152,0.161235,0.235565,0.832313},
	tx = 0.02, ty = 0., width = 800, height = 600);
  char s[80];
  sprintf (s, "t = %.7f", t);
  draw_string (s, size = 80);
  squares ("eta", linear = true, z = "eta", min = -0.3, max = 0.3);
  char name4[250];
  sprintf (name4, "profile_gamma_with_new_eta-h0-%g-L0-%g-nl-%d-N-%d-sigma0-%g-nu-%g-CFL-%g-d-%g-dh-%g-a-%g-c-%g-etas-%g-IMPOSE_US_POISSON_NON_SQUARE.mp4", h_0, L0, nl, N, sigma0, nu, CFL, d, delta_h, a, c, eta_s);
  save (name4);
}
