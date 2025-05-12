/**
# Clogging of a dual-plate channel by a droplet
*/
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "contact-embed.h"
#include "view.h"
#include "tag.h"
#include "maxruntime.h"
#include "dgeometry.h"

#define fnorm(v) (sqrt(sq(v.x) + sq(v.y)))

// Simulation parameters
double  tend      = 20.,            // End time
        uemax     = 0.1;            // Error on u field

int     iniLvl    = 4,              // Initial level of refinement
        iniDrpLvl = 5,              // Initial LOR for droplet
        iniGrdLvl = 6,              // Initial LOR for grid
        cellsInL  = 40,             // Min. cell number in channel width L
        maxLvl;

char    *restFile = "restore.dump", // Restore file name
        *mode     = "w";            // Writing mode: "w" or "a"

// Flow parameters
double  rhoRatio  = 829.15,         // Density ratio rho1/rho2
        muRatio   = 63.82,          // Dynamic viscosity ratio mu1/mu2
        Mo        = 1.94e-12,       // Morton, Mo = g*mu1^4/(rho1*sig^3)
        theta0    = 65,             // Contact angle (°)
        uair      = -0.1,           // Inflow velocity
        Bo        = 1;              // Bond, Bo = rho1*g*D^2/sig

// Geometric parameters
double  HoL       = 6.38,           // Aspect ratio H/L
        roL       = 0.3,            // Aspect ratio r/L
        DoL       = 2.,             // Aspect ratio D/L
        WoL       = 2.;             // Aspect ratio W/L


double channel (coord p) {
  coord p1 = {-1/2. - WoL,  -HoL};
  coord p2 = {-1/2.,        0.};
  coord p3 = {1/2.,         -HoL};
  coord p4 = {1/2. + WoL,   0.};
  return -1*union(curvedRectangle(p, p1, p2, roL),
      curvedRectangle(p, p3, p4, roL));
}

void generateDroplet (coord p1, double D) {
  /** This event create a droplet, test its volume
  and fuse it with the volume fraction field f.*/
      
  // Refine the area around the droplet
  for (int lev = iniLvl+1; lev<=iniDrpLvl; lev++)
    refine (sphere((coord) {x, y, z}, p1, max(Delta, D/2.)) > 0.
        && level < lev);

  // Refine the area at the boundary of the droplet
  scalar drpField[];
  drpField.refine = drpField.prolongation = fraction_refine;
  fraction (drpField, sphere((coord) {x, y, z}, p1, D/2.));

  for (int lev=iniDrpLvl+1; lev<=maxLvl; lev++) {
    refine ((drpField[] > 0. && drpField[] < 1.) && level < lev);
    fraction (drpField, sphere((coord) {x, y, z}, p1, D/2.));
  }

  foreach () {
    if (drpField[] > 0.)
      f[] = min(f[]+drpField[], 1.);
  }
}

void generateSolid () {
  // Solid creation
  // Initial refinement
  for (int lev=iniLvl+1; lev<=iniGrdLvl; lev++)
    refine (channel((coord) {x, y, z}) < 0. && level < lev);

  // Initial fraction    
  vertex scalar phi[];
  foreach_vertex()
    phi[] = channel((coord) {x, y, z});
  boundary ({phi});
  fractions (phi, cs, fs);

  // Refinement on the interface + fraction
  for (int lev=iniGrdLvl+1; lev<=maxLvl; lev++) {
    refine ((cs[] > 0. && cs[] < 1.) && level < lev);

    vertex scalar phi[];
    foreach_vertex()
      phi[] = channel((coord) {x, y, z});
    boundary ({phi});
    fractions (phi, cs, fs);
  }      
}

void plot (char name[]) {
  scalar f2[], unorm[];
  foreach () {
    f2[] = (f[] > 1e-2 && cs[] > 0.) ? f[] : 0.;
    unorm[] = norm(u);
  }

  translate (y = -Y0 - L0/2.) {
    draw_vof (c = "cs", filled = -1, 
        fc = {0.5215686274509804,0.3411764705882353,0.3411764705882353});
    draw_vof (c = "f2", filled = -1, fc = {1,1,1}, lc = {1,1,1});
    draw_vof (c = "f2", lw = 1);
    squares (color = "unorm", max = 10, min = 0.);   
    box ();
  }
  
  save (name);
}

// Inflow conditions
u.n[top] = dirichlet(uair);
p[top]   = neumann(0.);
pf[top]  = neumann(0.);
f[top]   = 0.;

// Outflow conditions
u.n[left]   = (u.x[] < 0) ? neumann(0.) : dirichlet(0.);
p[left]     = dirichlet(0.);
pf[left]    = dirichlet(0.);
f[left]     = (u.x[] < 0) ? neumann(0.) : 0.;

u.n[right] = (u.x[] > 0) ? neumann(0.) : dirichlet(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
f[right]   = (u.x[] > 0) ? neumann(0.) : 0.;

u.n[bottom] = (u.y[] < 0) ? neumann(0.) : dirichlet(0.);
p[bottom]   = dirichlet(0.);
pf[bottom]  = dirichlet(0.);
f[bottom]   = (u.y[] < 0) ? neumann(0.) : 0.;

// Embed boundary
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

int main(int argc, char *argv[]) {
  maxruntime(&argc, argv);

  // Modify default parameter values
  if (argv[1] != NULL && argv[2] != NULL 
      && argv[3] != NULL && argv[4] != NULL) {
    DoL       = 1./strtod(argv[1], NULL); 
    Bo        = strtod(argv[2], NULL);
    theta0    = strtod(argv[3], NULL);
    cellsInL  = atoi(argv[4]);
  }
  if (argv[5] != NULL) {
    restFile  = argv[5];
    mode      = "a";
  }

  const scalar c[] = theta0*pi/180.;
  contact_angle = c;

  // Geometric variable
  WoL = max(2.*DoL, 4.);    // Solid thickness in mm

  // Global variable
  DT = HUGE [0];
  L0 = 1.21[0]*max(2*WoL + 1, DoL + HoL + 2*DoL);
  X0 = -L0/2.;
  Y0 = -(HoL + 2*DoL + 0.1*L0);
  # if dimension == 3
  Z0 = -L0/2.;
  #endif
  maxLvl = (int) ceil( log(cellsInL*L0)/log(2.) );
  
  rho1 = 1.;
  rho2 = rho1/rhoRatio;
  G.y = -1.;
  f.sigma = (rho1 - rho2)*fnorm(G)*sq(DoL)/Bo;
  double La = sqrt(Bo/Mo);    // Laplace number
  mu1 = sqrt(f.sigma*rho1*DoL/La);
  mu2 = mu1/muRatio;

  init_grid (1 << iniLvl);

  // Run the simulation
  run();
}

event init(t = 0.) {
  // Unsuccessful restoration
  if (!restore (file = restFile)) {
    generateSolid ();

    // Droplet creation
    coord pDrp = {0., DoL/2.};
    generateDroplet(pDrp, DoL);
  }
  // Successful restauration
  else {
    generateSolid ();
  }
}

event adapt (i++) {
  scalar f1[], cs1[];
  foreach () {
    f1[] = f[];
    cs1[] = cs[];
  }	

  adapt_wavelet ((scalar*){f1,u,cs1}, 
      (double[]){0.001, uemax, uemax, 0.001}, maxLvl);
}

event plotEvent (t += 1e-1) {
  char file[100];
  sprintf (file, "%07.3f.jpg", t);
  plot (file);
}

event snapshot (t += 1) {
  char name[100];
  sprintf (name, "%07.3f.dump", t);
  dump (file = name);
}

event output (t += 1e-1, t <= tend) {
  scalar m[];
  foreach ()
    m[] = f[] > 1e-3;
  int n = tag (m);

  /** v: drop volume, ymax: drop max y, pc: mass center
   * ut: speed on trailing edge, ul: speed on leading edge
  */
  double v[n], ymax[n], um[n];
  int touchRight[n], touchLeft[n];
  double pcx[n], pcy[n], pcz[n];
  for (int j = 0; j < n; j++) {
    v[j] = pcx[j] = pcy[j] = pcz[j] = um[j] = 0.;
    ymax[j] = -HUGE;
    touchRight[j] = touchLeft[j] = 0;
  }

  // Determine some parameters on each droplet (volume, positions...)
  foreach ( reduction(+:v[:n]) reduction(max:ymax[:n]) reduction(+:pcx[:n]) reduction(+:pcy[:n]) reduction(+:um[:n]) reduction(max:touchRight[:n]) reduction(max:touchLeft[:n]) )
    if (m[] > 0) {
      int j = m[] - 1;
      if (cs[] > 1e-3) {
        v[j] += dv()*f[];
        um[j] += norm(u)*dv()*f[];
        ymax[j] = max(ymax[j], y);
        pcx[j] += dv()*f[]*x;
        pcy[j] += dv()*f[]*y;
        #if dimension == 3
        pcz[j] += dv()*f[]*z;
        #endif

        if (cs[] < 1.-1e-3) {
          if (x >= 1/2.) // touch the right solid
            touchRight[j] = 1;
          else if (x <= -1/2.) // touch the left solid
            touchLeft[j] = 1;
        }
      }
    }

  double vtmp = 0.;
  int jDrp = -1;
  int cloggingForAll = 0;
  for (int j = 0; j < n; j++) {
    // Identify the biggest but centered droplet or below -HoL (jDrp)
    if (v[j] > 1e-2) {
      cloggingForAll = max(cloggingForAll, touchLeft[j] && touchRight[j]);
      if ((sq(pcx[j]/v[j]) < 1.5*sq(1/2.) || pcy[j]/v[j] < -0.7*HoL || (touchLeft[j] && touchRight[j])) && v[j] > vtmp) {
        vtmp = v[j];
        jDrp = j;
      }
    }
  }

  if (jDrp != -1) {
    // Give mass center its true value
    double xc = pcx[jDrp]/v[jDrp], yc = pcy[jDrp]/v[jDrp];
    double umm = um[jDrp]/v[jDrp];

    // Clogging condition
    int clogging = cloggingForAll || (ymax[jDrp] >= -HoL);

    fprintf (fout, "%g %d %g %g %g %g\n",
      t, clogging, vtmp, umm, xc, yc);
    fflush (fout);

    // Stop simulation based on drop trailing edge position
    if (!clogging) {
      plot ("imgend.jpg");
      dump (file = "snapshotend.dump");
      exit(0);
    }
  }
}