/**
# Clogging of a dual-plate channel by a droplet

We aim to model the descent of a droplet, initially motionless, between two 
parallel plates with parameters: a mean separation length $L$, height $H$, and 
corner radius $r$. The droplet, initially positioned at $(0, D/2)$, consists of
water in an airflow. We manipulate three key factors: the Bond number 
$\mathrm{Bo} = \Delta\rho g D^2 / \sigma$, where $\Delta\rho = \rho_1 - 
\rho_2$ with $\rho_1$ and $\rho_2$ denoting water and air densities 
respectively, $g$ representing gravity's acceleration, and $\sigma$ indicating 
surface tension. We also vary the aspect ratio $L/D$ and the contact angle 
$\theta$.

![](../fig/channel_clogging_boundary.png
){ style="display: block; margin: auto;" }
<p style="text-align: center;">
  Illustration of the initial state
</p>

We solve the Navier-Stokes equations for a two-phase flow, employing the 
momentum-conserving approach. The channel is constructed using an embedded 
boundary technique. Instead of employing the arithmetic definition of 
viscosity, we adopt a harmonic formulation.
*/
#include "embed.h"
#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "navier-stokes/conserving.h"

/**
We need functionalities related to surface tension, gravity, and contact angle 
(refer to contact-embed from Popinet) as well as visualization functions.
*/
#include "tension.h"
#include "reduced.h"
#include "contact-embed_popinet.h"
#include "view.h"

/**
We require the *tag()* function to precisely track the droplet's movement, 
especially if it undergoes breakup. *maxruntime* is employed to manage the 
maximum runtime on a supercomputer. We assess performance using *perfs*.
*/
#include "tag.h"
#include "maxruntime.h"
#include "navier-stokes/perfs.h"
#include "mygeometry.h"

#define fnorm(v) (sqrt(sq(v.x) + sq(v.y)))

/**
We define some parameters.
*/
// Simulation parameters:
double  tend      = 20.,            // End time
        uemax     = 0.1;            // Error on u field

int     iniLvl    = 4,              // Initial level of refinement
        iniDrpLvl = 5,              // Initial LOR for droplet
        iniGrdLvl = 6,              // Initial LOR for grid
        cellsInL  = 10,             // Min. cell number in channel width L
        maxLvl;

char    *restFile = "restore.dump", // Restore file name
        *mode     = "w";            // Writing mode: "w" or "a"

// Flow parameters:
double  rhoRatio  = 829.15,         // Density ratio rho1/rho2
        muRatio   = 63.82,          // Dynamic viscosity ratio mu1/mu2
        Mo        = 1.94e-12,       // Morton, Mo = g*mu1^4/(rho1*sig^3)
        theta0    = 65,             // Contact angle (°)
        uair      = -0.1,           // Inflow velocity
        Bo        = 5;              // Bond, Bo = rho1*g*D^2/sig

// Geometric parameters:
double  HoL       = 6.38,           // Aspect ratio H/L
        roL       = 0.3,            // Aspect ratio r/L
        DoL       = 2.,             // Aspect ratio D/L
        WoL       = 2.;             // Aspect ratio W/L


/**
The channel is defined as the combinaison of two curved rectangles.
*/
double channel (coord p) {
  coord p1 = {-1/2. - WoL,  -HoL};
  coord p2 = {-1/2.,        0.};
  coord p3 = {1/2.,         -HoL};
  coord p4 = {1/2. + WoL,   0.};
  return -1*union(curvedRectangle(p, p1, p2, roL),
      curvedRectangle(p, p3, p4, roL));
}

/**
The *generateDroplet()* function generates a droplet by first coarsely refining
the area containing the droplet. Subsequently, it generates a coarse droplet. 
Then, it refines only the interface and defines a new, more refined droplet. 
This process of refining the interface and generating the droplet is repeated 
iteratively until the maximum level is reached.
*/
void generateDroplet (coord p1, double D) {
  // Refine the area around the droplet
  for (int lev = iniLvl+1; lev<=iniDrpLvl; lev++)
    refine (sphere((coord) {x, y, z}, p1, max(Delta, D/2.)) > 0.
        && level < lev);

  // Refine the area at the interface of the droplet
  scalar drpField[];
  drpField.refine = drpField.prolongation = fraction_refine;
  fraction (drpField, sphere((coord) {x, y, z}, p1, D/2.));

  for (int lev=iniDrpLvl+1; lev<=maxLvl; lev++) {
    refine ((drpField[] > 0. && drpField[] < 1.) && level < lev);
    fraction (drpField, sphere((coord) {x, y, z}, p1, D/2.));
  }

  foreach ()
    if (drpField[] > 0.)
      f[] = min(f[]+drpField[], 1.);
}

/**
The *generateSolid()* function performs a similar process to generateDroplet(),
but for the solid.
*/
void generateSolid () {
  // Refine the area around the solid
  for (int lev=iniLvl+1; lev<=iniGrdLvl; lev++)
    refine (channel((coord) {x, y, z}) < 0. && level < lev);

  // Refine the area at the interface of the solid
  vertex scalar phi[];
  foreach_vertex()
    phi[] = channel((coord) {x, y, z});
  boundary ({phi});
  fractions (phi, cs, fs);

  for (int lev=iniGrdLvl+1; lev<=maxLvl; lev++) {
    refine ((cs[] > 0. && cs[] < 1.) && level < lev);

    vertex scalar phi[];
    foreach_vertex()
      phi[] = channel((coord) {x, y, z});
    boundary ({phi});
    fractions (phi, cs, fs);
  }      
}

/**
*plot()* generates an image depicting the current state of the channel and 
saves it under the designated name.
*/
void plot (char name[]) {
  scalar f2[], unorm[];
  foreach () {
    f2[] = (f[] > 1e-2 && cs[] > 0.) ? f[] : 0.;
    unorm[] = norm(u);
  }

  translate (y = -Y0 - L0/2.) {
    draw_vof (c = "cs", filled = -1, 
        fc = {0.522, 0.341, 0.341});
    draw_vof (c = "f2", filled = -1, fc = {1,1,1}, lc = {1,1,1});
    draw_vof (c = "f2", lw = 1);
    squares (color = "unorm", max = 10, min = 0.);   
    box ();
  }
  
  save (name);
}

/**
We have enforced a minor inflow at the top boundary. The left, right, and 
bottom boundaries adhere to outflow conditions (Neumann), which transition into
a wall condition (Dirichlet) in case of backflow.
*/
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

  /**
  We utilize command line arguments to adjust default parameters such as the 
  aspect ratio $L/D$, the Bond number $\mathrm{Bo}$, the contact angle 
  $\theta$, and the minimum number of cells in the gap $L$. The final parameter
  specifies a snapshot name, enabling the simulation to restart.
  */
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

  /**
  We adjust the dimension of `DT` to adhere to dimensional analysis principles.
  Additionally, we compute the length of `L0` to ensure it is sufficiently wide
  to accommodate the largest droplet.
  */
  DT = HUGE [0];
  WoL = max(2.*DoL, 4.);    // Solid thickness in mm
  L0 = 1.21[0]*max(2*WoL + 1, DoL + HoL + 2*DoL);
  X0 = -L0/2.;
  Y0 = -(HoL + 2*DoL + 0.1*L0);
  maxLvl = (int) ceil( log(cellsInL*L0)/log(2.) );
  
  init_grid (1 << iniLvl);
  
  /**
  Surface tension is determined based on the Bond number, while viscosity is 
  calculated using both the Bond and Morton numbers.
  */
  rho1 = 1.;
  rho2 = rho1/rhoRatio;
  G.y = -1.;
  f.sigma = (rho1 - rho2)*fnorm(G)*sq(DoL)/Bo;
  double La = sqrt(Bo/Mo);    // Laplace number
  mu1 = sqrt(f.sigma*rho1*DoL/La);
  mu2 = mu1/muRatio;

  run();
}

/**
## Initialisation

In case of a simulation failure, we attempt to restore it by regenerating the 
channel and then the droplet. If restoration is successful, we regenerate the 
channel because information on the `fs` field is absent in the snapshot (for 
now).
*/
event init(t = 0.) {
  // Unsuccessful restoration
  if (!restore (file = restFile)) {
    generateSolid ();

    coord pDrp = {0., DoL/2.};
    generateDroplet(pDrp, DoL);
  }
  // Successful restauration
  else {
    generateSolid ();
  }
}

/**
Mesh adaptation is activated, refining the mesh in both interfacial regions 
(fluid and solid) as well as in areas of high velocity.
*/
event adapt (i++) {
  scalar f1[], cs1[];
  foreach () {
    f1[] = f[];
    cs1[] = cs[];
  }	

  adapt_wavelet ((scalar*){f1, cs1, u}, 
      (double[]){0.001, 0.001, uemax, uemax}, maxLvl);
}


/**
## Ouput

Movies and snapshots (dump files) are generated, serving purposes such as 
illustrating, debugging, restoring a previous simulation.
*/
event plotEvent (t += 1e-1) {
  plot ("channel_clogging_movie.mp4");
}

event snapshot (t += 1) {
  char name[100];
  sprintf (name, "%07.3f.dump", t);
  dump (file = name);
}

/**
Output is generated and stored in the `out` file. The simulation continues 
until it reaches a time limit specified by `tend`.
*/
event output (t += 1e-1, t <= tend) {
  /**
  Droplets are enumerated using the *tag()* function. Any connected region 
  for which f[] > 1e-3 (i.e. a droplet) will be identified by a unique “tag” 
  value between 0 and n-1.
  */
  scalar m[];
  foreach ()
    m[] = f[] > 1e-3;
  int n = tag (m);

  /** 
  With each cell now tagged with a unique number, we proceed to calculate 
  various data regarding all droplets, including their volume and position.
  */
  double  v[n],             // Droplet volume
          ymax[n],          // Maximum ordinate of the droplet
          pcx[n],           // Position of mass center (on x)
          pcy[n];           // Position of mass center (on y)
  bool    touchSolid[n];    // Boolean, true if the droplet touch the solid
  for (int j = 0; j < n; j++) {
    v[j] = pcx[j] = pcy[j] = 0.;
    ymax[j] = -HUGE;
    touchSolid[j] = false;
  }

  foreach ( reduction(+:v[:n]) reduction(max:ymax[:n]) reduction(+:pcx[:n]) 
      reduction(+:pcy[:n]) reduction(max:touchSolid[:n]) )
    if (m[] > 0) {
      int j = m[] - 1;
      if (cs[] > 1e-3) {
        v[j] += dv()*f[];
        ymax[j] = max(ymax[j], y);
        pcx[j] += dv()*f[]*x;
        pcy[j] += dv()*f[]*y;

        if (cs[] < 1.-1e-3 && x >= 1/2. && x <= -1/2.)
          touchSolid[j] = true;
      }
    }

  /**
  We now aim to identify the droplet we wish to track. This droplet is 
  determined as the largest droplet, which is either centered or located below 
  $-H/L$.
  */
  double vtmp = 0.;
  int jDrp = -1;
  bool cloggingForAllDrp = false;
  for (int j = 0; j < n; j++) {
    if (v[j] > 1e-2) {
      cloggingForAllDrp = max(cloggingForAllDrp, touchSolid[j]);
      if ((sq(pcx[j]/v[j]) < 1.5*sq(1/2.) || pcy[j]/v[j] < -0.7*HoL || 
            touchSolid[j]) && v[j] > vtmp) {
        vtmp = v[j];
        jDrp = j;
      }
    }
  }

  /**
  In conclusion, we output the volume and position of the identified droplet, 
  along with the state of the channel (whether clogged or free).
  */
  if (jDrp != -1) {
    double xc = pcx[jDrp]/v[jDrp], yc = pcy[jDrp]/v[jDrp];
    bool notClogged = !cloggingForAllDrp && ymax[jDrp] < -HoL; 
    fprintf (fout, "%g %d %g %g %g\n",
      t, notClogged, vtmp, xc, yc);
    fflush (fout);

    /**
    The simulation halts early if the channel is not clogged and if the droplet
    we are tracking has descended below $-H/L$.
    */
    if ( notClogged ) {
      plot ("channel_clogging_movie.mp4");
      char name[100];
      sprintf (name, "%07.3f.dump", t);
      dump (file = name);
      exit(0);
    }
  }
}

/**
## Results

![Animation of droplet sliding into a channel.](channel_clogging_movie.mp4)(
width="800"height="600")

![](../fig/channel_clogging_results.png
){ style="display: block; margin: auto;"}
<p style="text-align: center;">
Simulations of the penetration of a water droplet inside two-parallel plates
for an aspect ratio $L/D = 0.5$. In cases (a) and (b), a contact angle of 
$65°$ is employed, while cases (c) and (d) involve a contact angle of 
$150°$. The Bond number is specified as $0.278$ for cases (a) and 
(c), and $2.154$ for cases (b) and (d). The colors represent the 
magnitude of the velocity field.
</p>

# References
*/