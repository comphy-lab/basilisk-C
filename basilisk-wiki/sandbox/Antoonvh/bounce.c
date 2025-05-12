/**
![Astronaut Scott Kelly played ping pong in space. Video courtesy of [NASA Johnson](https://www.youtube.com/watch?v=TLbhrMCM4_0) ($\leftarrow$ click to watch the full movie on youtube).](https://media.giphy.com/media/y2lIYf0EmMlSU/giphy.mp4)

# A Bouncing Droplet
While doing some online 'research' [here](https://9gag.com/gag/a3KMwPN), it was decided to investigate the rebound of a water droplet from a collision with a hydrophobic wall in the absence of gravity. We identify that the most prominent physical features of this system can be solved for with the tools that are provided by Basilisk. We therefore call upon the most popular Navier-Stokes solver, the two-phase extension module and furthermore it seems wise to include the effects of surface tension as well. The curiousity is triggered as i wanted to check if it is easy enough to implement a hydrophobic boundary that gives somewhat consistent results. And I was curious to see if the movement of the hands in the movie is required for a rebound.     
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
/**
## The case set-up
First the boundary conditions are set; a no-slip and non-stick wall (i.e. the contact angle associated with the water-air-wall interface is $180^o$).  
 */
f[bottom]=0;
u.t[bottom]=dirichlet(0.);
/**
The physical system consists of a liquid droplet with a radius $R$ that travels trough the air towards the hydrophobic wall with a velocity $U$. Gauged from the movie, we approximate $R \approx 2\ \text{ cm}$ and $U\approx 3\ \text{ cm/s}$. Together with the fluid properties of the air and water ($\rho_a,\rho_w,\mu_a,\mu_w, \sigma$) we can identify some dimensionless groups that describe the system. 

 $$ Re = \frac{\rho_a UR}{\mu_a} \approx 30, $$ 

 $$ We = \frac{\rho_w U^2R}{\sigma} \approx 0.2, $$

 $$ \Pi_1 = \frac{\mu_w }{\mu_a} \approx 100, $$

 $$ \Pi_1 = \frac{\rho_w }{\rho_a} \approx 1000 \rightarrow 100, $$

 $$ D = 3 \rightarrow 2, $$

where $D$ stands here for the number of spatial dimensions of the system and the '$\rightarrow$' symbols indicate dimensionless groups where a consession is made with respect to the reality onboard the space station. This is done in order to keep the numerical experiment feasable on a single-core system. It would be relatively straight forward to increase these numbers to the more realistic values. 

We set the model parameter such that we obtain the approximated ratios, aproximately. It is done in such a way that we have a normalized velocity scale ($U$) and length scale ($R$). Furthermore, we set the maximum grid resolution to correspond to a $512^2$ equidistant grid.
 */

double R = 1.;
double U = 1.;
int maxlevel = 9;

int main(){
  mu2=1./30.; //gas phase
  mu1=100./30.; //liquid phase
  rho2=1.; //gas phase
  rho1=100;//liquid phase
  f.sigma=500.;
  init_grid(1<<7);
  L0=10;
  X0=Y0=-L0/2;
  run();
}

/**
## initialization
After we have refined the grid with a ring of high resolution mesh, the droplet is initialized an it is targeted at the bottom wall. 
 */

event init(i=0){
  refine(sq(x)+sq(y)<sq(R+0.5) && sq(x)+sq(y)>sq(R-0.25) && level<maxlevel);
  fraction(f,sq(R)-sq(x)-sq(y));
  foreach()
    u.y[]=-f[]*U;
}
/**
## Adaptation
Since the advective interface tracking and resolving for the surface tencile waves only requires a high resolution mesh at the locations of the interface, it seems smart to use grid adaptivity, we will check it this was indeed the case.  
*/
event adapt(i++)
  adapt_wavelet((scalar *){u,f},(double []){0.02,0.02,0.001},maxlevel);

/**
## Output
The general dynamics are visualized in a movie that shows the rebound of the droplet.
 */
event gfsviewq (t+=0.05;t<=10){
static FILE * fp5 = popen("gfsview-batch2D bounce.justphase.gfv |ppm2mp4 justphase.mp4","w");
output_gfs(fp5);
fprintf (fp5, "Save stdout { format = PPM width = 800 height = 450}\n");
}
/**

<video width="800" height="450" controls>
<source src="bounce/justphase.mp4" type="video/mp4">
</video> 

   Furtheremore, 'slow motion' movies are rendered where the velocities inside the bubble are visualized. This is supplemented with a movie that dislays the evolution of the grid. Finally, we output a file that logs the reduction factor, quantifying what fraction of grid cells we have used compared to a fixed and equidistant grid.
*/
event gfsview(t+=0.02;t<=10){
scalar omg[],ux[],uy[];
int n = 0;
foreach(){
    n++;
    ux[]=u.x[]*(f[]==1);
    uy[]=u.y[]*(f[]==1);
    omg[]=((u.x[0,1]-u.x[0,-1]) - (u.y[1,0]-u.y[-1,0]))/(2*Delta);
  }
  boundary({omg,ux,uy});
  static FILE * fp2 = popen("gfsview-batch2D bounce.arrows.gfv |ppm2mp4 arrows.mp4","w");
  output_gfs(fp2);
  fprintf (fp2, "Save stdout { format = PPM width = 800 height = 450}\n");

  static FILE * fp3 = popen("gfsview-batch2D bounce.cellsandvort.gfv |ppm2mp4 cv.mp4","w");
  output_gfs(fp3);
  fprintf (fp3, "Save stdout { format = PPM width = 800 height = 450}\n");

  static FILE * fpo = fopen("Pi","w");
  fprintf(fpo,"%g\t%d\t%d\t%d\n",t,i,n,(1<<(maxlevel*dimension))/n);
  printf("%g\t%d\t%d\t%d\n",t,i,n,(1<<(maxlevel*dimension))/n);
}
/**
You can view the results:

<video width="800" height="450" controls>
<source src="bounce/arrows.mp4" type="video/mp4">
</video> 
</hr>
<video width="800" height="450" controls>
<source src="bounce/cv.mp4" type="video/mp4">
</video> 

~~~gnuplot
set xlabel 'Iteration'
set ylabel 'Pi'
set key off
plot 'Pi' u 2:4 w l lw 3
~~~

## The dynamics from an energy perspective
We also want to understand the bounce event from an energy perspective. The droplet seems to be nearly stationary at some point. On may ask: *How can it than find the energy to regain momentum and travel bouncy back up again?*. Therefore, we log the kinetic energy of the flow. Furthermore, at the risk of spoiling the suspence of the reader, we define a helper function that helps determine the total circumference of the 2D droplet.
**/

void find_facets(Point point,scalar f,double xy[4]){
  coord n;
  n = mycs (point, f);
  double alpha = plane_alpha (f[], n);
  coord segment[2];
  if (facets (n, alpha, segment) == 2){
    xy[0] = x + segment[0].x*Delta;
    xy[1] = y + segment[0].y*Delta;
    xy[2] = x + segment[1].x*Delta;
    xy[3] = y + segment[1].y*Delta;
  }else{
    printf("Warning:\nCould not find facets; expect unexpected behaviour.\n");
  }
}
/**
This is used to help calculate the potential energy associted with the so-called "line energy" ($E_l$) with units $Jm^{-1}$, that is the 2D/planar analogy of [surface energy](https://en.wikipedia.org/wiki/Surface_energy). $E_l$ is defined as,

$$E_l = l \sigma,$$

where $l$ is the total length associated with the water-air interface. A Potential energy $P$ is defined as the difference between $E_l$ and its minimum geometrically possible value for a simply connected shape in 2D with surface area $A$ according to,

$$P = E_l - 2 \sigma \sqrt{\pi A} = \sigma (l-l_{min}),$$

where $l_{min}$ is the minimum circumference of a 2D shape with area $A$, i.e. based on a circle. 
*/
event energy(i+=5){
  static FILE * fpe = fopen("energy","w");
  double e=0;
  double l=0;
  double xyf[4];
  foreach(reduction(+:e) reduction(+:l)){
    e+=0.5*sq(Delta)*(sq(u.x[])+sq(u.y[])) * ((rho1*f[]) + (rho2*(1-f[])));
    if (f[]>0.00001 && f[] < 0.9999){
      find_facets(point,f,xyf);
      l+=pow(sq(xyf[0]-xyf[2])+sq(xyf[1]-xyf[3]),0.5);
    }
  }
  fprintf(fpe,"%g\t%d\t%g\t%g\t%g\t%g\n",t,i,e,(l-2*M_PI*R)*f.sigma,l,((l-2*M_PI*R)*f.sigma) + e);
}

/**
We can now view the evolution of the kinetic, potential and total energy. 

~~~gnuplot
set xlabel 'time U/R []'
set ylabel 'Energy '
set key on
plot 'energy' u 1:3 w l lw 3 title 'Kinetic',\
       'energy' u 1:4 w l lw 3 title 'Line Potential',\
       'energy' u 1:6 w l lw 3 title 'Total'
~~~

Corresponding to our visually-based suspicion, we can see the footprint of the collision where the kinetic energy is greatly reduced at $t\approx5 R/U$. After that, the kinetic energy literally bounces back up, to reach about 70% of the original value. Appearently the system has found a way to store its kinetic energy during the bounce. The plot of the evolution of the energy associated with the bubble interface and the total energy tells us that during the collision, the circumference of the droplet has increased, and thereby stores energy. As the droplet comes to a hold, the surface tension tries to minimize the droplet circumference again and thereby pushes the droplet away from the wall. Over the course of the simulation, energy is lost due to viscous effects. Also energy appears to be created in the initial approach, this is a numerical effect that may be reduced by employing a higher numerical precision or use even smarter numerical formulations (e.g. by adapting a [momemtum conserving advection scheme](http://basilisk.fr/src/navier-stokes/conserving.h)). 

## The next step
It seems attractive to study the dynamics in an axis symeteric limit such that we may get some preliminary insights on the dynamics of the full 3D system. The results of that study are presented [here](axibounce.c).   
*/
