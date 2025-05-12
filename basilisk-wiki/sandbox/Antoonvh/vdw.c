/**
# A lumped-parameter bottom-boundary condition 
On this page the methodology and first results for a so-called Van-de-Wiel-Boundary condition are presented. 

We will do a two-dimensional Navier-Stokes simulation for a boussinesq fluid. Gravity is therefore modelled with a buoyancy field ($b$), that is advected and diffused.   
*/

#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
scalar b[];
scalar * tracers= {b};

/**
## The system and its parameters
The fluid is initialized with a constant stratification whos strength is specified by the [Brunt Vaisala frequency](https://en.wikipedia.org/wiki/Brunt%E2%80%93V%C3%A4is%C3%A4l%C3%A4_frequency) ($N^2$) related to the buoyancy field ($b$) according to,
   
   $$b=yN^2,$$
   
Where $y$ is the coordinate in the direction of the buoyancy-force vector, i.e. also known as the upward vertical direction. To drive a flow, a buoyancy flux forcing at the bottom boundary ($S$) is modelled with an horizontally homogeneous periodic-in-time function according to,

$$S=B_0 \mathrm{max}\left(\mathrm{sin}(\frac{2\pi t}{T}),0\right) \right)+B_1,$$
*/
#define Sbflx (max(B0*sin(2*M_PI*t/T),0)+B1)
/**
where $B_0$ and $B_1$ are buoyancy flux scales and $T$ is a timescale associated with the aforementioned periodicity. Furthermore, a *dynamic* linearized negative feedback ($G$) for the bottom buoyancy ($b_{bottom}$) is introduced, according to the following additional buoyancy flux term;

   $$G=\Lambda (b_{eq}-b_{bottom}),$$     
   
*/
#define Gbflx (Lambda*(beq-(b[]+((b[]-b[0,1])/2.))))
/**
where $b_{eq}$ is an equilibrium bottom buoyancy value associated with the $G$ term and $\Lambda$ is the lumped parameter that represents the coupling strength associated with the feed-back processes in the system. Furthermore, momentum and buoyancy are diffused with the fluid's diffusivity $\kappa$ (i.e. $\mathrm{Pr}=1$). 

In order to limit the physical system's degrees of freedom, we set $b_{eq}=0$, consistent with the initial stratification. With the remaining system parameters; $\{N^2, T, B_0, B_1, \Lambda\, \kappa\}$ we can identify (6 parameters - 2 dynamical units=) four indipendent dimensionless groups. There exist infinetly many sets of such groups, an example set of groups is as follows;

$$\Pi_1=\frac{B_0}{B_1}$$

$$\Pi_2=TN^2,$$  

$$\Pi_3=\frac{B_0T^2}{\kappa},$$ 


$$\Pi_4=\frac{B_0T}{\Lambda}.$$

Where $\Pi_{1,2,3,4}$ is the only group to include $B_1, N^2, \kappa$, $\Lambda$, respectively. The physical interpretation of these groups may become appearent to me in a discussion with collegues, whom will most probably propose new sets of groups. What is the change that I pick their favorite set out of the infinitly many possibilities? Anyway, we digress: For the numerics on this page it does not really matter and without further ado, we choose $\Pi=-\pi$, so that the time-integrated buoyancy flux scale $\int_{0}^{t_{end}}S(t^*)\mathrm{d}t^*$ will not divergence for $t_{end} \rightarrow \infty$. Furthermore me choose;

$$\Pi_2=30, \Pi_3=360000, \Pi_4=120.$$

Actually, the domain height and width are also important physical parameters that should be included in the analysis as well, and would warrant the introduction of two new dimensionless groups. However, we are only interested in domains with an aspect ratio of unity (thats $\Pi_5$=1!) and a width of $L0$, where $L0$ is large enough such that the presence of the boundaries, apart from the bottom, does not affect to relevant statistics of the flows.
*/
double beq= 0.0; // beq is not used
double B0 = 0.04;
double B1 = -0.04/M_PI;
double sqN=1.;
double T = 30.; 
double kappa = 0.0001;
double Lambda = 0.01;

b[bottom]=neumann(0.);// Explicitly omit diffusive fluxes for buoyancy.
u.t[bottom]=dirichlet(0.); // no-slip
b[top]=neumann(sqN);   // Keep the initialized stratification strength at the top boundary. 
// b[top]=dirichlet(sqN*(L0)); //<- A fixed top buoyancy, consistent with the initial stratification, could also be chosen. 
face vector av[];
const face vector muc[]={kappa,kappa};
double tend;
FILE * fp1;
FILE * fp2;
/**
## Set-up
We use a square box of size $3\times 3$, discretization starts with 128 grid points in both directions. We tell the solver that we have a non-zero viscosity (*muc*) and store the acceleration due to buoyancy in the previously declared *av* field. Furthermore, the simulation is run until the physical time $t_{end}=3.5T$.  
*/
int main(){
  L0=3.;
  Y0=X0=0.;
  init_grid(128);
  a=av;
  mu=muc;
  tend=3.5*T;
  run();
}
/**
   We have learned [here](internalwacesAMR.c) and [there](internalwacesAMR.c) that, when dealing with a linearly stratified fluid on tree grids, the consistency of the solver benefits from paying some special attention to the pressure field's attributes for prolongation and refinement. Also we have shown [else where](afluidatrest.c) that it is wise to start the simulation with small values for the tolerance on the poisson solver's maximum residual (*TOLERANCE*) and maximum timestep (*DT*) when using a volume force that modifies the pressure. The stratification is initialized, and a noisy pertubation is introduced to the dynamic fields to minic a random pertubation.    
*/
event init(t=0){
  p.prolongation = refine_linear; //3-rd-order accurate for grid-alligned variation
  p.refine = refine_linear;  //3-rd-order accurate for grid-alligned variation
  DT=0.0001;
  TOLERANCE = 10E-9;
  foreach(){
    b[]=(y*sqN)+0.005*noise();
    u.x[]=0.005*noise();
    u.y[]=0.005*noise();
  }
  boundary(all);
  fp1 = popen("gfsview-batch2D vdw.pv.gfv |ppm2mp4 > fstfield.mp4 ","w");
  fp2 = popen("gfsview-batch2D vdw.ps2.gfv |ppm2mp4 > sndfield.mp4 ","w");
}
/**
## Implementation
The effects of gravity can convieniently be introduced by an acceleration term that is calculated from the definition of the buoyancy field. 
*/
event acceleration(i++){
  coord del = {0,1.}; // y-dir unit vector
  foreach_face()
    av.x[]=del.x*(b[]+b[-1])/2.;
  boundary_flux({av});
}
/**
   The buoyancy field is diffused in the event below. Also we calculate and apply the tendency due to the buoyancy flux $B$ and the dynamical feed-back term $G$. In order to obtain the tendency from the flux, we realize that for a simulation in $D$ dimensions, the flux is distributed over a cell of size $\Delta^D$, whereas the buoyancy "fluxes" into the cell, via its bottom face of size $\Delta^{D-1}$. Hence, the tendency will be, $\beta = (S+G)/\Delta$.     
*/
event tracer_diffusion(i++){
  scalar beta[];
  foreach(){
    beta[]=0.; // This may include additional tendencies, e.g. due to large-scale advection. For now we set them to be zero.
    if(y<Delta) // That are the Bottom cells
      beta[]=((Sbflx+Gbflx)/Delta);
  }
  if (i%20==0){ //diagnosis of (the implementation of) the total buoyancy flux and its components, $S$ and $G$.
    static FILE * fpp = fopen("integrated_fluxes","w");
    double bf=0.,S=0,G=0;  
    foreach_boundary(bottom){
      bf+=beta[]*pow(Delta,(double)dimension);
      S+=Sbflx*pow(Delta,(double)(dimension-1));
      G+=Gbflx*pow(Delta,(double)(dimension-1));
    }
    fprintf(fpp,"%g\t%d\t%g\t%g\t%g\n",t,i,bf,S,G);
    fflush(fpp);
  }
  
  /**
Apply the calculated tendency and do the diffusion.
  */
  diffusion(b,dt,mu,r=beta);
}
  /**
     Once every $T/2$, the spatially variable S and G fluxes are outputted to a file together with the near surface buoyancy, total flux and local grid-cell size.   
  */
event spatial(t+=T/2.){
  char name[100];
  sprintf(name,"spatial%g",t);
  FILE * fp = fopen(name,"w");
  foreach_boundary(bottom)
    fprintf(fp,"%g\t%g\t%g\t%g\t%g\t%g\n",x,Sbflx,Gbflx,b[],(Sbflx+Gbflx),Delta);
  fclose(fp);
}
/**
   After 100 small initial timesteps a consistent pressure field should have been found by the solver, and hence we can relax back to more sensible values for *TOLERANCE* and *DT* so that the simulation will actually start to progress in the physical time.
*/
event defs(i=100;i++){
  CFL=0.8;
  TOLERANCE=min(TOLERANCE*1.1,10E-3);
  DT=min(DT*1.001,0.1); 
}
/**
Adaptation is based on the estimated discretization error in the representation of the velocity component fields and the buoyancy field, and allow the algorithm to refine upto 9 levels of refinement.
*/
event adapt(i=120;i++)
  adapt_wavelet((scalar *){u,b},(double[]){0.02,0.02,0.02},9);
/**
   For some visual reference we output a movie showing us the evolution of the $\| \nabla b \|$ field on a logaritmic colour scale and the buoyancy $b$ itself. 
*/

event movie(t+=T/200.;t<=tend){
  double bf=0.;
  foreach_boundary(bottom)
    bf+=-mu.y[]*(b[]-b[ghost]);
  scalar nb[],xp[];
  foreach(){
    xp[]=x/100.-(L0/200);
    nb[]=0.;
    foreach_dimension()
      nb[]+=fabs((b[1]-b[-1])/(2*Delta));
    nb[]=log(fabs(nb[])+0.00001);
  }
  boundary({nb});
  //static FILE * fptest = popen("gfsview2D vdw.pv.gfv","w");
  //output_gfs(fptest); //<- For configuring the gfv-file
  output_gfs(fp1);
  fprintf(fp1, "Save stdout {format = PPM width = 700 height = 400}\n");
  output_gfs(fp2);
  fprintf(fp2, "Save stdout {format = PPM width = 700 height = 400}\n");
  printf("%d %g\n",i,t);
}

/**
## Results
One may want to watch the movie for the $\mathrm{log}\left( \| \nabla b \| \right)$ field (top left) and the buoyancy field (top right and bottom).

<video width="700" height="400" controls>
<source src="vdw/fstfield.mp4" type="video/mp4">
</video>
</br>
<video width="700" height="400" controls>
<source src="vdw/sndfield.mp4" type="video/mp4">
</video>

That looks OK-ish. The fun part is in the fluxes: 

~~~gnuplot Fluxes over time
set xlabel 't'
set ylabel 'bottom integrated fluxes'
plot 'integrated_fluxes' u 1:3 w lines t 'Total' lw 3 ,\
      'integrated_fluxes' u 1:4 w lines t 'S' lw 3,\
      'integrated_fluxes' u 1:5 w lines t 'G' lw 3
~~~

The feed-back term $G$ seems to play an important role in the dynamics. Also notice the minimum in the total flux at the onset of each stable episode. Interestingly, because we have not prescribed the bottom buoyancy *nor* the flux(!), we can plot their spatial structure, e.g. at $t=1.5T$, the buoyancy looks like; 

~~~gnuplot Spatial structure of the near surface buoyancy
set xlabel 'x'
set ylabel 'Buoyancy'
plot 'spatial45' u 1:4 title 'Near surface buoyancy'
~~~ 

  and the fluxed look like;

~~~gnuplot Spatial structure of the buoyancy fluxes at t=45
set ylabel 'flux'
set key box 3 
set key center
plot 'spatial45' u 1:5 title 'Total flux' ,\
     'spatial45' u 1:2 title 'S' ,\
     'spatial45' u 1:3 title 'G'
~~~

It all appears to have worked well.
*/