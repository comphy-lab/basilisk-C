/**
![A snapshot of a electro-magnetically generated dipole in a shallow fluid layer, visualzed by dye. Image courtesy of Chun Hou Wong.](/sandbox/Antoonvh/chundip.jpg)

# The Evolution of a point-vortex pair in a viscous fluid.
[C.H. Wong](http://repository.tue.nl/794399) studied the electro-magnetic generation of a dipolar vortex in a shallow fluid layer in a lab. Using PIV techniques, he found that the initial flow structure could accurately be discribed with a so-called point-vortex-dipole model. He also found that due to diffusion of vorticity, the structure evolved into a Lamb-Chapligin-type vortex via a stage that was best characterized by a so-called Super-Smooth-dipole model. In his analysis he focussed on the functional relation describing the vorticity and streamfunction inside the dipole's atmosphere, where the stream function is required to be evaluated in the co-moving frame. Such a function relation exists for all coherent vortex structures. Furthermore, the excentricity of the dipole's atmosphere (or separatrix) itself was studied. This closed curve encloses the part of the fluid that is entrained by the dipolar vortex structure. On this page we aim to see if we can reproduce his results using a 2-dimensional (2D) flow solver. 

## Set-up
In this simulation we omit the electromagnetic forcing phase of the experiment. Rather, we initialize the flow with two small so-called Rankine vortices, that form a point-vortex-like dipole. 
*/

#include "navier-stokes/centered.h"
#include "fractions.h"
/**
Two oppositly-signed point vortices with strength $\Gamma$ and separation distance $d$ advect trough a fluid with veloicity $U=\frac{\Gamma}{2\pi d}$. We set the strength of our initial vortices such that we can use normalized distance and velocicity scales. 
*/
double yp=0.5;
double xp=0.;
double Gamma = 2*M_PI; 
double rs = 0.02;
int maxlevel = 12;
/**
We define some usefull global variables. 
*/
double an=0.1;
int j=0;
double xd[3];
double sf[9999],vort[9999],ypp[9999]; //arrays for streamfunction, vorticity and y-coordinate.
scalar psi[],omega[];
FILE * fp;
char name[100];
/**
## Boundary conditions
We solve the flow in a frame that moves with $U$ in the x-direction. This way the solver can use larger (CFL-limited) time steps and the vortex structure does not cross the domain boundary. This requires us to set a inflow velocity at the boundaries.   
*/
u.n[left]=dirichlet(-1);
u.n[right]=dirichlet(-1);
/**
Since we will using and evaluating stream functions ($\psi$) to discribe our incompressible flow field, we need to set consistent boundary conditions according to the inflow velocity. Also, we want our $\psi$ field to enherit the anti-symetric nature of the dipole vortex.  
*/
psi[left]=dirichlet(y);
psi[right]=dirichlet(y);
psi[top]=dirichlet(X0+L0); //Free slip
psi[bottom]=dirichlet(X0); //Free slip

/**
## Initialization
We also define a function that refines the grid and initializes the vorticity field according to the Rankine vortex model. 
*/
void init_Rankine(double xp,double yp, double gamma, double rs, int ml){
  refine(sq(x-xp)+sq(y-yp)<sq(rs*10) && level <= ml-2);
  refine(sq(x-xp)+sq(y-yp)<sq(rs*5) && level <= ml-1);
  refine(sq(x-xp)+sq(y-yp)<sq(rs*3) && level <= ml);
  scalar f[];
  fraction(f,sq(x-xp)+sq(y-yp)<sq(rs));
  foreach()
    omega[]+=gamma/(M_PI*sq(rs)) * f[];
}

int main(){
  L0=20;
  X0=Y0=-L0/2;
  init_grid(512);
  run();
}
/**
We use a three-stage initialization process where the vorticity field is initialized first, then the corresponding stream function is evaluated by solving the associated Poisson problem. Finally the velocity components can be found by using the definition of the stream function.   

Also, since diffusion of vorticity, due to the fluids viscousity ($\nu$), is imporant for the physics of this study, we can identify a Reynolds number ($Re$) with system parameters {$U,d,\nu$}, according to; 

$$Re = \frac {Ud}{\nu}.$$

We choose the parameters such that $Re = 2000$. 
*/
event init(t=0){
  const face vector muc[]={0.005,0.005};
  mu=muc;
  init_Rankine(xp,yp,-Gamma,rs,maxlevel);
  init_Rankine(xp,-yp,Gamma,rs,maxlevel);
  boundary({omega});
  poisson(psi,omega);
  boundary({psi});
  foreach() {
    u.x[] = (psi[0,1]-psi[0,-1])/(2*Delta);
    u.y[] = -(psi[1,0]-psi[-1,0])/(2*Delta);
  }
  boundary(all);
}
/**
## Output
We output .mp4 movies to provide some visual reference.
*/
event movies(t+=0.05;t<=20){
  scalar lev[];
  foreach() {
    omega[]=((u.y[1,0]-u.y[-1,0]) - (u.x[0,1]-u.x[0,-1]))/(2*Delta);
    lev[]=level;
  }
  boundary({omega});
  static FILE * fp1 = popen ("ppm2mp4 g.mp4 ", "w");
  output_ppm (omega, fp1,n=pow(2,10),min=-10,max=10,linear = true);
  static FILE * fp2 = popen ("ppm2mp4 l.mp4 ", "w");
  output_ppm (lev, fp2,n=pow(2,10),min=0,max=maxlevel);
}
/**
Here is the result for the vorticity field. 

<video width="512" height="512" controls>
<source src="wong2015/g.mp4" type="video/mp4">
</video> 

We see that the dipole only starts to slow down when the positive and the negative vorticity starts 'eating' each other. Also it does not help the propagation speed that circulation starts to leak out of the dipole's atmosphere.

## The excentricity and the $\omega - \psi$ relation
If one can find a function $f$ such that, $\omega = f(\psi) $, the vortex strucute does not deform due to advection. We only have a chance of finding such a function in the frame of reference that co-moves with the dipole. Therefore we need to find the velocity of the dipole in the solver's frame of reference so that we can perform the required transformation of $\psi$.     
*/
event output(t=0.4;t+=an){
  xp=0;
  double xm=0;
  double ens=0,momg=0;
  /**
  We find the x-coordinate of the location of maximum vorticity (*xm*) and the 'centre of mass' of the enstrophy (*xd*). 
  */
  foreach(){
    omega[]=((u.y[1,0]-u.y[-1,0])- (u.x[0,1]-u.x[0,-1]) )/(2*Delta);
    xp+=x*sq(omega[])*sq(Delta);
    ens+=sq(omega[])*sq(Delta);
    if (fabs(omega[])>momg){
      momg=fabs(omega[]);
      xm=x;
    }
  }
  xd[j]=xp/ens;
  /**
  Store $\omega$ and $\psi$ for next (t+=0.1) output event. We will modify the stream function later, once we have found a second-order accurate approximation of the dipoles translation velocity. 
*/
  if (j==1){ 
    int g = 0;
    while (fabs(vort[g])>0.){
      vort[g]=0.;
      sf[g]=0.;
      ypp[g++]=0.;
    }
    boundary({omega});
    poisson(psi,omega);
    g=0;
    foreach(){
      if (fabs(omega[])>(momg/100.) && g<9999){
        vort[g]=omega[];
        sf[g]=psi[];
        ypp[g++]=y;
      }
    }
    /**
    For an estimation of the excentricity ($\epsilon$) of the diople's atmosphere, a first-order accurate approximation of the velocity has to suffice.
    */
    double vel =  ((xd[1]-xd[0])/0.1);
    foreach()
      psi[]+=vel*y;
    /**
    For some additional visual reference, we output an animated .gif of the $\psi=0$ contour, evaluated in the co-moving frame, and the vorticity field. 
    */
    static FILE * sfp = popen ("gfsview-batch2D wong2015.psi.gfv |ppm2gif --delay 200>gg.gif","w");
    output_gfs(sfp);
    fprintf (sfp, "Save stdout { format = PPM width = 500 height = 500}\n");
  /**
  We can observe an ellipsoid separatrix that due to the case set-up corresponds to $\psi=0$. 
  
  ![Vorticity field and the $\psi \approx 0$ contour for the 19 snapshots of our analysis.](wong2015/gg.gif)
  
  */
    double hy = 0.5;
    double ps1=-1.;
    double ps2 = -1.;
    /** 
    We find the half-length of the long(?) axis by looking for the y-coordinate where $\psi$ changes sign (for $y>0.5).
    */
    while (ps1*ps2>=0.){ 
      ps1=ps2;
      ps2=interpolate(psi,xm,hy);
      hy+=0.001;
    }
    hy-=0.001;
    /**
    The half-length of the short(?) axis is found by locating the x-coordinate of the stagnation point, that should be on the separatrix.   
    */
    ps1=1.;
    ps2=1.;
    double hx = 0;
    while (ps1*ps2>=0.){
      ps1=ps2;
      ps2=interpolate(u.x,hx+xm,0.)-vel;
      hx+=0.001;
    }
    hx-=0.001;
    /**
    And we write the corresponding excentricity to a file.
    */
    static FILE * fps = fopen("excent","w");
    fprintf(fps,"%g\t%g\t%g\n",t,hy/hx,hy);
  }
  /**
    Here is the result:
    
    ~~~gnuplot From a statistics' perspective,the excentricity is decreasing over time. 
    set xlabel 'Time'
    set ylabel 'Excentricity'
    plot 'excent' using 1:2
    ~~~
    
    Notice that $\epsilon=1.21$ corresponds to the value for the atmosphere of a point-vortex dipole.
    
    Next, the $\omega$ and corresponding $\psi$ values inside the dipole are written to a file. 
    */
  if (j==2){ //Print the stored omega and Psi in co-moving frame
    double vel = ((xd[2]-xd[0])/0.2); // relative velocity of co-moving frame.
    fprintf(ferr,"%g\t%g\t\n",t-0.1,vel+1);
    int g = 0;
    sprintf(name,"op%g",t-0.1);
    fp = fopen(name,"w");
    while(fabs(vort[g])>0){
      fprintf(fp,"%g\t%g\t%g\t%g\n",sf[g]+(vel*ypp[g]),vort[g],sq(sf[g]+(vel*ypp[g])),sq(vort[g]));
      g++;
    }
    fclose(fp);
    an=0.8;
    j=0;
  }else{// Wait for t+=0.1;
    an=0.1;
    j++;
  }
}
/**
Here is the result:

~~~gnuplot t=1.5 and t=3.5 
set xlabel 'Psi'
set ylabel 'Omega'
plot 'op1.5' using 1:2 title 't=1.5' ,\
'op3.5' using 1:2 title 't=3.5' 
~~~

Remarkably, there does indeed appear to be a one-to-one relation between $\psi$ and $\omega$. The initial dipole appears to be characterized by a high-order-uneven-power polynominal function for $f$, (i.e. as in $\omega = f(\psi)$). 

~~~gnuplot t=10.5 and t=15.5 
set xlabel 'Psi'
set ylabel 'Omega'
plot 'op10.5' using 1:2 title 't=10.5' ,\
'op15.5' using 1:2 title 't=15.5' 
~~~

For $t>10$, some of the vorticity has diffused outside the separatrix. Inside the dipole, there appears to be a linear relation between $\omega$ and $\psi$, a feature that characterized the Lamb-Chaplygin dipole model ([That tends to collide with solid objects on its trajectory](lamb-dipole.c)).


## Adaptation?

yes, grid adaptation was used:
    
    */
event adapt(i++){
  adapt_wavelet((scalar*){u},(double[]){0.05,0.05},maxlevel);
}
/**
Video of the evolution of the level of refinement:

<video width="512" height="512" controls>
<source src="wong2015/l.mp4" type="video/mp4">
</video> 

The maximum resolution was only applied in a small fraction of the spatial and temporal(!) domain. So it appears to be quite a good idea to use the adaptive-grid approach for this problem.   

## Further reading
The MSc. thesis of C.H. Wong provides a complete analysis of 2D vortex structures as observed with high-resolution PIV in a lab. The analysis was made possible by the unprecedented fidelity of novel particle-seeding techniques and dipole-generation mechanisms developed and explored by C.H. Wong.

The [report](http://repository.tue.nl/794399) is made available via the TU/e repository. It not only served as inspiration for this page, it also is a great starting point for further reading on this topic. 
*/