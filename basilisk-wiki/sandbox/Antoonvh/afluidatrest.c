/**
# A fluid at rest
As a simple test case, a fluid at rest is simulated. The physical system consists of a square box that is filled with a fluid with a density $\rho$ that is therefore affected by the acceleration due to gravity (*gr*).

## Set-up
The Navier-Stokes *flow* solver is used.
*/
#include "navier-stokes/centered.h"

face vector av[];
double uscale = 10E-8;
int j,ii;
double gr = 9.81;

int main(){
  Y0=0.;
  L0=10.;
  a=av;
  /**
  The test is run for five different cases of initializations.
  */
  for (j=0;j<=4;j++)
    run();
}

event init(t=0){
  init_grid(1<<6);
  /**
  ## Case 0
  We naively start the simulation with a limited timestep.
  */
  
  DT=1.;
  /**
  ## Case 1
  Here we do not want to rely on the poisson solver to find the pressure field in the first time step, therefore we help it by initializing the hydrostatic pressure field and setting consistent boundary conditions.
  */
  if (j==1){
  p[bottom]=dirichlet(L0*gr);
  p[top]=dirichlet(0.);
   foreach()
     p[]=gr*(L0-y);
   boundary({p});
 }
   /**
## Case 2
Same as case 1 but now the domain is set-up such that the domain-averaged pressure is zero. This way, the initial guess of a zero pressure field at the root cell is exactly correct. 
*/ 
  if (j==2){
    Y0=-L0/2.; // Not really needed?
    p[bottom]=dirichlet(L0*gr/2.);
    p[top]=dirichlet(-L0*gr/2.);
    foreach()
      p[]=-gr*y;
    boundary({p});
  }
  
   /**
   ## Case 3 & 4
   We now again do rely on the Poisson solver to find a correct pressure field for initialization, but compared to the naive case 0, we decrease the tolerated residual of the iterative approach, and decide to help the solver a bit by decreasing the maximum timestep. We only do this for a few timesteps.
*/
  if  (j>2){
    Y0=0.; // Default
   p[top]    = neumann(a.n[ghost]*fm.n[ghost]/alpha.n[ghost]); // Default
   p[bottom] = neumann(-a.n[]*fm.n[]/alpha.n[]); // Default
   DT = 0.1;
   TOLERANCE=10E-6;
   }
} 
  
event acceleration (i++){
  /**
  We aply a body force that will modify the pressure such that $\partial p / \partial y = -gr$. 
  */
  ii++;
  coord del={0,-gr};
  foreach_face()
    av.x[]=del.x;
  /**
  For case 3 and 4, the timestep *DT* and *TOLERANCE* are set back to their default values after 10 integration timesteps. For case 4 we also do the same a few timesteps after restoring.
  */
  if (j>2 && (i==10||ii==10)){
    DT=1.;
    TOLERANCE=10E-3;
  }
}   
/**
## Dump and restore
Apart from different initialization, we also check how the methods fare when dump en restore are used. Thereupon, the simulation is dumped and restored at $t=50$, i.e. halfway the total simulation physical time. 
*/
event dump_and_restore(t=50){
  char name[100];
  sprintf(name,"dumpert%d",j);
  dump(name);
  restore(name);
  DT=1.;
  /**
  ## Case 4
  For the fourth case, the tolerance and timestep are again decreased. After ten timesteps these are set back to their default values again in the acceleration event. 
  */
  if (j==4){
    ii=0;
    TOLERANCE=10E-6;
    DT=0.1;
  }
  dtnext(DT);
}
/**
## Output
The total kinetic energy is logged and we make movies of the evolution of the *u.x* velocity fields for each case. 
*/

event t_is_hundred(t=100){}
event resul(i+=1;t<=100.){
  double e=0.;
  foreach(reduction(+:e))
    e+=(sq(u.x[])+sq(u.y[]))*sq(Delta);
  if (j==0){
    char name[100];
    sprintf(name,"%d.dat",j);
    static FILE * fpp = fopen(name,"w");
    fprintf(fpp,"%g\t%g\n",t,e+10E-15);
    static FILE * fp = popen ("ppm2mp4 jj.mp4", "w");
    output_ppm (u.x, fp,min=-uscale,max=uscale, n = 256, linear = false );
  }
  if (j==1){
    char name[100];
    sprintf(name,"%d.dat",j);
    static FILE * fpp1 = fopen(name,"w");
    fprintf(fpp1,"%g\t%g\n",t,e+10E-15);
    static FILE * fp1 = popen ("ppm2mp4  j1.mp4", "w");
    output_ppm (u.x, fp1,min=-uscale,max=uscale, n = 256, linear = false );
  }
  if (j==2){
    char name[100];
    sprintf(name,"%d.dat",j);
    static FILE * fpp2 = fopen(name,"w");
    fprintf(fpp2,"%g\t%g\n",t,e+10E-15);
    static FILE * fp2 = popen ("ppm2mp4  j2.mp4", "w");
    output_ppm (u.x, fp2,min=-uscale,max=uscale, n = 256, linear = false );
  }
  if (j==3){
    char name[100];
    sprintf(name,"%d.dat",j);
    static FILE * fpp3 = fopen(name,"w");
    fprintf(fpp3,"%g\t%g\n",t,e+10E-15);
    static FILE * fp3 = popen ("ppm2mp4 j3.mp4", "w");
    output_ppm (u.x, fp3,min=-uscale,max=uscale, n = 256, linear = false );
  }
  if (j==4){
    char name[100];
    sprintf(name,"%d.dat",j);
    static FILE * fpp4 = fopen(name,"w");
    fprintf(fpp4,"%g\t%g\n",t,e+10E-15);
    static FILE * fp4 = popen ("ppm2mp4 j4.mp4", "w");
    output_ppm (u.x, fp4,min=-uscale,max=uscale, n = 256, linear = false );
  }
}

/**
## Results:
From the movies, we can observe that there are spurious currents within the domain. Only case 1 and 2 do not have them after initialization, but do so after the *dump and restore* event.  

Case 0:

<video width="256" height="256" controls>
<source src="afluidatrest/jj.mp4" type="video/mp4">
</video>


Case 1 and 2, left and right, respectively:

<video width="256" height="256" controls>
<source src="afluidatrest/j1.mp4" type="video/mp4">
</video> 
<video width="256" height="256" controls>
<source src="afluidatrest/j2.mp4" type="video/mp4">
</video> 

case 3 and 4, left and right, respectively:

<video width="256" height="256" controls>
<source src="afluidatrest/j3.mp4" type="video/mp4">
</video> 
<video width="256" height="256" controls>
<source src="afluidatrest/j4.mp4" type="video/mp4">
</video> 


Lets have a more quantified look at the kinetic energy:

~~~gnuplot
set yr [0.000000000000005: 0.0000001]
set logscale y 
set format y "10^{%L}"
set xlabel 'Time'
set ylabel 'Energy + 10^{-14}'
plot "0.dat" using 1:2 title "Case 0" with lines lw 3 ,\
     "1.dat" using 1:2 title "Case 1"  with lines lw 3 ,\
     "2.dat" using 1:2 title "Case 2" with lines lw 3,\
     "3.dat" using 1:2 title "Case 3" with lines lw 3 ,\
     "4.dat" using 1:2 title "Case 4" with lines lw 3
~~~

Use strategy 4 to mitigate this issue. Never use strategy 1 or 2!

## Bug or Feature?
I think this page shows that there is an error associated with finding an accurate pressure field during a single timestep when starting from a poor initial guess. However, the severity of this issue depends on the TOLERANCE, grid-resolution (not shown) and timestepping parameters (now shown) and *gr* in this case. Therefore, it remains just the effect of finite (but controllable) discretization errors in the simulation, and therefore is not a bug. However, I would *guess* that dumping and restoring the pressure field as well might solve this issue when using *restore*.  
*/
