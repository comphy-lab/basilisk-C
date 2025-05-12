/**
# An atmospheric boundary-layer filled with a smoke cloud.
We simulate a boudary layer that is filled with a so-called smoke cloud. This cloud cools at the top due to the emission of longwave radiation. The case was introduced as an intercomparison case by Bretherton et al. (1999). This page shows the implementation of the case in Basilisk using an LES formulation. However, on this page we run the case in 2 dimensions to that we can view the results. 

## Set-up
In order to invoke the LES formulation we include the Navier-stokes solver and the *SGS* headerfile. 
*/

//#include "grid/octree.h" // <- uncomment for 3D as was used in Van Hooft et al. (2018)
#include "navier-stokes/centered.h"
#include "SGS.h"
/**
The maximum and minimum levels of refinment are set. We also update the radiation tendencies once every 5 timesteps. This way we limit the effect of the fact that the current implementation of the tendency term due to radiation is not very efficient and should be improved in the future. Also we set some other global variables that according to the case description by Bretherton et al. (1999).   
*/
int maxlevel=8,minlevel=5;
int Radupdate = 5;
double T0 = 291.5, F0 = 60., p0=1.1436, Ka=0.02, Cp=1004.;

double temp = 675*sq(dimension+1)/(4-dimension); // ~1 hour in 2D or 3 hour in 3D
/**
We declare fields for the smoke cloud fraction (*S*), the (virtual) potential temperature (*T*), the radiative flux field (*F*), and the field that stores the temperature tendency (*co*). We also tell the solver that we want to advect *S* and *T* with the flow and diffuse it according to the sub-grid scale  model. Finally, an helper acceleration field *av* is declared that helps with the implementation of the 'gravity'.
*/
scalar S[],T[],F[],co[];
scalar * tracers = {T,S}; // Temperature and smoke fraction get both advected and diffused according to the SGS model.
face vector av[];

int main(){
  F0=F0/sq((4-dimension));
  L0 = 3200;
  a=av; 
  X0=Y0=Z0=0.;
  init_grid(1<<(minlevel));
  co.refine=refine_linear;
  run();
}
/**
## Initialization
The fields *T* and *S* field are initialized with an inversion at $z=700m$ (i.e. $y=700$ in our simulation) with a width of 25 metres. 
*/
event init (t = 0){
  DT=0.01;
  TOLERANCE=1e-8; 
  for (scalar s in tracers)
    s.gradient=minmod2;
  foreach(){
    S[] = 0.5+0.5*tanh((700-y)/12.5);
    T[] = T0+3.5*tanh((y-700)/12.5)+(0.1*noise()*(y<700));
    T[] += 0.001*(y-712.5)*(y>712.5);  
  }
  boundary(all);
  while(adapt_wavelet({S},(double[]){0.05}, maxlevel,minlevel).nf>10){
    foreach(){
      S[] = 0.5+0.5*tanh((700-y)/12.5);
      T[] = T0+3.5*tanh((y-700)/12.5)+(0.01*noise()*(y<700));
      T[] += 0.001*(y-712.5)*(y>712.5);  
    }
    boundary(all);
  }
}  
/**
Gravity is implemented via a buoyancy formulation based on the virtual potential temperature.
*/
event acceleration(i++){
  foreach_face(y) 
    av.y[]=9.81*((T[]+T[0,-1]-T0))/(2*T0);
}
/**
## Radiation
Bretherton et al. (1999) discribed how the radiation tendencies should be evaluated in a staggered-grid manner. We neglect that and do it in our own way so that we do not need to take proper care of face fluxes etc. 

Based on the smoke cloud fraction field ($S$) we can evaluate the readiative fluxes ($F$), according to:

$$F = F_0 e^{-K_a S_l},$$

where $F_0$ and $K_a$ are specified constants, and $S_l$ is defined for a point at height $z$ as follows:

$$S_l = \int_z^{z_{top}} \rho\ S\ dz,$$

where $\rho$ is taken as the reference density $\rho_0$, and with our grid orientation and set-up; $z_{top}=Y0+L0=L0$. 

The tendencies due to the radiation are applied each timestep. However the tendency field *co*, is only updated each *Radupdate* (i.e. 5) timesteps.
*/
event radiation(i++){
  if (i%Radupdate==0){
    boundary({S});
    foreach(){  //Calculate for each grid cell, there is room to improve this!
      double sl = 0;
      double xp = x;
      double yp = y;
      double zp = z;
      while (yp<L0){ // This is the not smart part of evaluating the radiation scheme. 
        Point point = locate(xp,yp,zp); 
        yp=y;// jump to the found cell centered height.
        sl+= p0*interpolate(S,xp,yp,zp)*Delta;
        yp+= Delta/1.5;
      }
      F[] = F0*exp(-Ka*sl);
    }
    boundary({F});
    /**
    We have found the flux, now we take the vertical derivative (1D radiative divergence) to obtain the tendency.
    */
    foreach(){
      co[] =-(F[0,1]-F[0,-1])/(2*p0*Cp*Delta);
    }
  }
  /**
  We copy the tendency into a helper field (*cooling*). We do this because the *cooling* field will be modified by the solver such that we cannot use it in case we do not update it next time step nor in the output, see [here](www.basilisk.fr/src/diffusion.h).
  */
 
  boundary({co});
  scalar cooling[]; 
  foreach()
    cooling[]=co[];
  boundary({cooling});
  diffusion(T,dt,zerof,beta = cooling);// Do cooling, not diffusion!
}
/**
##Adaptation
We use the well-known wavelet-based technique for error estimation and perform the corresonding grid adaptation.
*/
event adapt(i++){
  double eS = 0.1/(dimension-1);
  double eT = 0.7/(dimension-1);    
  adapt_wavelet({S,T},(double[]){eS,eT}, maxlevel,minlevel,all);
}
/**
This *back_to_defaults* event is added so that we can start with a very small timestep (*DT*) and tolerance for the Poisson problems (*TOLERANCE*), so that the iterative solver has some low-impact iterations before it finds the not-consistently-initialized pressure field. 
*/
event back_to_defaults(i=50){
  DT=1.;
  CFL = 0.7;
  TOLERANCE=1e-3; 
}
/**
## 2D output
In 2D we output three movies. 
*/
event mov(t+=10;t<=temp){
  if (dimension==2){
    double mc = 4./(3600*sq(4-dimension)); //Cooling scale in K/sec
    scalar lev[];
    foreach()
      lev[]=level;
    boundary({S,T}); //Linear interp for movie
    static FILE * fpg = popen("ppm2mp4 S.mp4","w"); 
    output_ppm(S,fpg,linear=true,n=pow(2,maxlevel),min=0,max=1,map=gray);
    
    static FILE * fpc = popen("ppm2mp4 c.mp4","w"); 
    output_ppm(co,fpc,linear=true,n=pow(2,maxlevel),min=-mc,max=mc);
    
    static FILE * fpl = popen("ppm2mp4 l.mp4","w"); 
    output_ppm(lev,fpl,n=pow(2,maxlevel),min=minlevel,max=maxlevel);
  }
  else if (dimension==3){
    //no movie, no nothing
    
  }
  fprintf(ferr,"%d %g\n",i,t);
}

/**
## results

We can view the movies. First we choose to show the smoke-cloud fraction field. 

<video width="512" height="512" controls>
<source src="smoke/S.mp4" type="video/mp4">
</video> 

We can observe the emergence of subsiding shells in the boundary layer. Next we have a look at the process that drives the flow; the radiative cooling:

<video width="512" height="512" controls>
<source src="smoke/c.mp4" type="video/mp4">
</video> 

Looks OK. Notice that the depth of the cooling layer is controlled by the "optical thickness" parameter *Ka* that was tuned to help the participating models by supressing the scale separation. Finally we look at the evolution of the grid's level of refinement:

<video width="512" height="512" controls>
<source src="smoke/l.mp4" type="video/mp4">
</video> 

Well done *adapt_wavelet* function!

## Reference:
  
Bretherton, Christopher S., et al. *"An intercomparison of radiatively driven entrainment and turbulence in a smoke cloud, as simulated by different numerical models."* Quarterly Journal of the Royal Meteorological Society 125.554 (1999): 391-423.

Read,print,copy and obtain it via [this link](https://pdfs.semanticscholar.org/55d5/7ce0eb7d89a2d87a47948a717af82f076c8a.pdf).
*/
