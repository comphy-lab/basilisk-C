/**
Lets see how we can generate and display a mp4 movie in the sandbox using GFS view. I try to learn from this [page](http://www.basilisk.fr/sandbox/popinet/step.c) in Stephane's sandbox. Therefore we need to generate some frames:
*/
#include "navier-stokes/centered.h"
#include "view.h".

scalar omg[],r[],s[],psi[];
int maxlevel = 8;
double xo= -2.5, yo= 0;
double k=3.83170597;
double temp= 5;

int main(){
  L0=15;
  X0=Y0=-L0/2;
  init_grid(1<<(maxlevel-1));
  run();
}

event init(t=0){
  refine(pow(pow((x-xo),2)+(pow((y-yo),2)),0.5) < 1.5 && level<maxlevel);
  foreach() {
    r[] = pow(pow((x-xo),2)+(pow((y-yo),2)),0.5);
    s[] = (y-yo)/r[];
    psi[] = ((r[]>1)*((1/r[]))*s[]) + ((r[]<1)*((-2*j1(k*r[])*s[]/(k*j0(k)))+(s[]*r[])));
  }
  boundary(all);
  foreach() {
    u.x[] = ((psi[0,1]-psi[0,-1])/(2*Delta));
    u.y[] = -(psi[1,0]-psi[-1,0])/(2*Delta);
  }
  boundary(all);
}

event adapt(i++){
  adapt_wavelet((scalar *){u},(double []){0.05,0.05},maxlevel,maxlevel-4); 
}

event movie(t+=0.05;t<=temp){
  foreach() {
    omg[]=(u.x[0,1]-u.x[0,-1] - (u.y[1,0]-u.y[-1,0]))/(2*Delta);
  }
  boundary({omg});
   /**
   Using GFS view for a MP4 movie
   */
  static FILE * fpmp4 =
    popen ("gfsview-batch2D mp4movie.outzoom-vor.gfv | "
           "ppm2mp4  lamb.mp4", "w");
  output_gfs (fpmp4);
  fprintf (fpmp4, "Save stdout { format = PPM width = 512 height = 512}\n");
  /**
  Or the ouputppm function for an animated gif
  */
  static FILE * fpgif = popen ( "ppm2gif > lambgif.gif", "w");
  output_ppm (omg,fpgif,n=512,min=-10,max=10,linear=true);
  /**
  And with Bview we also make an MP4.
  */
  view(width=512,height=512);
  clear();
  squares("omg",min=-10,max=10,linear=true);
  cells();
  save("bviewed.mp4");
    
}

/**
Display the movie as animated gif:

![This works fine](mp4movie/lambgif.gif)

Now we try to display the movie in the .mp4 format rendered with gfsview:

![The GFSview movie](mp4movie/lamb.mp4)

and with Bview,  

![The Bview movie](mp4movie/bviewed.mp4)

The gif is about 1.5 Mb, whereas the MP4's are about 120 Kb. 
*/