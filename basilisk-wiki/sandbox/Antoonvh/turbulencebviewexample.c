/**
# My bview example
Click to see the result via youtube:

[![Click here to view the movie with youtube](https://img.youtube.com/vi/yV8Vf8iQ4Xo/0.jpg)](https://www.youtube.com/watch?v=yV8Vf8iQ4Xo)

*/

#include "grid/octree.h"
#include "utils.h"
#include "view.h"
vector u[];
scalar err[];
scalar c[];

int main(){
  init_grid (N);
  restore("dumpiso");
  X0=Y0=Z0=-L0/2.;// The domain in centered around {x,y,z}={0,0,0}, so that the theta and phi angles can be readily used to change the perspective. 
  scalar dis[];
  foreach(){
    foreach_dimension()
      dis[]+= dv()*(sq(u.x[1] - u.x[-1]) +
                    sq(u.x[0,1] - u.x[0,-1]) +
		    sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    if (fabs(x)>120 || fabs(y)>120 || fabs(z)>120)
      dis[]=0.;
  }
  while(adapt_wavelet({dis},(double[]){0.01},8).nc);
  coord d = {0,1,0};
  for (double i=1.;i<=512;i++){
    double p =(i+50)/250;
    double t = 0.5*(sin(i/50)); 
    view(fov=40,tx=-0.,ty=-0.,theta=t,phi=p,);
    double a = ((int)i % 256) - 128.;
    if (i<=255){
      squares("dis",linear=true,map=gray,min=0,max=0.05,alpha=a);
      cells(d,a);
    }
    if (i>=256){
      squares("dis",min=0,max=0.05,n=d,alpha=a);
      cells(alpha=a);
    }
    box();
    static FILE * fp = NULL;
    if (!fp) fp = popen ("ppm2mp4 movie.mp4", "w");
    save (fp = fp);
    clear();
    fprintf(ferr,"frame #%g\n",i); //Track the rendering process in the terminal
  }
}


/**
*/
