
#include "grid/quadtree.h"
#include "vof-simple.h"
#include "run.h"

scalar f[], * interfaces = {f};
f[left]    = neumann(0);
f[right]   = neumann(0);
f[top]     = neumann(0);
f[bottom]  = neumann(0);
vector uf[];

int main (int argc, char * argv[])
{
  init_grid (64);
  L0 = 1.0;
  origin (-L0/2, -L0/2, 0);
  run();
}

event init (i = 0) {
  scalar dist[]; vertex scalar phi[];
  foreach() 
    dist[]=0.2-sqrt(sq(x)+sq(y));
  foreach_vertex()
    phi[]=0.25*(dist[]+dist[-1]+dist[0,-1]+dist[-1,-1]);
  
  fractions(phi,f);
  boundary ({f});

  FILE * fp = fopen("normal.dat", "w");
  foreach()
    if(f[]!=1&&f[]!=0){
       double r = sqrt(sq(x)+sq(y));
       double theta = (y>=0?acos(x/r)*180/pi:360-acos(x/r)*180/pi);

       // Normal vector calculated by Mixed-Youngs-Centerd scheme
       coord nmyc = mycs(point,f); 

       // Normal vector calculated by Youngs normal approximation
       coord nyng = youngs_normal (point, f);

       // Normal vector calculated from the distance function
       coord ndst;
       foreach_dimension()
         ndst.x = -(dist[1]-dist[-1])/(2*Delta);

       fprintf(fp,"%5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n",
               theta,nmyc.x,nmyc.y,nyng.x,nyng.y,ndst.x,ndst.y);
       }
   fclose(fp);
}

![Configuration of theta along the droplet interface ](configuration.png)

![X component of normal vector](nx.png)

![Y component of normal vector](ny.png)


