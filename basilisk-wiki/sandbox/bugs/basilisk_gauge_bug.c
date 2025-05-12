/**
There were missing brackets. */

#include "saint-venant.h"
#define MAXLEVEL 5

int main()
{
  init_grid (1 << MAXLEVEL);

  run();

}

event init (i=0)
{
  foreach(){
    zb[]=-1;
    h[]=max(0.,0.1*x-zb[]);
  }
  boundary ({zb,h});
}

Gauge gauges1[] = {
  // file   lon      lat         description
  {"centre1", 0.5,0.5, "centre"},
  {NULL}
};

event outgauges1 (i++;i<=10) {
  output_gauges(gauges1,{eta,u.x,u.y});
}

Gauge gauges2[] = {
  // file   lon      lat         description
  {"centre2", 0.5,0.5, "centre"},
  {NULL}
};

event outgauges2 (i++;i<=10) {
  output_gauges(gauges2,{eta,u.x,u.y});
}
