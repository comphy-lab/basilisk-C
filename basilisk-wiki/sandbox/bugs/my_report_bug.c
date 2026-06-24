#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "lambda2.h"
#include "navier-stokes/perfs.h"
#include "maxruntime.h"
const double MU = 0.01;
face vector av[];
double RUNID;            // For restarting purpose                                 
int maxlevel = 7;

int main (int argc, char * argv[])
{
  maxruntime (&argc, argv);
  if (argc > 1)
    RUNID = atoi(argv[1]);

  L0 = 2.*pi [1];
  dimensions(nx=4,ny=1,nz=1);
  foreach_dimension()
    periodic (right);

  mu[] = {MU,MU,MU};
  a = av;
  N = 2048;
  run();
}

event init (i = 0) {
  double u0 = 1., k = 1.;
  if (!restore (file = "restart"))
    foreach() {
      u.x[] = u0*(cos(k*y) + sin(k*z));
      u.y[] = u0*(sin(k*x) + cos(k*z));
      u.z[] = u0*(cos(k*x) + sin(k*y));
    }
}

event end(i = 5 + (RUNID * 10)){
  char name[80];
  sprintf (name, "restart");
  dump (file = name);
  return 1;
}

~~~~~~~
## multigrid3D dump/restart grid depth do not match ## 

When using multigrid3D header, at a specific amount of processors, dump/restore produces the following error: 

```text
grid depths do not match. Aborting.
```

For example, this can happened when,

- dimensions(1,1,1):  using 64 CPUs on $512 \times 512 \times 512$ grid point.
- dimensions(2,1,1):  using 1024 CPUs on $1024 \times 512 \times 512$ grid point.
- dimensions(4,1,1):  using 256 CPUs on $1024 \times 256 \times 256$ grid point.
- dimensions(4,1,1):  using 256 CPUs on $2048 \times 512 \times 512$ grid point.
  
This error started appearing after I updated to the 2026 version. I updated because I was trying to fix the issue described here:(https://basilisk.fr/sandbox/bugs/dump_integer_overflow.c)




 
            
