/**
There is a missing metric in the update_tracer of the embed.h header 
which causes a wrong flux computation for axisymmetric + embed cases
using the bcg scheme for advection of tracers.

At line 828:
```
f[] += dt*(flux.x[] - flux.x[1])/(Delta);
```
should be 
```
f[] += dt*(flux.x[] - flux.x[1])/(Delta*cm[]);
```

This has an impact only in the domain (not at cut-cell boundaries)
and it only happens if embed is used in combination with axi 
(there might be issues with other metrics from other solvers too):
- If embed is used without any other metrics, it works fine because cm[] = 1.
- If axi is used alone, this function is not called and the correct update is used.

Here is a minimal pipe example taken from a [forum topic](https://groups.google.com/g/basilisk-fr/c/UtUOJfVlSXM/m/l-T47pDvAAAJ) reporting the problem.
Embed is just in the header, but not initializaed such that 
it can be compared to a full domain simulation without cut-cell.
pipecase = 0, 1 or 2 can be used to test the reference case, the source emb.h and the fixed emb.h

Comment: I noticed it afterward by compiling the code on the server, but it seems that dimensional analysis is not respected with the source embed.h while it is retrieved using either the modified embed_axi_tracers.h or deactivating the embedded boundaries. This is another hint that the metric coefficient should be there.
*/

#define pipecase 1

#if pipecase == 1
#include "embed.h"
#elif pipecase == 2
#include "boniouv/fixes/embed_axi_tracers.h"
#endif
#include "axi.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "view.h"

#if !EMBED
#include "fractions.h"
#endif

#define LEVEL 7

scalar f[];
scalar * tracers = {f};

double Reynolds = 1.;
face vector muv[];
int ivtk = 0.;

int main()
{
  L0 = 8. [1];
  init_grid (1 << (LEVEL));
  mu = muv;
  run(); 
}

double R = 0.999, U0 = 1.;

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*2.*R*U0/Reynolds;
}

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < R);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

#if EMBED
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
#endif

u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

event init (t = 0)
{
  fraction(f, x < 2. && y < R);
}

/** Output of fields*/
event png_output (t = 3.) {
  output_ppm (u.x, file = "velocity.png", n = 512);
  output_ppm (f, file = "tracer.png", n = 512);
}

/**
* With the embed.h, you see that the tracer $f$ is not correctly transported
* With the fixed embed_axi_tracers.h, the solution compares well with the solution without any of the embed headers (reference solution).

![Tracer is not correctly transported with embed.h](axi_embed_tracers/tracer.png)

*/


/**
However, there are still some points of confusion:

* This bug does not seem to affect velocity, even though it is calling this function at some point
  in the prediction step to compute the advective term for $u$.
* The small CFL fix might also contain bugs in the presence of a metric, as some inconsistencies seem
  to appear (sometimes `cm[]` is used, then `cs[]` is used for the update of neighbour cells), but I
  did not get far enough to call it a "bug", as I did not have issues from this fix.
*/