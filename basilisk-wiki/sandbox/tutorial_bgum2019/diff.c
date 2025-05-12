/**
# Tutorial: Diffusion of a Gaussian pulse

![Challenge: Code a diffusion problem in Basilisk](diff/s.mp4)

Here a simple diffusion problem is solved. One could have found
[here](/src/README) that the [`diffusion.h` solver](/src/diffusion.h)
exists. Further, from an [example](/src/examples/brusselator.c), we
learn that we need to include the time loop ourselves, unlike for the
popular [Navier-Stokes solver](/src/navier-stokes/centered.h), which
already includes it for the user.
*/

#include "diffusion.h"
#include "run.h"

 /**
The evolution of a scalar field `s` is considered. It is decleared
like so:
 */
scalar s[];

/**
The time loop is started by calling `run()`. The value of `DT` is set,
which will serve as a time-stepping parameter. It became available by
`#include`-ing `run.h`.
 */
int main() {
  DT = 0.05;
  run();
}
/**
The solution is initialized at `t = 0`. Apart from taking care of the
upcoming `events`, the `run()` command also initializes a grid. By
default Basilisk runs in two dimensions on a quadtree grid, and uses
`N = 64` cells in both directions. Further, the default domain size
and origin are; `L0 = 1., {X0,Y0} = {0.,0.}`). The `foreach()` loop
iterates over all cells and sets the cell-centered coordinates (`x,y`)
in the background.
 */
event init (t = 0) {
  foreach() 
    s[] = exp(-(sq(x - 0.5) + sq(y - 0.5))*10.); 
}
/**
A `.mp4` movie is generated that displays the evolution of `s`, using
the `output_ppm()` function. With the relavant converters
[installed](/src/INSTALL), we only need to set the `file`
name. Further, the resolution and the colorbar range are set.
 */
event mov (t += 0.1)  
  output_ppm (s, file = "s.mp4", n = 256, min = -1, max = 1);

/**
   In the following event the time integration is performed. Again,
   inspired by a relevant [example](/src/examples/brusselator.c), we
   tell the time-loop to set the actual timestep (`dt`), based on our
   maximum value (`DT`). The documentation for the [diffusion
   solver](/src/diffusion.h) reveals the proper sequence of the
   arguments to the `duffusion()` function (albeit via a
   `struct`ure). A constant diffusivity field (`kap`) is declared and
   initialized. The values are defined on cell faces. The faces have
   distinct directions and hence we need to set the corresponding
   vector components in both dimensions.
 */
event diff (i++) {
  dt = dtnext (DT);
  const face vector kap[] = {0.01, 0.01};
  diffusion (s, dt, kap);
}
/**
   Finally, a `data` file is writen to check if the scalar field `s`
   is conserved.
 */

event lot (i += 5) {
  static FILE * fp = fopen ("data", "w");
  fprintf (fp, "%g %d %.8g\n", t, i, statsf(s).sum);
}
/**
~~~gnuplot We can plot the results with `gnuplot`
set xlabel 'time'
set ylabel 'total s'
set grid
plot 'data' u 1:3
~~~

The scalar (`s`) is not exactly conserved due to the finite
`TOLERANCE` for the iterative solver. Note that additional errors are
present in the solution field due to the discrete representation of
the solution (i.e. due to the limited resolution in space and time),
and due to the binary representation of floating point numbers.

The time loop stops when it "sees" no further events. Therefore, we
request it to continue to go to the `stop` event at `t = 10`. This
event could have many other names.
 */
event stop (t = 10);
