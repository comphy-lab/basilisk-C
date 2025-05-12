/** 
# Testing `initial_conditions_2Dto3D.h`

In this example use the results from a 2D two-phase simulation to initialize a
3D grid. 2D results were saved inside a file
[snapshot_t6.50000.h5](convert/snapshot_t6.50000.h5) format using
[output_xdmf.h](../output_fields/xdmf/output_xdmf.h.html). This file contains:
```
  >>h5ls -r snapshot_t6.50000.h5
  /                        Group
  /Cells                   Group
  /Cells/f                 Dataset {2335438, 1}
  /Cells/p                 Dataset {2335438, 1}
  /Cells/u.x               Dataset {2335438, 3}
  /Geometry                Group
  /Geometry/Points         Dataset {2762309, 3}
  /Topology                Dataset {2335438, 4}
```
Results correpond to a square domain with $L_x=L_x=0.24$ using 4096x4096
points. We may load results from a region of size $L'_x=L'_x=0.03$ using
$512x512$ points. The file is converted into a gnuplot-compatible format (in
double precision) using
[convert_snapshots_to_binary.py](convert_snapshots_to_binary.py) as follows:
```
  >> python convert_snapshots_to_binary.py snapshot_t6.50000.h5 10 0.03
```
which should create a set of binary files (`snapshot_t6.50000_l.bin`,
`snapshot_t6.50000_f.bin`, and so on.)

An example of the 2D volume fraction field is shown on the left, and the
corresponding 3D `vof` interface is shown on the right:

<center><img src="initial_conditions_2D.png" alt="drawing" width="400"/>
<img src="initial_conditions_3D.png" alt="drawing" width="400"/>
<figcaption>An example of the intial condition in 2D (left) 
and the corresponding 3D interface (right). </figcaption>
</center>
<br/><br/>

*/


#include "view.h"

#define MAXLEVEL 9  
#define _mindel (L0/N)
#define D0 (L0/8.)
#if dimension == 3  
  #define rectanglebox(extra)                                     \  
    intersection((D0 / 2. + extra - y), (D0 / 2. + extra + y))
#endif
const char *file_restart_path = "../convert/snapshot_t6.50000";
#include "initial_conditions_2Dto3D.h"

int main()
{
  L0 = 1.00;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << MAXLEVEL;
  init_grid(N);

  scalar f[], p[];
  vector u[];
  foreach(){
    f[] = 0;
    p[] = 0;
    foreach_dimension()
      u.x[] = 0;
  }
  
  initial_condition_2Dto3D(f, u, p, D0/2., -D0/2.);
  restriction({f, u});
  
  #if dimension == 3
  /** If 3D, we set to zero the fields outside the region of interest  */
  foreach(){
    f[]*= (rectanglebox(1.*_mindel) > 0.);
    p[]*= (rectanglebox(1.*_mindel) > 0.);
    foreach_dimension(){
      u.x[]*= (rectanglebox(1.*_mindel) > 0.);
    }
  }
  #endif

  /** Then, visualize the results just to make sure  */

  scalar l[];
  foreach()
    l[] = level;

  #if dimension == 2
  {
    draw_vof("f");
    squares("f", linear = false);
    save("init_f.png");

    box();
    squares("u.x", linear = false);
    save("init_u.png");

    box();
    squares("u.y", linear = false);
    save("init_v.png");

    box();
    squares("l", linear = false);
    save("init_grid.png");
  }
  #endif

  #if dimension == 3
  {
    view(camera="bottom");
    box();
    draw_vof("f");
    save("init_f.png");

    view(camera="bottom");
    box();
    squares("u.x", linear = false, n = {1, 0, 0}, alpha = 0.);
    squares("u.x", linear = false, n = {0, 1, 0}, alpha = 0.);
    squares("u.x", linear = false, n = {0, 0, 1}, alpha = 0.);
    save("init_u.png");

    view(camera="bottom");
    box();
    squares("u.z", linear = false, n = {1, 0, 0}, alpha = 0.);
    squares("u.z", linear = false, n = {0, 1, 0}, alpha = 0.);
    squares("u.z", linear = false, n = {0, 0, 1}, alpha = 0.);
    save("init_w.png");

    view(camera="iso");
    box();
    squares("l", linear = false, n = {1, 0, 0}, alpha = 0.);
    squares("l", linear = false, n = {0, 1, 0}, alpha = 0.);
    squares("l", linear = false, n = {0, 0, 1}, alpha = 0.);
    save("init_grid.png");
  }
  #endif
  
}