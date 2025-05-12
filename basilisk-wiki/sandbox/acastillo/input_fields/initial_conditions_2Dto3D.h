/**
# Using 2D results to initialize a 3D simulation

<center><img src="initial_conditions_2D.png" alt="drawing" width="400"/>
<img src="initial_conditions_3D.png" alt="drawing" width="400"/>
<figcaption>An example of the intial condition in 2D (left) 
and the corresponding 3D interface (right). </figcaption>
</center>
<br/><br/>

This function reads 2D simulation results from a binary file in a format
compatible with the gnuplot binary matrix format in double precision, see
[auxiliar_input.h](auxiliar_input.h). The 2D results are then used to initialize
a 3D simulation, which may be useful to reduce computational cost by avoiding
long transients, or to focus on the development of 3D instabilities from a 2D
base state.

### *initial_condition_2Dto3D()*: Function to read 2D simulation results.

This function reads a series of files containing the grid level
(`example_l.bin`), volume fraction (`example_f.bin`), in-plane velocities
(`example_u.bin` and `example_v.bin`), and pressure (`example_p.bin`), and
stores them in the corresponding fields. If used in 2D, data is loaded in the
$(x,y)$-plane. If used in 3D, data is loaded in the $(x,z)$-plane, and the
velocity $v_y$ is set to zero. If there are walls in the transversal direction, 
we may load the fields only between $y\in[d_1,d_2]$. All files should start 
with the same prefix stored inside a variable `file_restart_path` 

The arguments and their default values are:

*f*
: scalar field to be initialized.

*u*
: vector field to be initialized.

*p*
: scalar field to be initialized.

*d1*
: upper bound in the $y$-coordinate

*d2*
: lower bound in the $y$-coordinate

#### Example Usage

```c
  #define D0 (L0/2.)
  scalar f[], p[];
  vector u[];
  const char *file_restart_path = "example";
  #include "initial_conditions_2Dto3D.h"
  initial_condition_2Dto3D(f, u, p, D0/2., -D0/2.);
```

see, also [example 1](test_init_2Dto3D.c)

*/


#include "auxiliar_input.h"
// Define some useful macros
#define cond1(y,d1,d2) ((l[] < level) && ((y < d1) && (y > d2)))
#define cond2(y,d1,d2,del) ((y > d1+del) || (y < d2-del))
#define wallbox(d, extra) intersection((d + extra - y), (-d + extra + y))
void initial_condition_2Dto3D(scalar f, vector u, scalar p, double d1, double d2){
  // Initialize transveral velocity 
	foreach()
	  u.y[] = 0.;
  
  // Now, we match the refinement level
  scalar l[];
  read_matrix(file_restart_path, "_l", l);      

#if dimension == 3
  // If 3D, we may unrefine outside of the region of interest
  int maxlevel = MAXLEVEL;
  for (int li = maxlevel; li >= 4; li--){
    unrefine( (cond1(y,d1,d2) || cond2(y,d1,d2,16*_mindel)) && level > li);
  }
#endif

  // Then, we read the corresponding fields
  read_matrix(file_restart_path, "_f", f);
  read_matrix(file_restart_path, "_u", u.x);
  #if dimension == 2
    read_matrix(file_restart_path, "_v", u.y);
  #else
    read_matrix(file_restart_path, "_v", u.z);
  #endif
  read_matrix(file_restart_path, "_p", p);
}
// Delete the macros
#undef wallbox
#undef cond1
#undef cond2