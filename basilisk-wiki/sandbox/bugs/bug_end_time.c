/**
This small example demonstrates the inability to change the *date* of the final event of the simulation during execution.
This happens with both `event end (t = end_time)` ("time" condition) and `event end (i = end_iteration)` ("iteration" condition) types of criteria.

In both cases, the simulation ends when reaching the first definition of the final date (`end_time` or `end_iteration`), regardless if we increased or decreased it.
*/

#include "run.h"

/**
# Time condition
*/
#if 1
/**
Let `end_time` be the initial date of the final event, which ends the run.
*/

double end_time = 10.;


int main ()
{  
  run();
}


event init (i = 0)
{
  printf ("Before init, end_time is %g\n", end_time);
/**
For some reason, let's decide to change the final time.
This can be done for instance if a convergence criterion hasn't been reached and the run needs to last a little longer.
*/
  end_time += 1.;
    
  printf ("After init, end_time is %g\n", end_time);
}


event end (t = end_time) {
  printf ("End event reached at t = %g, whereas end_time is %g\n", t, end_time);
}


/**
The results we get in `out` are as follows :\
`Before init, end_time is 10`\
`After init, end_time is 11`\
`End event reached at t = 10, whereas end_time is 11`\

`# Quadtree, 1 steps, 4.6e-05 CPU, 3.1e-05 real, 1.32e+08 points.step/s, 0 var`
*/


/**
# Iteration condition
*/
#else

/**
Now let's try the same but this time the run ends when the number of iterations `i` reaches `end_iteration`.
*/
int end_iteration = 10;

int main ()
{  
  run();
}

event init (i = 0) {
  printf ("Before init, end_iteration is %i\n", end_iteration);
  end_iteration -= 1;
  printf ("After init, end_iteration is %i\n", end_iteration);
}

event end (i = end_iteration) {
  printf ("End event reached at i = %i, whereas end_iteration is %i\n", i, end_iteration);
}

#endif

/**
We get similar results :\
`Before init, end_iteration is 10`\
`After init, end_iteration is 9`\
`End event reached at i = 10, whereas end_iteration is 9`\

`# Quadtree, 10 steps, 5.3e-05 CPU, 3.8e-05 real, 1.08e+09 points.step/s, 0 var`
*/
