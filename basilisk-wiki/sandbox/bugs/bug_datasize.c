/**
The bug we address here is one encountered when (unfortunately) declaring a variable named `datasize`.
Here is an small example.
*/


/**
Let's just include some simple 1D headers.
*/
#include "grid/multigrid1D.h"
#include "run.h"

/**
We only put 4 points but this bug appears for more values of `N`, and in 2D as well.
*/

int main()
{
  N = 4;
  run();  
}


/**
We setup a simple `init` event where we declare a variable named `datasize` as an integer and a scalar named `toto`.
*/

event init (i = 0) {
  int datasize = 7;
  scalar toto[];
  /**
  We fill `toto` with random values and display them in the same loop.
  */
  printf ("Value of toto in the first loop :\n");
  foreach() {
    toto[] = noise();
    printf("%g ", toto[]);
  }

  /**
  We display `toto` again, in another loop.
  */
  printf ("\nValue of toto in the second loop :\n");
  foreach()
    printf("%g ", toto[]);
}

/**
We make it so the programm ends after only one time step, no need to go further to see the bug.
*/
event end (t = 0.00001) {
  printf ("\n");
}

/**
We get the following results :

`Value of toto in the first loop :`
`-0.680375 0.211234 -0.566198 -0.59688 `

`Value of toto in the second loop :`
`6.92254e+235 -3.03422e-189 1.26266e+43 -0.59688`
*/


/**
The problem here is that when compiling (even with the `-Wall -Wextra` options), no error or warning is raised.
*/

/**
Important note :
This bug arises only when `datasize < 8`.
*/

