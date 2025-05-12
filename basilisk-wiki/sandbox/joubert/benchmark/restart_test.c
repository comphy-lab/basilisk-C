/**
# Simple test of Basilisk View

This is used to see the evolution of the size and time necessary to
read a restart file and generate a bigger dump file 
when increasing maxlevel. */

#include "fractions.h"
#include "utils.h"
#include "output.h"

int maxlevel = 4;

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);

  /**
  We first define a volume fraction field. */
  
  init_grid (1 << maxlevel);
  origin (-0.5,-0.5,-0.5);
  scalar f[];
  fraction (f, sq(x) + sq(y) + sq(z) - sq(0.3));
  if (restore (file = "restart")) {
    scalar pid[];
    foreach()
      pid[] = fmod(pid()*(npe() + 37), npe());
    boundary ({pid});
    /**
      Testing of the dump generation. */
    dump ();
  }
}
