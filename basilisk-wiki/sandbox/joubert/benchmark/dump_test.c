/**
# Simple test of generation of dump

This is used to check the evolution of the size and time necessary to
generate a dump file when increasing maxlevel. */

#include "fractions.h"
#include "utils.h"
#include "output.h"

int maxlevel = 5;

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
  /**
  Testing of the dump generation. */
  dump ();
}
