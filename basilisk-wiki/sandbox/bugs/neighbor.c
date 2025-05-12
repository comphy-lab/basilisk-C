/**
# foreach_neighbor() is not properly nested */

scalar a[];

int main()
{
  init_grid (1);
  scalar * interfaces = {a};
  foreach() {
    
    /**
    This works... */
    
    foreach_neighbor() {
      for (scalar s in interfaces)
	fprintf (stderr, "a %g %g\n", x, y);
    }
    
    /**
    ... but this doesn't */
    
    foreach_neighbor()
      for (scalar s in interfaces)
	fprintf (stderr, "b %g %g\n", x, y);
  }
}

/**
~~~gnuplot
set size ratio -1
set xrange [-2:3]
set yrange [-2:3]
set key outside
plot '< grep a log' u 2:3 w p pt 5 t 'with braces', \
     '< grep b log' u 2:3 w p pt 7 t 'without braces'
~~~
*/