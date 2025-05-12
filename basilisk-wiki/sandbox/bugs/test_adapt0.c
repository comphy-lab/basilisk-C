/**
# Testing adaption 


When I run test_adapt0.c with max level 16 this is what I get:

~~~
laneem@niwa-36410:~/basilisk/update/basilisk/src/examples> make test_adapt0.tst
qcc -O2 -fopenmp -g -Wall -o test_adapt0/test_adapt0 test_adapt0.c -lm
[test_adapt0.tst]
laneem@niwa-36410:~/basilisk/update/basilisk/src/examples> more test_adapt0/log 
t i h.min h.max h.sum u.x.rms u.x.max dt
0 0 10 10 1.2365e+11 0 0 0
# refined 32 cells, coarsened 0 cells
274.72 1 10 10 1.2365e+11 0 0 274.72
# refined 84 cells, coarsened 0 cells
412.079 2 10 10 1.2365e+11 0 0 137.36
# refined 108 cells, coarsened 0 cells
480.759 3 10 10 1.2365e+11 0 0 68.6799
# refined 374 cells, coarsened 0 cells
515.099 4 10 10 1.2365e+11 0 0 34.34
# refined 532 cells, coarsened 0 cells
532.269 5 10 10 1.2365e+11 0 0 17.17
# refined 1448 cells, coarsened 8 cells
540.854 6 10 10 1.2365e+11 0 0 8.58499
# refined 2028 cells, coarsened 0 cells
545.147 7 10 10 1.2365e+11 0 0 4.2925
# refined 5846 cells, coarsened 0 cells
547.293 8 10 10 1.2365e+11 0 0 2.14625
# refined 8212 cells, coarsened 0 cells
548.366 9 10 10 1.2365e+11 0 0 1.07312
# refined 23336 cells, coarsened 8 cells
548.903 10 10 10 1.2365e+11 1.45462e-20 1.40433e-17 0.536562
# refined 32748 cells, coarsened 0 cells
549.171 11 10 10 1.2365e+11 2.18193e-20 2.1065e-17 0.268281
# refined 26236 cells, coarsened 0 cells
549.439 12 10 10 1.2365e+11 2.90924e-20 2.80867e-17 0.268281
# refined 0 cells, coarsened 0 cells
549.707 13 10 10 1.2365e+11 3.63655e-20 3.51083e-17 0.268281
# refined 0 cells, coarsened 0 cells
549.976 14 10 10 1.2365e+11 4.36386e-20 4.213e-17 0.268281
# refined 0 cells, coarsened 0 cells
550.244 15 10 10 1.2365e+11 5.29489e-20 2.24693e-16 0.268281
# refined 0 cells, coarsened 0 cells
~~~

And the results is what I would expect - refined to level 16 around the edges of the specified box.

However if I change the maxlevel to 17 I get this:

~~~
laneem@niwa-36410:~/basilisk/update/basilisk/src/examples> make test_adapt0.tst
qcc -O2 -fopenmp -g -Wall -o test_adapt0/test_adapt0 test_adapt0.c -lm
[test_adapt0.tst]
t i h.min h.max h.sum u.x.rms u.x.max dt
0 0 10 10 1.2365e+11 0 0 0
# refined 32 cells, coarsened 0 cells
274.72 1 10 10 1.2365e+11 0 0 274.72
# refined 84 cells, coarsened 0 cells
412.079 2 10 10 1.2365e+11 0 0 137.36
# refined 108 cells, coarsened 0 cells
480.759 3 10 10 1.2365e+11 0 0 68.6799
# refined 374 cells, coarsened 0 cells
515.099 4 10 10 1.2365e+11 0 0 34.34
# refined 532 cells, coarsened 0 cells
532.269 5 10 10 1.2365e+11 0 0 17.17
# refined 1448 cells, coarsened 8 cells
540.854 6 10 10 1.2365e+11 0 0 8.58499
# refined 2028 cells, coarsened 0 cells
545.147 7 10 10 1.2365e+11 0 0 4.2925
# refined 5846 cells, coarsened 0 cells
547.293 8 10 10 1.2365e+11 0 0 2.14625
# refined 8212 cells, coarsened 0 cells
548.366 9 10 10 1.2365e+11 0 0 1.07312
# refined 23336 cells, coarsened 8 cells
548.903 10 10 10 1.2365e+11 1.45462e-20 1.40433e-17 0.536562
# refined 32748 cells, coarsened 0 cells
549.171 11 10 10 1.2365e+11 2.18193e-20 2.1065e-17 0.268281
/home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652:error: [Thread debugging using libthread_db enabled]
/home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652:error: Using host libthread_db library "/lib64/libthread_db.so.1".
/home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652:error: [New Thread 0x2aaaabaf9700 (LWP 9638)]
/home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652:error: [New Thread 0x2aaaabcfa700 (LWP 9639)]
/home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652:error: [New Thread 0x2aaaabefb700 (LWP 9640)]
/home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652:error: [New Thread 0x2aaaac0fc700 (LWP 9641)]
/home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652:error: [New Thread 0x2aaaac2fd700 (LWP 9642)]
/home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652:error: [New Thread 0x2aaaac4fe700 (LWP 9643)]
/home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652:error: [New Thread 0x2aaaac6ff700 (LWP 9644)]
/home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652:error: Program received signal SIGSEGV, Segmentation fault.
make: *** [test_adapt0.tst] Error 1
~~~

Oddly it seems to have crashed at a lower level of refinement than the early run.

It seems to have crashed while trying to allocate children - this is what I get from the core dump:

~~~
laneem@niwa-36410:~/basilisk/update/basilisk/src/examples/test_adapt0> gdb ./test_adapt0 core
GNU gdb (GDB) SUSE (7.5.1-0.7.29)
Copyright (C) 2012 Free Software Foundation, Inc.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.  Type "show copying"
and "show warranty" for details.
This GDB was configured as "x86_64-suse-linux".
For bug reporting instructions, please see:
<http://www.gnu.org/software/gdb/bugs/>...
Reading symbols from /home/laneem/basilisk/update/basilisk/src/examples/test_adapt0/test_adapt0...done.
Core was generated by `./test_adapt0'.
Program terminated with signal 11, Segmentation fault.
#0  0x000000000040b595 in alloc_children (j=<optimized out>, i=<optimized out>, point=...)
    at /home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652
652              assert (!CHILD(k,l));
(gdb) where
#0  0x000000000040b595 in alloc_children (j=<optimized out>, i=<optimized out>, point=...)
    at /home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:652
#1  increment_neighbors (point=...) at /home/laneem/basilisk/update/basilisk/src/grid/quadtree.h:678
#2  0x000000000040ba7f in refine_cell (point=..., list=0x62ed90, flag=4, refined=0x623950)
    at /home/laneem/basilisk/update/basilisk/src/grid/quadtree-common.h:50
#3  0x000000000041d620 in adapt_wavelet (p=...) at /home/laneem/basilisk/update/basilisk/src/grid/quadtree-common.h:188
#4  0x000000000041db7d in adapt () at test_adapt0.c:53
#5  0x000000000041dbbb in do_adapt (i=-1, t=0, _ev=0x4000000) at test_adapt0.c:111
#6  0x000000000040358f in event_do (action=<optimized out>, t=<optimized out>, i=<optimized out>, ev=<optimized out>)
    at /home/laneem/basilisk/update/basilisk/src/grid/events.h:140
#7  events (i=11, t=549.17088372626267, action=true) at /home/laneem/basilisk/update/basilisk/src/grid/events.h:187
#8  0x000000000041c491 in run () at /home/laneem/basilisk/update/basilisk/src/predictor-corrector.h:56
#9  0x000000000041c62f in main () at test_adapt0.c:32
~~~
      
It seems to be related to memory but it is not obviously maxing it out at least from a cursory watching top while it is running point of view. If I make the square it refines bigger I can only run to level 15 before it crashes. But is other situations I have run simulations with far more actual points in them.

Using the time command gives
      
~~~
laneem@niwa-36410:~/basilisk/update/basilisk/src/examples/test_adapt1> \time ./test_adapt1 
Command terminated by signal 11
2.96user 2.55system 0:04.46elapsed 123%CPU (0avgtext+0avgdata 2222404maxresident)k
0inputs+4443448outputs (0major+578832minor)pagefaults 0swaps
~~~
  
or with -p flag
  
~~~
real 7.73
user 2.96
sys 2.84
~~~
  
Using valgrind
  
~~~
laneem@niwa-36410:~/basilisk/update/basilisk/src/examples/test_adapt1> valgrind --main-stacksize=16000000 ./test_adapt1 
==14543== Memcheck, a memory error detector
==14543== Copyright (C) 2002-2012, and GNU GPL'd, by Julian Seward et al.
==14543== Using Valgrind-3.8.1 and LibVEX; rerun with -h for copyright info
==14543== Command: ./test_adapt1
==14543== 
# refined 32 cells, coarsened 0 cells
# refined 84 cells, coarsened 0 cells
# refined 108 cells, coarsened 0 cells
# refined 374 cells, coarsened 0 cells
# refined 532 cells, coarsened 0 cells
# refined 1448 cells, coarsened 8 cells
# refined 2028 cells, coarsened 0 cells
# refined 5846 cells, coarsened 0 cells
# refined 8212 cells, coarsened 0 cells
# refined 23336 cells, coarsened 8 cells
# refined 32748 cells, coarsened 0 cells
==14543== Invalid read of size 8
==14543==    at 0x40B515: increment_neighbors (quadtree.h:652)
==14543==    by 0x40B9FE: refine_cell (quadtree-common.h:50)
==14543==    by 0x41D17F: adapt_wavelet (quadtree-common.h:188)
==14543==    by 0x41D6DC: adapt (test_adapt1.c:50)
==14543==    by 0x41D71A: do_adapt (test_adapt1.c:109)
==14543==    by 0x40358E: events (events.h:140)
==14543==    by 0x41BFF0: run (predictor-corrector.h:56)
==14543==    by 0x41C18E: main (test_adapt1.c:29)
==14543==  Address 0x8ccd0 is not stack'd, malloc'd or (recently) free'd
==14543== 
==14543== 
==14543== Process terminating with default action of signal 11 (SIGSEGV): dumping core
==14543==  Access not within mapped region at address 0x8CCD0
==14543==    at 0x40B515: increment_neighbors (quadtree.h:652)
==14543==    by 0x40B9FE: refine_cell (quadtree-common.h:50)
==14543==    by 0x41D17F: adapt_wavelet (quadtree-common.h:188)
==14543==    by 0x41D6DC: adapt (test_adapt1.c:50)
==14543==    by 0x41D71A: do_adapt (test_adapt1.c:109)
==14543==    by 0x40358E: events (events.h:140)
==14543==    by 0x41BFF0: run (predictor-corrector.h:56)
==14543==    by 0x41C18E: main (test_adapt1.c:29)
==14543==  If you believe this happened as a result of a stack
==14543==  overflow in your program's main thread (unlikely but
==14543==  possible), you can try to increase the size of the
==14543==  main thread stack using the --main-stacksize= flag.
==14543==  The main thread stack size used in this run was 16003072.
==14543== 
==14543== HEAP SUMMARY:
==14543==     in use at exit: 6,372,878,906 bytes in 15,153 blocks
==14543==   total heap usage: 54,237 allocs, 39,084 frees, 27,198,217,954 bytes allocated
==14543== 
^C==14543== LEAK SUMMARY:
==14543==    definitely lost: 0 bytes in 0 blocks
==14543==    indirectly lost: 0 bytes in 0 blocks
==14543==      possibly lost: 4,737,741,004 bytes in 13,556 blocks
==14543==    still reachable: 1,635,137,902 bytes in 1,597 blocks
==14543==         suppressed: 0 bytes in 0 blocks
==14543== Rerun with --leak-check=full to see details of leaked memory
==14543== 
==14543== For counts of detected and suppressed errors, rerun with: -v
==14543== ERROR SUMMARY: 1 errors from 1 contexts (suppressed: 4 from 4)
Killed
~~~


I need to work out what this means
*/

#include "spherical.h"
#include "saint-venant.h"

/**
We then define a few useful macros and constants. */

#define MAXLEVEL 20
#define MINLEVEL 5
#define ETAE     1e-2 // error on free surface elevation (1 cm)
#define LON0  0.
#define LAT0   0.
#define DOMAIN_SIZE 1.

int main()
{

  Radius = 6371220.;
  // the domain is 1 degrees squared
  size (DOMAIN_SIZE);
  // centered on 0,0 longitude,latitude
  origin (LON0 - L0/2.,LAT0  - L0/2.);

  init_grid (1 << MINLEVEL);

  /**
  We then call the *run()* method of the Saint-Venant solver to
  perform the integration. */

  run();
}

scalar lim[];

/**
## Adaptation
*/
int adapt() {
  foreach(){
    double bound = L0*0.05;
    double xlim = fabs(x-LON0) < bound ? 1. : 0.;
    double ylim = fabs(y-LAT0) < bound ? 1. : 0.;
    lim[] = xlim * ylim ;}
  boundary ({lim});

  /**
  We can now use wavelet adaptation 
  The function then returns the number of cells refined. */

  astats s = adapt_wavelet ({lim}, (double[]){ETAE},
			    MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
}


event init (i = 0){
  foreach()
  zb[]=-10;
  conserve_elevation();

  /**
  The initial still water surface is at $z=0$ so that the water depth
  $h$ is... */

  foreach()
    h[] = max(0., - zb[]);
  boundary ({h});

  }

/**
## Outputs

### At each timestep

We output simple summary statistics for *h* and *u.x* on standard
error. */

event logfile (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g\n", t, i, s.min, s.max, s.sum, 
	   n.rms, n.max, dt);

}

/**
### Snapshots
*/

event snapshots (i++; i <= 15) {
  /**
  We save snapshot files along the way. */

  char *outfile = NULL;
  outfile = (char *) malloc(sizeof(char) * 16);
  sprintf(outfile, "snapshot-%d.gfs", i);
  output_gfs (file = outfile, t = t);
}

/**
## Adaptivity

And finally we apply our *adapt()* function at every timestep. */

event do_adapt (i++) adapt();