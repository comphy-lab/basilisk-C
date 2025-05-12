/**
#AMR with up/down-sampling


Attempt to recreate figures made by Antoon [here](http://basilisk.fr/sandbox/Antoonvh/the_adaptive_wavelet_algorithm)

Seems succesful.

*/

/**
~~~gnuplot Comparison coarse/fine grid on two finest levels of refinement
f(x) = exp(-x*x)
plot 'log1' u 1:2 w p lc 'black' pt 7 t 'Coarse grid',\
     'log2' u 1:2 w p lc 'blue' pt 7  t 'Leaves',\
      '' u 1:4 w p lc 'red' pt 7 t 'Interpolated fine value',\
     f(x) w l dt 2 lc 'black'
~~~


~~~gnuplot Error
plot 'log1' u 1:3,\
    'log2' u 1:3
~~~

##The code
Code is extracted from the tree-common.h

*/


#include "view.h"


int main(){
  L0 = 8;
  origin (-4., -4.);
  int MAXLEVEL = 4;
  init_grid(1<<MAXLEVEL);

  scalar s[];
  s.prolongation = refine_bilinear;
  // s.prolongation = refine_linear;
  s.restriction = restriction_average;
  foreach(){
    // s[] = 3*x*x*x+x*x+x+1; // 3x³ + x² + x + 1 
    s[] = exp(-sq(x)); 
  }
  // boundary({s});
  restriction({s});

  double refine_max = 1.; // this will mimic a threshold on the Hessian for the
  // refinement.


  scalar e[], e1[], e2[], e3[];

/** 
A set of variables needed for the algo...
*/
  astats st = {0, 0};
  scalar * listc = NULL;
  static const int refined = 1 << user;

  foreach_cell() {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
        if (cell.flags & too_coarse) {
          cell.flags &= ~too_coarse;
          refine_cell (point, listc, refined, &tree->refined);
          st.nf++; /*what is that ?*/
        }
        continue;
      }
      else { // !is_leaf (cell)
        if (cell.flags & refined) {
    // cell has already been refined, skip its children
          cell.flags &= ~too_coarse;
          continue;
        }
  // check whether the cell or any of its children is local
        bool local = is_local(cell);
        if (!local)
          foreach_child()
            if (is_local(cell)) {
              local = true; break;
            }
        if (local) {

          double sc[1 << dimension];
          int c = 0;
          foreach_child()
            sc[c++] = s[];

          s.prolongation (point, s);
          c = 0;
          foreach_child() {
            e[]  = fabs(sc[c] - s[]); // calculation of the error
            e1[] = s[]; // here we save the interpolated value using only
            // the values from the coarsest grid
            e2[] = (e[] - refine_max); // threshold for refinement
            e3[] = (e[] - refine_max/1.5); // threshold for coarsening    
            // 
            s[] = sc[c++];
          }
        }
      }
    }
    else // inactive cell
      continue;
  }

/**
One must notice that the *e1* field values come from the parent cells
because it is modified in a foreach_child() structure.
*/


/**
##Let's write some stuff !

First we write values on the grid where level == (MAXLEVEL-1)
*/
  char filename [100];
  snprintf(filename, 100,  "log1");
  FILE * fp1;    
  fp1 = fopen (filename,"w");
  foreach_cell(){
    if(y==Delta/2. && level == MAXLEVEL -1){
      fprintf(fp1, "%g %g %g %g %g %g %g %g\n", x, s[],e[], e1[], e2[], e3[],
       refine_max,refine_max/2.);
    }
  }
  fclose(fp1);
  snprintf(filename, 100,  "log2");
  fp1 = fopen (filename,"w");

/**
And then values on the finest grid
*/
  foreach(){
    if(y==Delta/2.){
      fprintf(fp1, "%g %g %g %g %g %g %g %g\n", x, s[],e[], e1[], e2[], e3[],
       refine_max,refine_max/2.);
    }
  }
  fclose(fp1);
  squares("e");
  save("e.png");
  squares("e2", min = 0);
  save("e2.png");
  squares("e3", min = 0);
  save("e3.png");

  return 0;
}