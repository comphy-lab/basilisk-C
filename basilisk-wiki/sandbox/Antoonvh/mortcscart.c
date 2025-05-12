/**
# Morton versus Cartesian indexing
On this page we test the socalled "locality" of a morton-style interator versus a regular Cartesian iterator. We test it for a $32^3$ grid and assume that our results are somehow scaleable to other grids. We focus the analysis on the distance in iteration-sequence space between cells that are neighbors in the 3D-grid space   
*/
#include "grid/octree.h"

scalar m[],xyz[];
int main(){
  int maxlevel = 5;
  int cells = pow(2,3*maxlevel); 
  double cartarr[cells];
  double mortarr[cells];
  int i=0;
  while (i<cells){
    cartarr[i]=0.;
    mortarr[i]=0.;
    i++;
  }
  /**
  The grid is initialized and we use the `foreach()` loop for the Morton-style iteration and store the indixes in field *s*.  
  */
  init_grid(1<<maxlevel);
  int a=1;
  foreach()
    m[]=a++;
  a=1;
  int o = 1+BGHOSTS;
  /**
 For the Cartesian-style indexing, we define a $x-y-z$-sequence iterator. The result is stored in the *xyz* field.
  */
  for (int k= o; k<N+o; k++){
    for (int j= o; j<N+o; j++){
      for (int i= o; i<N+o; i++){
	Point point;
	point.i=i; point.j=j; point.k=k; point.level=maxlevel;
	xyz[]=a++;
      }
    } 
  }
  double distcart=0;
  double distmort=0;
  /**
  Below we perform our analysis. Foreach cell we log the distance to three of its face-sharing neighbors. The boundary-ghost-cell values are set using the default scalar-field boundary condtion. This way, the index distance to ghost cells is 0 and does not 'pollute' the results. We define a total distance that is the sum of all individual neighbours' distances.  
  */
  boundary(all);
  foreach(){
    double cart=xyz[];
    double mort=m[];
    foreach_dimension(){
      double cd = fabs(cart-xyz[1,0,0]);
      double md = fabs(mort-m[1,0,0]);
      distcart += cd;
      distmort += md;
       /**
    Remarkably(?), the total 'index distance' to neighbors is exactly equal for both approaches (i.e. $N^2(N^2(N-1))+N(N^2(N-1))+(N^2(N-1))\propto N^5)$ for a $N^3$ grid). Therefore we check the underlying distribution of the indexing distances. 
      */
      cartarr[(int)(cd+0.5)]++;
      mortarr[(int)(md+0.5)]++;
    }
  }
  FILE * fp = fopen("hist","w");
  i=1;
  while (i<cells){
    fprintf(fp,"%d\t%g\t%g\n",i,cartarr[i],mortarr[i]);
    i++;
  }
  /**
  Below, the histrogram of the index distances between neighbors is shown:
~~~gnuplot Notice the logaritmic x and y axes
set xr [0.5 :100001]
set yr [0.5 :80000]
set logscale y
set logscale x
set xlabel 'Indexing distance'
set ylabel 'Number of occurences'
  plot 'hist' u 1:2 w lines lw 3 t 'Cartesian-style' ,\
  'hist' u 1:3 w lines lw 3 t 'Morton-style'
  ~~~
  
   We can see that the Cartesian style indexing has resulted in three-values for its index distances. ($1,N$ and $N^2$). The Morton-style curve has index distance values for each power of two. Compared to the Cartesian-style index distances; there are more neigbors with an index distance smaller than $N$, but the Morton-style indexing pays with more cells at an index distance larger than $N^2$. Remember, the (total) first order moment associated with each histrogram is equal!   
  */
}