/**
# Drawing a space-curve

In this example, we define a space-curve $\mathcal{C}(\xi,t)$ and 
compute a local Frenet-Serret basis ($\bf\hat{T}, \hat{N}, \hat{B}$)

<table>
<tr>
<td><center>![prescribed curve](test_draw_filaments2/prescribed_curve.png){ width="75%" }</center></td>
<td><center>![Frenet-Serret frame](test_draw_filaments2/prescribed_curve_with_vectors.png){ width="75%" }</center></td>
</tr>
</table>

*/

#include "grid/multigrid3D.h"
#include "view.h"
#include "filaments.h"
#include "draw_filaments.h"

int main()
{
  L0 = 2*pi;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << 6;
  init_grid(N);

  int turns = 3;
  int nseg_per_turn = 32;
  int nseg = nseg_per_turn*turns+1;
  double R=1.0;
  double H=2.0;
  double dtheta = 2*pi/((double)nseg_per_turn);
  double theta[nseg];
  double a[nseg];
  coord C[nseg], dC[nseg], d2C[nseg];

  // Define a curve 
  for (int i = 0; i < nseg; i++){
    theta[i] = dtheta * (double)i - 2*pi;
    C[i].x = R * cos(theta[i]);
    C[i].y = R * sin(theta[i]);
    C[i].z = (H/(2*pi)) * theta[i] - H;
    a[i] = 0.1 + 0.01*cos(3.0*theta[i]);
  } 

  {
    view (camera="iso");
    draw_tube_along_curve(nseg, C, a);
    save ("prescribed_curve.png");
  }

  // Find the first and second derivatives to compute a Frenet-Serret Frame
  coord xshift = {0, 0, turns*H}, dxshift = {0, 0, 0};
  fd_derivative(nseg, dtheta,  xshift,  C,  dC);
  fd_derivative(nseg, dtheta, dxshift, dC, d2C);

  double sigma[nseg], ell[nseg];
  coord Tvec[nseg], Nvec[nseg], Bvec[nseg];
  
  for (int i = 0; i < nseg; i++){
    sigma[i] = 0.;    
    foreach_dimension(){      
      sigma[i] += sq(dC[i].x);
      Tvec[i].x =  dC[i].x/sqrt(vecdot( dC[i],  dC[i]));
      Nvec[i].x = d2C[i].x/sqrt(vecdot(d2C[i], d2C[i]));
    }
    sigma[i] = sqrt(sigma[i]);    
    Bvec[i] = vecdotproduct(Tvec[i], Nvec[i]);
  }

  ell[0] = 0.;
  for (int i = 0; i < nseg-1; i++){
    ell[i+1] = ell[i] + sigma[i+1]*dtheta;
  }

  view (camera="iso");  
  draw_space_curve_with_vectors(nseg, C, Tvec, Nvec, Bvec, scale=0.33);   
  save ("prescribed_curve_with_vectors.png");

  

  FILE *fp = fopen("curve.txt", "w"); 
  for (int i = 0; i < nseg; i++){
    fprintf (fp, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", i, theta[i], 
         C[i].x,    C[i].y,    C[i].z, 
        dC[i].x,   dC[i].y,   dC[i].z, 
       d2C[i].x,  d2C[i].y,  d2C[i].z,       
      Tvec[i].x, Tvec[i].y, Tvec[i].z,
      Nvec[i].x, Nvec[i].y, Nvec[i].z,
      Bvec[i].x, Bvec[i].y, Bvec[i].z);
  }
}

/** 

## Outputs

~~~gnuplot Coordinates of the prescribed curve
set term pngcairo enhanced size 600,480 font ",12"
set output 'curve.png'
plot 'curve.txt' u 1:3 w lp title "C.x",\
              "" u 1:4 w lp title "C.y",\
              "" u 1:5 w lp title "C.z"
set xlabel "index"
set ylabel "value"
~~~

~~~gnuplot First derivative
set term pngcairo enhanced size 600,480 font ",12"
set output 'first_derivative.png'
plot 'curve.txt' u 1:6 w lp title "dC.x",\
              "" u 1:7 w lp title "dC.y",\
              "" u 1:8 w lp title "dC.z"
set xlabel "index"
set ylabel "value"
~~~

~~~gnuplot Second derivative
set term pngcairo enhanced size 600,480 font ",12"
set output 'second_derivative.png'
plot 'curve.txt' u 1:9  w lp title "d2C.x",\
              "" u 1:10 w lp title "d2C.y",\
              "" u 1:11 w lp title "d2C.z"
set xlabel "index"
set ylabel "value"
~~~

~~~gnuplot Unit Tangent vector
set term pngcairo enhanced size 600,480 font ",12"
set output 'tangent_vector.png'
plot 'curve.txt' u 1:12  w lp title "T.x",\
              "" u 1:13 w lp title "T.y",\
              "" u 1:14 w lp title "T.z"
set xlabel "index"
set ylabel "value"
~~~

~~~gnuplot Unit Normal vector
set term pngcairo enhanced size 600,480 font ",12"
set output 'normal_vector.png'
plot 'curve.txt' u 1:15  w lp title "N.x",\
              "" u 1:16 w lp title "N.y",\
              "" u 1:17 w lp title "N.z"
set xlabel "index"
set ylabel "value"
~~~

~~~gnuplot Unit Binormal vector
set term pngcairo enhanced size 600,480 font ",12"
set output 'binormal_vector.png'
plot 'curve.txt' u 1:18  w lp title "B.x",\
              "" u 1:19 w lp title "B.y",\
              "" u 1:20 w lp title "B.z"
set xlabel "index"
set ylabel "value"
~~~
*/