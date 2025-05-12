/**
# Drawing a space-curve

In this example, we define a space-curve $\mathcal{C}(\xi,t)$ and draw it.

<table>
<tr>
<td><center>![curve 1: a ring](test_draw_filaments/prescribed_curve1.png){ width="75%" }</center></td>
<td><center>![curve 2: a section of an helix](test_draw_filaments/prescribed_curve2.png){ width="75%" }</center></td>
</tr>
</table>

*/

#include "grid/multigrid3D.h"
#include "view.h"
#include "draw_filaments.h"

int main()
{
  L0 = 2*pi;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << 6;
  init_grid(N);

  int nseg = 128;
  double R=1.0, H=1.0;
  double dtheta = 2*pi/((double)nseg-1);
  double theta[nseg];
  double a[nseg];
  coord C[nseg];

  for (int i = 0; i < nseg; i++){
    theta[i] = dtheta * (double)i - 2*pi;
    C[i].x = R * cos(theta[i]);
    C[i].y = R * sin(theta[i]);
    C[i].z = 0;
    a[i] = 0.1 + 0.025*cos(2.0*theta[i]);
  } 

  {
    view (camera="iso");
    draw_tube_along_curve(nseg, C, a);
    save ("prescribed_curve1.png");
  }

  for (int i = 0; i < nseg; i++){
    theta[i] = 2.0 * dtheta * (double)i - 2*pi;
    C[i].x = R * cos(theta[i]);
    C[i].y = R * sin(theta[i]);
    C[i].z = (H/(2*pi)) * theta[i];
    a[i] = 0.1 + 0.025*cos(2.0*theta[i]);
  }

  {
    view (camera="iso");
    draw_tube_along_curve(nseg, C, a);
    save ("prescribed_curve2.png");
  }
}