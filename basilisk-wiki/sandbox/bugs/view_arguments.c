/**                                                                           # View not recognizing tx,ty arguments                                          
When setting the camera position with the parameters tx and ty, if a zero value is specified no argument is taken for those parameters. When view is called multiple times this may generate errors in the camera position. If in a previous call the camera is displaced from the default centered position, in successive calls assigning $$tx = 0$$ does not move back the camera to the centered position.

<div class="message">
This has been fixed by [this patch](src/?changes=20230923152223).</div>
*/

#include "fractions.h"
#include "view.h"

int main() {

  /**                                                                        
  We first define a volume fraction field. */

  init_grid (16);
  origin (-0.5,-0.5,-0.5);
  scalar f[];
  fraction (f, sq(x) + sq(y) + sq(z) - sq(0.3));

  /**                                                                        
  We set a first view displacing the camera in the x direction. */
  
  view (tx = -0.2, ty = 0., width = 400, height = 400);
  box();
  draw_vof ("f");
  save ("out.png");

  /**
  ![](view_arguments/out.png)
   
  We set a second view centered on the origin. */
  
  view (tx = 0., ty = 0., width = 400, height = 400);
  box();
  draw_vof ("f");
  save ("out2.png");

  /**
  ![View does not take the new argument](view_arguments/out2.png)
   
  A simple workaround is to specify a very small number instead of zero. */
  
  view (tx = 1.e-12, ty = 0., width = 400, height = 400);
  box();
  draw_vof ("f");
  save ("out3.png");

  /**
  ![With the workaround](view_arguments/out3.png) */

}