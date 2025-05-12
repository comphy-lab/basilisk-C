/**
# A minimal "game" with Basilisk
   
   This boring game solves a poisson equation with an adaptive
   refinement criterion ($\zeta$). It uses `=` and `-` to tune it and
   `q` to stop the input.
*/

#include "poisson.h"
#include "view.h"
#include <GL/glut.h>
/**
We declare some function prototypes that are in cur.c
 */
void start_screen (void);
void stop_screen (void);
int getch (void);

double zeta = 0.005;
scalar f[];
f[left] = dirichlet (0.);
f[right] = dirichlet (0.);

f[bottom] = dirichlet (0.);
f[top] = dirichlet (0.);

/**
A function that is called for drawing on the screen is here called
`display()`.
*/

void display (void) {
  view (fov = 19);
  squares ("f");
  cells(lw = 0.3);
  glutSwapBuffers();
}

/**
A function that is called when the window in idle is here called
`update()`.
*/

void _update_(void) {
  scalar s[];
  foreach()
    s[] = exp (-(sq(x - 0.1) + sq(y + 0.05))*20);
  boundary ({s});
  
  char m = getch();
  if (m == '=')
    zeta *= 1.5;
  if (m == '-')
    zeta /= 1.3;
  if (m == 'q')
    stop_screen();
  poisson (f, s);
  adapt_wavelet ({f}, (double[]){zeta}, 99, 3, 10);
  /**
After this, we update the window using the `glut` syntax.
   */
  glutPostRedisplay();
}
/**
The main function sets up the windows for drawing and starts the
`glut` main loop.
 */

int main (int argc, char** argv ) {
  L0 = 2;
  X0 = Y0 = -L0/2.;
  init_grid (N);
  start_screen();
  glutInit (&argc, argv);
  glutInitDisplayMode (GLUT_DOUBLE);   
  glutInitWindowSize (700,700);        
  glutInitWindowPosition (100,100);    
  glutCreateWindow ("Boring Game"); 
  glutDisplayFunc (display);  
  glutIdleFunc (_update_);
  glutMainLoop();
}
