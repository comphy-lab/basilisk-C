/**
# Render Basilisk in realtime

It is possible to draw bview commands to a window at high framerates.

##`GLUT`

For this purpose we can use the `free GLUT` library. Perhaps this
works:

~~~literatec
$sudo apt install freeglut3-dev
~~~

## using GLUT with Basilisk's events 

Now the inconvienence starts: The all-important `glutMainLoop()` does
not seem to be compatible with Basilisk's `run()` loop. An example
code to work with this issue could look like this:

 */
#include "navier-stokes/centered.h"
#include "view.h"

#include <GL/glut.h>
/**
We define a function that determines what is being drawn in the window.
 */
void display() {
  scalar omega[];
  vorticity (u, omega);
  view (fov = 19);
  squares ("omega", min = -2, max = 2);
  cells (lw = 0.5);
  glutSwapBuffers(); //Rather than `save()`
}
/**
   Further, we define a function that calls events. Note that `t` is
   not updated and it will remain 0 during this simulation.
*/
void call_events() {
  events (true);       //Do events 
  iter = inext;        //Update iteration
  glutPostRedisplay(); //Redraw the solution in the window
}
/**
The main function is tasked with setting up the GLUT window and
calling the glut main loop. This loop does not exit untill the window
in closed.
 */
int main (int argc, char** argv ) {
  L0 = 10;
  X0 = Y0 = -L0/2;
  init_grid (N);
  glutInit (&argc, argv);
  glutInitDisplayMode (GLUT_DOUBLE);   
  glutInitWindowSize (900,900);        
  glutInitWindowPosition (100,100);    
  glutCreateWindow ("Vortex-wall collision"); 
  glutDisplayFunc (display);  
  glutIdleFunc (call_events);
  glutMainLoop();
  return 0;
}
/**
## Basilisk's events

We use regular events, apart from the fact that we do not use the time
parameter (e.g.: `event init (t = 0)`)
 */
u.t[left] = dirichlet (0.);
#define RAD (pow(pow((x - xo), 2)+(pow((y - yo), 2)), 0.5))
#define ST ((y - yo)/RAD)
double k = 3.83170597, xo = 0.1, yo = 0.2;
const face vector muc[] = {0.001, 0.001};

event init (i = 0) {
  DT = 0.1;
  mu = muc;
  scalar psi[];
  refine (RAD < 2.0 && level <= 8);
  refine (RAD < 1.0 && level <= 9);
  foreach() 
    psi[] = ((RAD > 1)*((1/RAD))*ST) +
    ((RAD < 1)*((-2*j1(k*RAD)*ST/(k*j0(k))) + (RAD*ST)));
  boundary ({psi});
  foreach() {
    u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta));
    u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta);
  }
  boundary ({u.x, u.y});
}

event adapt (i++) 
  adapt_wavelet ({u.x, u.y}, (double[]){0.01, 0.01}, 8);

/**
## Compilation

Compile using the relevant libraries:

~~~literatec 
$ qcc -Wall realtime-rendering.c -lm -L$BASILISK/gl -lglutils \
      -lfb_glx -lGLU -lGLEW -lGL -lX11 -lGL -lglut
$ ./a.out
~~~
 */
