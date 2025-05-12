/**
# Adaptive Snake

Play adaptive snake:

![A highscore!](snake.snakecut2.mp4) 

## How to install and play?

Opun an actual request:

On a Debian machine it may take a few minutes:
  
* Install Basilisk with `bview`: (see: [Instructions](/src/INSTALL))
* Obtain the `ncurses` and free`glut` libraries.

~~~literatec
$ sudo apt install libncurses5-dev freeglut3-dev
~~~

* Copy the [`cur.c`](cur.c) code to a file `cur.c`
* Copy the [raw `snake.c`](snake.c?raw) code to a file `snake.c`  
* Now create two object files with gcc *and* qcc:

~~~literatec
$ gcc -c -o cur.o cur.c
$ qcc -c -o snake.o snake.c
~~~

* Good! Now compile this using the relevant libraries

~~~literatec
$ gcc -O2 -Wall snake.o cur.o -lm -lncurses -L$BASILISK/gl \
         -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lglut
~~~

* Move the terminal to the right-hand side of your screen and do:

~~~literatec
$ ./a.out
~~~

* Now click the terminal to let it register input. 

## The next version:

Adaptive Snake 2.0 exists; [see here](snake2.c).
 
 */

#include "view.h"

#include <GL/glut.h>

scalar s[], age[];

int maxlevel = 6;

bool stop_game = false;

int length;
coord targ;
int dir[2];
Point head, target;
double seconds_left, tot_time = 120;
timer time_elapsed;
double speed = 20;

void display() {
  view (fov = 21);
  if (stop_game) {
    draw_string ("Game Over", 1);
    char str[99];
    sprintf (str, "You scored: %d", (int)(length - 5));
    draw_string (str, 2);
    draw_string ("Close screen to exit", 3);
  } else {
    squares ("s", min = -2, max = 2);
    cells (lw = 0.5);
    char str[99];
    sprintf (str, "score: %d", length - 5);
    draw_string (str);
    draw_string ("Press `q` to stop", 3);
    sprintf (str, "Timer (sec): %.3g", seconds_left);
    draw_string (str, 1);
    box (notics = true, lc = {0.8,0,0.8});
  }
  glutSwapBuffers();
}

void update_game() {
  if (!stop_game) {
    events (true);
    iter = inext;
  } else {
    system ("sleep 0.5");
  }
  glutPostRedisplay();
}

static inline void coarsen_age (Point point, scalar s) {
  int max = 0;
  foreach_child()
    if (s[] > max)
      max = s[];
  s[] = max;
}

static inline void refine_age (Point point, scalar s) {
  int pc = s[];
    foreach_child()
      s[] = pc;
}

void start_screen (void);
int getch (void);
void stop_screen (void);

int main(int argc, char** argv ) {
  srand (time(NULL));
  perf.gt = timer_start();
  age.coarsen = coarsen_age;
  age.refine = refine_age;
  start_screen();
  foreach_dimension()
    periodic (left);
  L0 = 2;
  X0 = Y0 = -L0/2;
  init_grid (1 << maxlevel);
  glutInit (&argc, argv);
  glutInitDisplayMode (GLUT_DOUBLE);   
  glutInitWindowSize (1000,1000);        
  glutInitWindowPosition (100,100);    
  glutCreateWindow ("Adaptive Snake"); 
  glutDisplayFunc (display);  
  glutIdleFunc (update_game);
  glutMainLoop();
  stop_screen();
  return 0;
}

event init (i = 0) {
  speed = 20;
  init_grid (1 << 6);
  foreach()
    age[] = 0;
  head = locate (0, 0);
  Point point = head;
  age[] = 1;
  dir[0] = 1, dir[1] = 0;
  length = 5;
  targ.x = noise(), targ.y = noise();
  target = locate (targ.x, targ.y);
}

event obtain_movement (i++) {
  char m = getch();
  if (m == 'w' && dir[1] != -1)
    dir[0] = 0, dir[1] = 1;
  if (m == 'a' && dir[0] != 1)
    dir[0] = -1, dir[1] = 0;
  if (m == 's' && dir[1] != 1)
    dir[0] = 0, dir[1] = -1;
  if (m == 'd' && dir[0] != -1)
    dir[0] = 1, dir[1] = 0;

  if (m == '=')
    speed *= 2;
  if (m == '-')
    speed /= 2;

  if (m == 'q') {
    stop_game = true;
    stop_screen();
  }
  
  //  You may beta test version 2.0
  
  if (m == 'u') {
    refine (age[] < 2 && age[] > 0 && level <= maxlevel);
    maxlevel++;
  }
  
  if (m == 'j') {
    unrefine (level >= maxlevel);
    maxlevel--;
  }
}

event move (i++) {
  foreach() {
    if (age[] == 1) {
      if (age[dir[0], dir[1]] > 1) { // You have hit yourself
	event ("init");
	break;
      }
      if (is_leaf(neighbor(dir[0], dir[1])) && point.i == head.i) {
	age[dir[0],dir[1]] = 1;
      	break;
      } else {     //This is not supposed to happen
      	exit (0);
      }
    }
  }
  foreach() {
    if (age[] >= 1.)
      age[]++;
    if (age[] > length)
      age[] = 0;
  }
  head.i += dir[0];
  head.j += dir[1];
  if (head.i >= (1 << head.level) + GHOSTS)
    head.i -= (1 << head.level);
  if (head.j >= (1 << head.level) + GHOSTS)
    head.j -= (1 << head.level);
  if ( head.i < GHOSTS )
    head.i += (1 << head.level);
  if ( head.j < GHOSTS )
    head.j += (1 << head.level);
  Point point = head;
  age[] = 1;
  char str[99];
  sprintf (str, "sleep %g", 1./speed);
  system (str);
}

event adapt (i++) {
  boundary ({age});
  adapt_wavelet ({age}, (double[]){0.01}, maxlevel);
}

event eat_target (i++) {
  target = locate (targ.x, targ.y);
  if (head.i == target.i && head.j == target.j
      && head.level == target.level) {
    targ.x = noise(), targ.y = noise();
    length++;
  }
}

event do_timer (i++) {
  perf.t = timer_elapsed (perf.gt);
  seconds_left = tot_time - perf.t;
  if (seconds_left <= 0) {
    stop_game = true;
    stop_screen();
  }
  putc('.', stdout);
  fflush (stdout);
}

event compute_s (i++, last) {
  foreach() {
    s[] = nodata;
    if (age[] == 1) 
      s[] = 2;
    if (age[] > 1)
      s[] = 1;
  }
  Point point = locate (targ.x, targ.y);
  s[] = point.level - target.level;
  point = head;
  s[] = 10;
}
