/**
# Adaptive Snake 2.0  

Players `Red` and `Blue`, respectively, use `w a s d` and `i j k l` for
up, left, down, right, respectively. If a snake head runs into a snake,
it resets and exchanges 25 percent of its remaining life with the
other. Addional life duration can be obtained by finding the target
(green), which may be hiding at a lower level.

![Good game Red, Well done Blue!](snake2.mp4)

Adaptive Snake 2.0 is the two-player version of [Adaptive Snake
(1.0)](snake.c)

## How to install and play?

On a Debian machine it may take a few minutes:
  
* Install Basilisk with `bview`: (see: [Instructions](/src/INSTALL))
* Obtain the `ncurses` and free`glut` libraries.

~~~literatec
$ sudo apt install libncurses5-dev freeglut3-dev
~~~

* Copy the [`cur.c`](cur.c) code to a file `cur.c`
* Copy the [raw `snake2.c`](snake2.c?raw) code to a file `snake2.c`  
* Now create two object files with gcc *and* qcc:

~~~literatec
$ gcc -c -o cur.o cur.c
$ qcc -c -o snake2.o snake2.c
~~~

* Good! Now compile this using the relevant libraries

~~~literatec
$ gcc -O2 -Wall snake2.o cur.o -lm -lncurses -L$BASILISK/gl \
  -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lglut
~~~

* Move the terminal to the right-hand side of your screen and do:

~~~literatec
$ ./a.out
~~~

* Now click the terminal to let it register input, then press a key to
  start.

## Implementation
*/
#include "view.h"
  
#include <GL/glut.h>
#pragma autolink -L. -lcur -lncurses -lglut

scalar s[], age[];
int maxlevel = 6;

bool stop_game = false, the_beginning = true;

int length[2] = {5, 5}, loser = -1;
coord targ;
int dir[2][2];
Point head[2], target;
double seconds_left[2], tot_time[2] = {60, 60};
double speed = 20;

void display2() {
  view (fov = 21);
  if (the_beginning) { //This does not work..
    draw_string ("Pres any key to start...", 1);
  }
  else if (stop_game) {
    draw_string ("Game Over", 1);
    if (loser == 0)
      draw_string ("Well done Blue!", 2, 22, {0.1, 0.1, 0.9}, 2);
    else if (loser == 1)
      draw_string ("Well done Red!", 2, 22, {0.9, 0.1, 0.1}, 2);
    draw_string ("Close screen to exit", 3);
  } else {
    squares ("s", min = -2, max = 2);
    cells (lw = 0.5);
    char str[99];
    sprintf (str, "timer Red: %d", (int) (seconds_left[0] + 0.9));
    draw_string (str, 1);
    sprintf (str, "timer Blue: %d", (int) (seconds_left[1] + 0.9));
    draw_string (str, 2);
    draw_string ("Press `q` to stop", 3);
    box (notics = true, lc = {0.4, 0, 0.4});
  }
  glutSwapBuffers();
}

void update_game() {
  if (the_beginning) {
    getchar();
    perf.gt = timer_start();
    the_beginning = false;
  }
  if (!stop_game) {
      events (true);
      iter = inext;
    } else {
      system ("sleep 0.5");
  }
  glutPostRedisplay();
}

static inline void coarsen_age (Point point, scalar s) {
  s[] = 0.;
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
  age.coarsen = coarsen_age;
  age.refine = refine_age;
  start_screen();
  foreach_dimension()
    periodic (left);
  L0 = 2;
  X0 = Y0 = -L0/2;
  init_grid (1 << maxlevel);
  glutInit (&argc, argv);
  glutInitDisplayMode (GLUT_SINGLE);   
  glutInitWindowSize (900,900);        
  glutInitWindowPosition (100,50);    
  glutCreateWindow ("Adaptive Snake 2.0"); 
  glutDisplayFunc (display2);  
  glutIdleFunc (update_game);
  glutMainLoop();
  stop_screen();
  return 0;
}

void reset_snake (int j, double xp, double yp) {
  refine (sq(x - xp) + sq(y - yp) < sq(4*Delta) && level < maxlevel);
  if (j == 0) {
    foreach() 
      age[] = age[] > 0 ? 0 : age[];
    head[j] = locate (xp, yp);
    dir[0][j] = -1, dir[1][j] = 0;
    Point point = head[j];
    age[] = 1;
  } else {
    foreach() 
      age[] = age[] < 0 ? 0 : age[];
    head[j] = locate (xp, yp);
    dir[0][j] = 1, dir[1][j] = 0;
    Point point = head[j];
    age[] = -1;
  }
  length[j] = 5;
}

event init (i = 0) {
  speed = 20;
  init_grid (1 << 6);
  foreach()
    age[] = 0;
  reset_snake (0,  0.5,  0.5);
  reset_snake (1, -0.5, -0.5);
  targ.x = noise(), targ.y = noise();
  target = locate (targ.x, targ.y);
}

event obtain_movement (i++) {
  int ii = 0;
  char m[99];
  while ((m[ii] = getch()) != -1) // Upto 99 keystrokes per step...
    ii++;    
  int j = 0;
  
  while (j < ii) { //Red
    if ((m[j] == 'w')  && dir[1][0] != -1) {
      dir[0][0] = 0, dir[1][0] = 1; break; }
    if (m[j] == 'a' && dir[0][0] != 1) {
      dir[0][0] = -1, dir[1][0] = 0; break; }
    if (m[j] == 's' && dir[1][0] != 1) {
      dir[0][0] = 0, dir[1][0] = -1; break;}
    if (m[j] == 'd' && dir[0][0] != -1) {
      dir[0][0] = 1, dir[1][0] = 0; break; }
    
    if (m[j] == 'q') {
      stop_game = true;
      stop_screen();
    }
    j++;
  } 
  j = 0;
  while (j < ii) { //Blue
    if ((m[j] == 'i')  && dir[1][1] != -1) {
      dir[0][1] = 0, dir[1][1] = 1; break; }
    if (m[j] == 'j' && dir[0][1] != 1) {
      dir[0][1] = -1, dir[1][1] = 0; break; }
    if (m[j] == 'k' && dir[1][1] != 1) {
      dir[0][1] = 0, dir[1][1] = -1; break;}
    if (m[j] == 'l' && dir[0][1] != -1) {
      dir[0][1] = 1, dir[1][1] = 0; break; }
    j++;
  }
}

event move (i++) {
  bool resets[2] = {false, false};
  foreach() {
    if (age[] == 1) {
      if (age[dir[0][0], dir[1][0]] != 0 ) { // Dead 
	resets[0] = true;
	tot_time[0] -= seconds_left[0]/4.;
	tot_time[1] += seconds_left[0]/4.;
	break;
      }
      if (is_leaf(neighbor(dir[0][0], dir[1][0]))) {
	age[dir[0][0],dir[1][0]] = 1;
	break;
      }
    }
  }
  foreach() {
    if (age[] == -1) {
      if (age[dir[0][1], dir[1][1]] != 0 ) { // Dead 
	resets[1] = true;
	tot_time[1] -= seconds_left[1]/4.;
	tot_time[0] += seconds_left[1]/4.;
	break;
      }
      if (is_leaf(neighbor(dir[0][1], dir[1][1]))) {
	age[dir[0][1],dir[1][1]] = -1;
	break;
      }
    }
  }
  
  for (int l = 0; l < 2; l++) 
    if (resets[l])
      reset_snake (l, noise(), noise());

  foreach() {
    if (age[] > 0.)
      age[]++;
    if (age[] > 0 && age[] > length[0])
      age[] = 0;
    if (age[] < 0.)
      age[]--;
    if (age[] < 0 && age[] < -length[1])
      age[] = 0;
  }
  for (int l = 0; l < 2; l++) {
    head[l].i += dir[0][l];
    head[l].j += dir[1][l];
    if (head[l].i >= (1 << head[l].level) + GHOSTS)
      head[l].i -= (1 << head[l].level);
    if (head[l].j >= (1 << head[l].level) + GHOSTS)
      head[l].j -= (1 << head[l].level);
    if ( head[l].i < GHOSTS )
      head[l].i += (1 << head[l].level);
    if ( head[l].j < GHOSTS )
      head[l].j += (1 << head[l].level);
  }
  
  Point point = head[0];
  age[] = 1;
  point = head[1];
  age[] = -1;

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
  for (int l = 0; l < 2; l++) {
    if (head[l].i == target.i && head[l].j == target.j
	&& head[l].level == target.level) {
      targ.x = noise(), targ.y = noise();
      length[l] += 4;
      tot_time[l] += 10.*exp (-(tot_time[l] - 60)/40.);
    }
  }
}

event do_timer (i++) {
  perf.t = timer_elapsed (perf.gt);
  for (int l = 0; l < 2; l++) {
    seconds_left[l] = tot_time[l] - perf.t;
    if (seconds_left[l] <= 0) {
      stop_game = true;
      stop_screen();
      loser = l;
    }
  }
  putc('.', stdout);
  fflush (stdout);
}

event compute_s (i++, last) {
  foreach() {
    s[] = nodata;
    if (age[] > 0)
      s[] = 1;
    if (age[] == 1) 
      s[] = 10;
    if (age[] < 0)
      s[] = -1;
    if (age[] == -1)
      s[] = -10;
  }
  Point point = locate (targ.x, targ.y);
  s[] = 0;
}
