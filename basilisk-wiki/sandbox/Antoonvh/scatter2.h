/**
# A Scattering bview function

the `scatter` function can be used to scatter
[`Particles`](particle.h) in bview. The formulation is inspired by
other bview functions in [draw.h](src/draw.h).
 */

trace
void scatter (Particles p, float s = 3, float pc[3] = {0,0,0}) {
  bview * view = draw();
#if dimension == 2
  glTranslatef (0., 0., view->lc*view->fov/24.); //from draw_lines()
#endif
  glColor3f(pc[0], pc[1], pc[2]);
  glPointSize(s*view->samples);
  glBegin (GL_POINTS);
  foreach_particle_in(p) {
#if dimension == 2    
    glvertex2d (view, x, y);
#else // dimension == 3
    glvertex3d (view, x, y, z);
#endif
  }
  glEnd();
  view->ni++; 
}
/**
## Curve
From [A. Castillo's](../acastillo/filaments/filaments.h) sandbox.
 */

void draw_curve(Particles p, float lw = 5, float lc[3] = {0,0,0}, bool noloop = false){
  glColor3f(lc[0], lc[1], lc[2]);
  glLineWidth(lw);
  particles pp = pl[(int)p];
  foreach_particle_in(p) {
    glBegin(GL_LINES);
    int in = _j_particle < pn[p] - 1 ? _j_particle + 1 : noloop ? _j_particle : 0;
    glVertex3f(x, y, z);
    glVertex3f(pp[in].x, pp[in].y, pp[in].z);
    glEnd();
  }
}
/**
## Usage

* [Visualizing flow tracer particles](tracer-particles.h)
* [A confetti cannon](confetti.c)

*/
