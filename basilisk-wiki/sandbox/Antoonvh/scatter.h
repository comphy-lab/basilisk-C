/**
# A Scattering bview function

The formulation of the function is inspired by the functions under `draw.h` and assumes a `global int` by the name of `n_part` exists. It works in 2D and 3D.
 */
#if dimension == 3 
void glPointParameterfv(GLenum pname, const GLfloat * params);
#endif

/**
   The function is inspired by the functions under `draw.h` and assumes a `global int` by the name of `n_part` exists.
 */

trace
void scatter (coord * loc, float s = 20, float pc[3] = {0,0,0}, float coefs[3] = {0.01,0.2,0.5}){
  bview * view = draw();
#if dimension == 2
  glTranslatef (0., 0., view->lc*view->fov/24.); //from draw_lines()
  //  glPointParameterfv(GL_POINT_DISTANCE_ATTENUATION, coefs);
#endif
  //  glEnable (GL_BLEND);
  // glEnable (GL_POINT_SMOOTH);
  glColor3f(pc[0], pc[1], pc[2]);
  glPointSize(s);
  glBegin (GL_POINTS);
  for (int j = 0; j < n_part; j++){
#if dimension == 2    
    glvertex2d (view, loc[j].x, loc[j].y);
#else // dimension == 3
    glvertex3d (view, loc[j].x, loc[j].y, loc[j].z);
#endif
  }
  //glDisable (GL_BLEND);
  //glDisable (GL_POINT_SMOOTH);
  glEnd();
  view->ni++; 
}
/**
## Usage

* [Visualizing flow tracers](particles.h)

## Test

* [Distance attenuation in 3D](distanceat.c)

*/
