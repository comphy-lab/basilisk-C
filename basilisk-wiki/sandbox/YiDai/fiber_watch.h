/**
# A fiber bview function

Similar to Antoon's [tracer-particles.h](/sandbox/Antoonvh/tracer-particles.h). */

trace 
void fiberW_node (lagfiber * fiber, float s = 4, float pc[3] = {0,0,0})
{
  bview * view = draw();
  #if dimension == 2
  glTranslatef (0., 0., view->lc*view->fov/24.); //from draw_lines()
  #endif
  glColor3f (pc[0], pc[1], pc[2]);
  glPointSize (s*view->samples);
  glBegin (GL_POINTS);
  for (int ii = 0; ii < fiber->nlp; ii++){
    #if dimension == 2
    glvertex2d (view, fiber->nodes[ii].post1.x, fiber->nodes[ii].post1.y);
    #else // dimension == 3
    glvertex3d (view, fiber->nodes[ii].post1.x, fiber->nodes[ii].post1.y, fiber->nodes[ii].post1.z);
    #endif
  }
  glEnd();
  view->ni++;
}
