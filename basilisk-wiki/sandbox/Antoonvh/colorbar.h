/**
# `Colorbar.h`
   
A colourbar for bview in the bottom left corner. Mind the inconsistent spelling... 
 */

bool colourbar (Colormap map = jet , float size = 15, float pos[2] = {-.95, -.95},
		char * label = "", double lscale = 1.5, double min = -HUGE,
		double max = HUGE, bool horizontal = false, bool draw_box = false,
		bool mid = false, float lc[3] = {0}, float lw = 3, float fsize = 50) {

  //Learn from `draw_string` to deal with screen-space coordinates
  bview * view = draw();
  glDisable (GL_LIGHTING);
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();             
  glLoadIdentity();
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  float fheight = gl_StrokeHeight();
  if (!size)
    size = 15;
  int imax = 50;
  float width  = 2./size;
  float  h = 0, height  = 4*width, dh = height/(imax);
  glTranslatef (pos[0], pos[1], 0);
  
  // The colorbar itself
  double cmap [NCMAP][3];
  (*map) (cmap);
  assert (imax > 1);
  double lwn = 1;
  if (horizontal)
    lwn = view->width*height/(1.5*imax);
  else
    lwn = view->height*height/(1.5*imax);
  glLineWidth (lwn > 1 ? lwn : 1);
  glBegin(GL_LINES);
  for (int i = 0; i < imax; i++) {
    Color c = colormap_color (cmap, (float)i/(imax - 1), 0, 1); 
    glColor3f (c.r/255., c.g/255., c.b/255.);
    if (horizontal) {
      glVertex2f (h + dh/2, 0);
      glVertex2f (h + dh/2, width);
    } else {
      glVertex2f (0, h + dh/2);
      glVertex2f (width, h + dh/2);
    }
    h += dh;
  }
  glEnd();
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));
  glColor3f (lc[0], lc[1], lc[2]);
  
  // A box around the color scale
  if (draw_box) {
    glBegin (GL_LINE_LOOP);
    glVertex2f (0,0);
    if (horizontal) {
      glVertex2f (0, width);
      glVertex2f (height, width);
      glVertex2f (height, 0);
    } else {
      glVertex2f (width, 0);
      glVertex2f (width, height);
      glVertex2f (0, height);
    }
    glEnd();
  }

  // Min and max values when specified
  float fwidth  = gl_StrokeWidth ('1');
  if (!fsize)
    fsize = 20;
  float hscale = 2./(fsize*fwidth), vscale = hscale*view->width/view->height;
  char str[99];
  glColor3f (lc[0], lc[1], lc[2]);
  if (horizontal) 
    glTranslatef (0, -(fheight/(view->height)), 0);
  else 
    glTranslatef (width, -(fheight/(3*view->height)), 0);
  glScalef (hscale, vscale, 1.);
  sprintf (str, "%g\n", min);
  if (min > -HUGE) {
    if (horizontal)
      glTranslatef (-fwidth*(strlen(str) - 1)/2, 0, 0);
    gl_StrokeString (str);
    glTranslatef (0, fheight, 0);
    if (horizontal)
      glTranslatef (fwidth*(strlen(str) - 1)/2, 0, 0);
   
  }
  if (horizontal)
    glTranslatef (height/hscale,0, 0);
  else
    glTranslatef (0, height/vscale, 0);
  sprintf (str, "%g\n", max);
  if (max < HUGE) {
    if (horizontal)
      glTranslatef (-fwidth*(strlen(str) - 1)/2, 0, 0);
    gl_StrokeString (str);
    glTranslatef (0, fheight, 0);
    if (horizontal)
      glTranslatef (fwidth*(strlen(str) - 1)/2, 0, 0);
   
  }
  // Add central value
  if (mid) {
    sprintf (str, "%g\n", (min + max)/2);
    if (horizontal) 
      glTranslatef (-height/(2*hscale) - fwidth*(strlen(str) - 1)/2,0, 0);
    else
      glTranslatef (0, -height/(2*vscale), 0);
    gl_StrokeString (str);
    glTranslatef (0, fheight, 0);
    if (horizontal)
      glTranslatef (height/(2*hscale) + fwidth*(strlen(str) - 1)/2, 0, 0);
    else
      glTranslatef (0, height/(2*vscale), 0);
  }
  // Add label
  if (horizontal)
    glTranslatef (-height/(2*hscale) - lscale*fwidth*(strlen(label) - 1)/2, width/vscale, 0);
  else
    glTranslatef (-width/hscale, 0, 0);
  
  glScalef (lscale, lscale, 1.);
  glTranslatef (0, fheight, 0);
  gl_StrokeString (label);
  
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();  
  return true;
}
