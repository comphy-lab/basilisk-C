/* Adapted from
   
   https://github.com/openscad/openscad/blob/master/src/OffscreenContextGLX.cc

   There are also versions for Windows (OffscreenContextWGL.cc) and for
   Apple (OffscreenContextCGL.mm), but they haven't been adapted yet.

   by S. Popinet, 2017
*/

/*

  Create an OpenGL context without creating an OpenGL Window. for Linux.

  See also

  glxgears.c by Brian Paul from mesa-demos (mesa3d.org)
  http://cgit.freedesktop.org/mesa/demos/tree/src/xdemos?id=mesa-demos-8.0.1
  http://www.opengl.org/sdk/docs/man/xhtml/glXIntro.xml
  http://www.mesa3d.org/brianp/sig97/offscrn.htm
  http://glprogramming.com/blue/ch07.html
  OffscreenContext.mm (Mac OSX version)

*/

/*
 * Some portions of the code below are:
 * Copyright (C) 1999-2001  Brian Paul   All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * BRIAN PAUL BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
 * AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

#include "system-gl.h"
#include <GL/glx.h>

#include "fbo.h"
#include "OffscreenContext.h"

struct _OffscreenContext {
  GLXContext openGLContext;
  Display *xdisplay;
  Window xwindow;

  int width, height;
  fbo_t *fbo;
  unsigned char *image;
  unsigned int *depth;
  GLuint depth_texture;
};

#include "OffscreenContextAll.h"

static XErrorHandler original_xlib_handler = NULL;
static bool XCreateWindow_failed = false;
static int XCreateWindow_error(Display *dpy, XErrorEvent *event)
{
  fprintf (stderr, "XCreateWindow failed: XID: %d request: %d minor: %d\n",
	   (int) event->resourceid, event->request_code, event->minor_code);
  char description[1024];
  XGetErrorText( dpy, event->error_code, description, 1023 );
  fprintf (stderr, " error message: %s\n", description);
  XCreateWindow_failed = true;
  return 0;
}

/*
  create a dummy X window without showing it. (without 'mapping' it)
  and save information to the ctx.

  This purposely does not use glxCreateWindow, to avoid crashes,
  "failed to create drawable" errors, and Mesa "WARNING: Application calling 
  GLX 1.3 function when GLX 1.3 is not supported! This is an application bug!"

  This function will alter ctx.openGLContext and ctx.xwindow if successfull
*/
static bool create_glx_dummy_window (OffscreenContext * ctx)
{
  int attributes[] = {
    //support all 3, for OpenCSG
    GLX_DRAWABLE_TYPE, GLX_WINDOW_BIT | GLX_PIXMAP_BIT | GLX_PBUFFER_BIT,
    GLX_RENDER_TYPE,   GLX_RGBA_BIT,
    GLX_RED_SIZE, 8,
    GLX_GREEN_SIZE, 8,
    GLX_BLUE_SIZE, 8,
    GLX_ALPHA_SIZE, 8,
    GLX_DEPTH_SIZE, 24,
    GLX_STENCIL_SIZE, 0,
    GLX_DOUBLEBUFFER, true,
    None
  };

  Display * dpy = ctx->xdisplay;

  int num_returned = 0;
  GLXFBConfig * fbconfigs = glXChooseFBConfig( dpy, DefaultScreen(dpy),
					       attributes, &num_returned );
  if (fbconfigs == NULL) {
    fprintf (stderr, "glXChooseFBConfig failed\n");
    return false;
  }

   XVisualInfo * visinfo = glXGetVisualFromFBConfig( dpy, fbconfigs[0] );
  if (visinfo == NULL) {
    fprintf (stderr, "glXGetVisualFromFBConfig failed\n");
    XFree(fbconfigs);
    return false;
  }

  // can't depend on xWin==NULL at failure. use a custom Xlib error
  // handler instead.
  original_xlib_handler = XSetErrorHandler(XCreateWindow_error);

  Window root = DefaultRootWindow(dpy);
  XSetWindowAttributes xwin_attr;
  int width = ctx->width;
  int height = ctx->height;
  xwin_attr.background_pixmap = None;
  xwin_attr.background_pixel = 0;
  xwin_attr.border_pixel = 0;
  xwin_attr.colormap = XCreateColormap( dpy, root, visinfo->visual, AllocNone);
  xwin_attr.event_mask = StructureNotifyMask | ExposureMask | KeyPressMask;
  unsigned long int mask = CWBackPixel | CWBorderPixel |
    CWColormap | CWEventMask;

  Window xWin = XCreateWindow( dpy, root, 0, 0, width, height,
			       0, visinfo->depth, InputOutput,
			       visinfo->visual, mask, &xwin_attr );

  // Window xWin = XCreateSimpleWindow( dpy, DefaultRootWindow(dpy),
  // 0,0,42,42, 0,0,0 );

  XSync(dpy, false);
  if (XCreateWindow_failed) {
    XFree(visinfo);
    XFree(fbconfigs);
    return false;
  }
  XSetErrorHandler(original_xlib_handler);

  // Most programs would call XMapWindow here. But we don't, to keep
  // the window hidden
  // XMapWindow( dpy, xWin );

  GLXContext context = glXCreateNewContext(dpy, fbconfigs[0],
					   GLX_RGBA_TYPE, NULL, true);
  if (context == NULL) {
    fprintf (stderr, "glXCreateNewContext failed\n");
    XDestroyWindow(dpy, xWin);
    XFree(visinfo);
    XFree(fbconfigs);
    return false;
  }

  //GLXWindow glxWin = glXCreateWindow( dpy, fbconfigs[0], xWin, NULL );

  if (!glXMakeContextCurrent( dpy, xWin, xWin, context )) {
    //if (!glXMakeContextCurrent( dpy, glxWin, glxWin, context )) {
    fprintf (stderr, "glXMakeContextCurrent failed\n");
    glXDestroyContext(dpy, context);
    XDestroyWindow(dpy, xWin);
    XFree(visinfo);
    XFree(fbconfigs);
    return false;
  }

  ctx->openGLContext = context;
  ctx->xwindow = xWin;

  XFree(visinfo);
  XFree(fbconfigs);

  return true;
}

static bool create_glx_dummy_context (OffscreenContext * ctx);

OffscreenContext * create_offscreen_context(int w, int h)
{
  OffscreenContext * ctx = calloc (1, sizeof (OffscreenContext));
  ctx->width = w, ctx->height = h;

  // before an FBO can be setup, a GLX context must be created
  // this call alters ctx->xDisplay and ctx->openGLContext 
  // and ctx->xwindow if successful
  if (!create_glx_dummy_context(ctx)) {
    free (ctx);
    return NULL;
  }
  
  return create_offscreen_context_common (ctx);
}

bool teardown_offscreen_context(OffscreenContext *ctx)
{
  if (ctx) {
    if (!teardown_offscreen_context_common(ctx))
      return false;
    fbo_unbind(ctx->fbo);
    fbo_delete(ctx->fbo);
    XDestroyWindow( ctx->xdisplay, ctx->xwindow );
    glXDestroyContext( ctx->xdisplay, ctx->openGLContext );
    XCloseDisplay( ctx->xdisplay );
    free (ctx);
    return true;
  }
  return false;
}

bool save_framebuffer(OffscreenContext * ctx, FILE * output)
{
  glXSwapBuffers(ctx->xdisplay, ctx->xwindow);
  return save_framebuffer_common(ctx, output);
}

unsigned char * get_framebuffer_pixels (OffscreenContext * ctx)
{
  glXSwapBuffers(ctx->xdisplay, ctx->xwindow);
  return get_framebuffer_pixels_common(ctx);
}

unsigned int * get_framebuffer_depth (OffscreenContext * ctx)
{
  if (!ctx) return NULL;
  return get_framebuffer_depth_common(ctx);
}

#pragma GCC diagnostic ignored "-Waddress"
static bool create_glx_dummy_context(OffscreenContext * ctx)
{
  // This will alter ctx->openGLContext and ctx->xdisplay and
  // ctx->xwindow if successfull
  int major;
  int minor;
  int result = false;

  ctx->xdisplay = XOpenDisplay(NULL);
  if (ctx->xdisplay == NULL) {
    fprintf (stderr, "Unable to open a connection to the X server.\n");
    char * dpyenv = getenv ("DISPLAY");
    fprintf (stderr, "DISPLAY=%s\n", dpyenv ? dpyenv : "");
    return false;
  }

  // glxQueryVersion is not always reliable. Use it, but then
  // also check to see if GLX 1.3 functions exist

  glXQueryVersion(ctx->xdisplay, &major, &minor);
  if (major==1 && minor<=2 && glXGetVisualFromFBConfig==NULL) {
    fprintf (stderr, "Error: GLX version 1.3 functions missing. "
	     "Your GLX version: %d.%d\n", major, minor);
  } else {
    result = create_glx_dummy_window(ctx);
  }

  if (!result) XCloseDisplay(ctx->xdisplay);
  return result;
}

#if TEST

#include <GL/glut.h>

int main(int argc, char * argv[]) {

  glutInit(&argc, argv);
  
  OffscreenContext * ctx = create_offscreen_context(800, 600);

  ///////////////////////////

  glClearColor( 0.0, 0.0, 0.0, 1. );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glEnable( GL_DEPTH_TEST );
  glShadeModel( GL_FLAT );
  glViewport( 0, 0, ctx->width, ctx->height );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity( );
  gluPerspective( 90., 1.,0.1, 1000. );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
  gluLookAt( 0., 0., 3., 0., 0., 0., 0., 1., 0. );
  //  glTranslatef( TransXYZ[0], TransXYZ[1], TransXYZ[2] );
  //  glMultMatrixf( RotMatrix );
  glRotatef (45, 1, 0, 0);
  //  glScalef( scale, scale, scale );
  glColor3f( 1., 1., 1. );
  glutWireTeapot( 1. );
  
  ////////////////

  assert (save_framebuffer (ctx, stdout));
  
  assert (teardown_offscreen_context (ctx));
  return 0;
}

#endif
