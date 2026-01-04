// Functions shared by OffscreenContext[platform].cc
// #include this directly after definition of struct OffscreenContext.

/*!
  Capture framebuffer from OpenGL and write it to the given ostream.
  Called by save_framebuffer() from platform-specific code.
*/
static bool save_framebuffer_common(OffscreenContext *ctx, FILE * output)
{
  if (!ctx) return false;
  int samplesPerPixel = 3; // R, G, B
  unsigned char pixels[ctx->width*ctx->height*samplesPerPixel];
  glReadPixels (0, 0, ctx->width, ctx->height,
		GL_RGB, GL_UNSIGNED_BYTE, pixels);
  fprintf (output, "P6 %d %d 255\n", ctx->width, ctx->height);
  fwrite (pixels, sizeof(unsigned char),
	  ctx->width*ctx->height*samplesPerPixel, output);
  return true;
}

//	Called by get_framebuffer_pixels() from platform-specific code.
static unsigned char * get_framebuffer_pixels_common (OffscreenContext *ctx)
{
  if (!ctx) return NULL;
  glReadPixels (0, 0, ctx->width, ctx->height,
		GL_RGBA, GL_UNSIGNED_BYTE, ctx->image);
  return ctx->image;
}

//	Called by get_framebuffer_depth() from platform-specific code.
static unsigned int * get_framebuffer_depth_common (OffscreenContext *ctx)
{
  if (!ctx) return NULL;

  if (!ctx->depth) {
    ctx->depth = malloc (sizeof(unsigned int)*ctx->width*ctx->height);
    glGenTextures (1, &ctx->depth_texture);
    glBindTexture (GL_TEXTURE_2D, ctx->depth_texture);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D (GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT,
		  ctx->width, ctx->height,
		  0, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, 0);
  }

  glCopyTexSubImage2D (GL_TEXTURE_2D, 0, 0, 0, 0, 0, ctx->width, ctx->height);
  glGetTexImage (GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT,
		 ctx->depth);
  
  return ctx->depth;
}

//	Called by create_offscreen_context() from platform-specific code.
static OffscreenContext *create_offscreen_context_common(OffscreenContext *ctx)
{
  if (!ctx) return NULL;
  GLenum err = glewInit(); // must come after Context creation and before FBO c$
  if (GLEW_OK != err) {
    fprintf (stderr, "Unable to init GLEW: %s\n", glewGetErrorString(err));
    return NULL;
  }

  ctx->fbo = fbo_new();
  if (!fbo_init(ctx->fbo, ctx->width, ctx->height)) {
    return NULL;
  }

  ctx->image = malloc (sizeof(unsigned char)*ctx->width*ctx->height*4);
  ctx->depth = NULL;
  
  return ctx;
}

//	Called by teardown_offscreen_context() from platform-specific code.
static bool teardown_offscreen_context_common(OffscreenContext *ctx)
{
  if (!ctx) return false;
  free (ctx->image);
  if (ctx->depth) {
    glDeleteTextures (1, &ctx->depth_texture);
    free (ctx->depth);
  }
  return true;
}
