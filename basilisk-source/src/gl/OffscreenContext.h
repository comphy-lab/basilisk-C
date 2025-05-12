#include <stdio.h>
#include <stdbool.h>

typedef struct _OffscreenContext OffscreenContext;

OffscreenContext *create_offscreen_context(int w, int h);
bool teardown_offscreen_context(OffscreenContext *ctx);
bool save_framebuffer(OffscreenContext *ctx, FILE * fp);
unsigned char * get_framebuffer_pixels (OffscreenContext *ctx);
unsigned int * get_framebuffer_depth (OffscreenContext *ctx);
