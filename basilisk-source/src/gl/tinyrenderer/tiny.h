typedef struct { // fixme: should be opaque
  unsigned char * image;
  int width, height;
  real * zbuffer;
  float * depth;
} framebuffer;

framebuffer * framebuffer_new (unsigned width, unsigned height);
void framebuffer_destroy (framebuffer * p);
unsigned char * framebuffer_image (framebuffer * p);
float * framebuffer_depth (framebuffer * p);

extern framebuffer * TinyFramebuffer;

void tiny_viewport (const int x, const int y, const int w, const int h);

typedef struct {
  unsigned char r,g,b,a;
} TinyColor;

typedef int (* TinyShader) (const void * shader, const vec3 bar, TinyColor * color);

void tiny_triangle (const vec4 clip_verts[3],
		    const void * shader, const TinyShader fragment,
		    int face, // -1: back, 0: front and back, 1: front
		    framebuffer * image);
void tiny_line (const vec4 clip_verts0, const vec4 clip_verts1,
		const TinyColor * color, float thickness,
		framebuffer * image);
void tiny_point (const vec4 clip_verts0, const TinyColor * color, float raidus,
		 framebuffer * image);
