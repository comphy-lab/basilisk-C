#include <limits>
#include <iostream>
#include "geometry.h"
#include "model.h"
extern "C" {
#include "tiny.h"
}

mat4 Projection, ModelView;

void projection(const double f) { // check https://en.wikipedia.org/wiki/Camera_matrix
  Projection = (mat4){{1,0,0,0}, {0,-1,0,0}, {0,0,1,0}, {0,0,-1/f,0}};
}

// check https://github.com/ssloy/tinyrenderer/wiki/Lesson-5-Moving-the-camera
void lookat(const vec3 eye, const vec3 center, const vec3 up) {
  vec3 z = vec3_normalized (vec3_sub (center, eye));
  vec3 x = vec3_normalized (cross(up, z));
  vec3 y = vec3_normalized (cross(z, x));
  mat4 Minv = (mat4){{x.x,x.y,x.z,0},   {y.x,y.y,y.z,0},   {z.x,z.y,z.z,0},   {0,0,0,1}};
  mat4 Tr   = (mat4){{1,0,0,-eye.x}, {0,1,0,-eye.y}, {0,0,1,-eye.z}, {0,0,0,1}};
  ModelView = mat4_mul4 (Minv, Tr);
}

constexpr int width  = 1600; // output image size
constexpr int height = 1600;
constexpr vec3 light_dir = {1,1,1}; // light source
constexpr vec3       eye = {1,1,3}; // camera position
constexpr vec3    center = {0,0,0}; // camera direction
constexpr vec3        up = {0,1,0}; // camera up vector
const TinyColor white = {255,255,255,255}, red = {0,0,255,255}, black = {0,0,0,255};

struct Shader {
  const Model &model;
  vec3 uniform_l;       // light direction in view coordinates
  mat23 varying_uv;  // triangle uv coordinates, written by the vertex shader, read by the fragment shader
  mat3 varying_nrm; // normal per vertex to be interpolated by FS
  mat3 view_tri;    // triangle in view coordinates

  Shader(const Model &m) : model(m) {
    // transform the light vector to view coordinates
    uniform_l = vec3_normalized (vec4_proj3 (mat4_mul (ModelView, vec4_embed(light_dir, 0.))));
  }

  static TGAColor sample2D(const TGAImage &img, vec2 &uvf) {
    return img.get(uvf.x*img.width(), uvf.y*img.height());
  }

  virtual void vertex(const int iface, const int nthvert, vec4& gl_Position) {
    set_col23 (&varying_uv, nthvert, model.uv (iface, nthvert));
    set_col3 (&varying_nrm, nthvert,
	      vec4_proj3(mat4_mul (mat4_invert_transpose (ModelView),
				   vec4_embed(model.normal(iface, nthvert), 0.))));
    gl_Position = mat4_mul (ModelView, vec4_embed (model.vert(iface, nthvert), 1.));
    set_col3 (&view_tri, nthvert, vec4_proj3 (gl_Position));
    gl_Position = mat4_mul (Projection, gl_Position);
  }
};

static int fragment (const void * shader, const vec3 bar, TinyColor * gl_FragColor)
{
  Shader * s = (Shader *) shader;
  vec3 bn = vec3_normalized (mat3_mul (s->varying_nrm, bar)); // per-vertex normal interpolation
  vec2 uv = mat23_mul (s->varying_uv, bar); // tex coord interpolation

  // for the math refer to the tangent space normal mapping lecture
  // https://github.com/ssloy/tinyrenderer/wiki/Lesson-6bis-tangent-space-normal-mapping
  mat3 AI = mat3_invert ({
      vec3_sub (mat3_col (s->view_tri, 1), mat3_col (s->view_tri, 0)),
      vec3_sub (mat3_col (s->view_tri, 2), mat3_col (s->view_tri, 0)),
      bn
    });
  vec3 i = mat3_mul (AI, {s->varying_uv.x.y - s->varying_uv.x.x, s->varying_uv.x.z - s->varying_uv.x.x, 0});
  vec3 j = mat3_mul (AI, {s->varying_uv.y.y - s->varying_uv.y.x, s->varying_uv.y.z - s->varying_uv.y.x, 0});
  mat3 B = mat3_transpose ({ vec3_normalized (i), vec3_normalized (j), bn });

  // transform the normal from the texture to the tangent space
  vec3 n = vec3_normalized (mat3_mul (B, s->model.normal(uv)));
  double diff = std::max(0., vec3_scalar (n, s->uniform_l)); // diffuse light intensity
  // reflected light direction, specular mapping is described here:
  // https://github.com/ssloy/tinyrenderer/wiki/Lesson-6-Shaders-for-the-software-renderer
  vec3 r = vec3_normalized (vec3_sub (vec3_mul (n, vec3_scalar (n, s->uniform_l)*2.), s->uniform_l));
  // specular intensity, note that the camera lies on the z-axis
  // (in view), therefore simple -r.z
  double spec = std::pow(std::max(-r.z, 0.), 5 + s->sample2D(s->model.specular(), uv)[0]);
    
  TGAColor c = s->sample2D(s->model.diffuse(), uv);
  for (int i : {0,1,2}) // (a bit of ambient light, diff + spec), clamp the result
    ((unsigned char *)gl_FragColor)[i] = std::min<int>(10 + c[i]*(diff + spec), 255);
  return false; // the pixel is not discarded
}

int main (int argc, char** argv) {
  if (2 > argc) {
    std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
    return 1;
  }
  framebuffer * image = framebuffer_new (width, height);
  lookat (eye, center, up);                            // build the ModelView matrix
  tiny_viewport (width/8, height/8, width*3/4, height*3/4); // build the Viewport matrix
  projection (vec3_norm (vec3_sub (eye, center)));     // build the Projection matrix

  for (int m = 1; m < argc; m++) { // iterate through all input objects
    Model model(argv[m]);
    Shader shader(model);
    for (int i = 0; i < model.nfaces(); i++) { // for every triangle
      vec4 clip_vert[3]; // triangle coordinates (clip coordinates), written by VS, read by FS
      for (int j : {0,1,2})
	shader.vertex(i, j, clip_vert[j]); // call the vertex shader for each triangle vertex
#if 1
      float thickness = 1.;
      tiny_line (clip_vert[0], clip_vert[1], &red, thickness, image);
      tiny_line (clip_vert[1], clip_vert[2], &red, thickness, image);
      tiny_line (clip_vert[2], clip_vert[0], &red, thickness, image);
#endif
#if 1	    
      tiny_triangle (clip_vert, &shader, fragment, 1, image); // actual rasterization routine call
#endif
    }
  }
#if 1 
  TGAImage tga(width, height, TGAImage::RGB);
  unsigned char * p = framebuffer_image (image);
  for (int x = 0; x < width; x++)
    for (int y = 0; y < height; y++) {
      unsigned char * i = p + 4*(x + y*width);
      tga.set (x, y, TGAColor{i[0], i[1], i[2], i[3]});
    }
  tga.write_tga_file("framebuffer.tga");
#endif
  framebuffer_destroy (image);  
  return 0;
}
