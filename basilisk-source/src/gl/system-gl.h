#ifndef NULLGL
#include <GL/glew.h>

#ifdef __APPLE__
 #include <OpenGL/OpenGL.h>
#else
 #include <GL/gl.h>
 #include <GL/glu.h>
#endif

#else // NULLGL
#define GLint int
#define GLuint unsigned int
inline void glColor4fv( float *c ) {}
#endif // NULLGL

bool report_glerror(const char * function);
