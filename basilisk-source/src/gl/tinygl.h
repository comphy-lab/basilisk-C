/**
# Minimal OpenGL header

See also [fb_tiny.c]().

It would be handy if the following constants were the same for all
OpenGL implementations (they are copied from the Mesa
implementation). */

#define GL_MODELVIEW				0x1700
#define GL_MODELVIEW_MATRIX			0x0BA6
#define GL_NO_ERROR 				0
#define GL_CULL_FACE				0x0B44
#define GL_CULL_FACE_MODE			0x0B45
#define GL_DEPTH_TEST				0x0B71
#define GL_NORMALIZE				0x0BA1
#define GL_DITHER				0x0BD0
#define GL_SMOOTH				0x1D01
#define GL_SMOOTH_POINT_SIZE_RANGE		0x0B12
#define GL_SMOOTH_POINT_SIZE_GRANULARITY	0x0B13
#define GL_SMOOTH_LINE_WIDTH_RANGE		0x0B22
#define GL_SMOOTH_LINE_WIDTH_GRANULARITY	0x0B23
#define GL_PERSPECTIVE_CORRECTION_HINT		0x0C50
#define GL_FASTEST				0x1101
#define GL_POLYGON_SMOOTH_HINT			0x0C53
#define GL_LIGHT_MODEL_TWO_SIDE			0x0B52
#define GL_TRUE					1
#define GL_POSITION				0x1203
#define GL_DIFFUSE				0x1201
#define GL_LIGHT0                               0x4000
#define GL_LIGHTING				0x0B50
#define GL_LIGHTING_BIT				0x00000040
#define GL_FRONT				0x0404
#define GL_BACK					0x0405
#define GL_FRONT_AND_BACK			0x0408
#define GL_AMBIENT				0x1200
#define GL_AMBIENT_AND_DIFFUSE			0x1602
#define GL_COLOR_MATERIAL			0x0B57
#define GL_COLOR_MATERIAL_FACE			0x0B55
#define GL_COLOR_MATERIAL_PARAMETER		0x0B56
#define GL_PROJECTION				0x1701
#define GL_PROJECTION_MATRIX			0x0BA7
#define GL_QUADS				0x0007
#define GL_INVALID_ENUM				0x0500
#define GL_INVALID_VALUE			0x0501
#define GL_INVALID_OPERATION			0x0502
#define GL_INVALID_FRAMEBUFFER_OPERATION        0x0506
#define GL_OUT_OF_MEMORY			0x0505
#define GL_STACK_UNDERFLOW			0x0504
#define GL_STACK_OVERFLOW			0x0503
#define GL_TEXTURE_2D                           0x0DE1
#define GL_VIEWPORT				0x0BA2
#define GL_VIEWPORT_BIT				0x00000800
#define GL_PROJECTION_MATRIX			0x0BA7
#define GL_COLOR_BUFFER_BIT			0x00004000
#define GL_DEPTH_BUFFER_BIT			0x00000100
#define GL_LINE_LOOP				0x0002
#define GL_LINES				0x0001
#define GL_LINE_STRIP				0x0003
#define GL_POLYGON				0x0009
#define GL_POLYGON_MODE				0x0B40
#define GL_POLYGON_SMOOTH			0x0B41
#define GL_POLYGON_STIPPLE			0x0B42
#define GL_POLYGON_OFFSET_FACTOR		0x8038
#define GL_POLYGON_OFFSET_UNITS			0x2A00
#define GL_POLYGON_OFFSET_POINT			0x2A01
#define GL_POLYGON_OFFSET_LINE			0x2A02
#define GL_POLYGON_OFFSET_FILL			0x8037
#define GL_POLYGON_TOKEN			0x0703
#define GL_POLYGON_SMOOTH_HINT			0x0C53
#define GL_POLYGON_BIT				0x00000008
#define GL_POLYGON_STIPPLE_BIT			0x00000010
#define GL_TRIANGLE_FAN				0x0006
#define GL_LINE_LOOP				0x0002
#define GL_LINES				0x0001
#define GL_LINE_STRIP				0x0003
#define GL_POLYGON				0x0009
#define GL_POLYGON_MODE				0x0B40
#define GL_POLYGON_SMOOTH			0x0B41
#define GL_POLYGON_STIPPLE			0x0B42
#define GL_POLYGON_OFFSET_FACTOR		0x8038
#define GL_POLYGON_OFFSET_UNITS			0x2A00
#define GL_POLYGON_OFFSET_POINT			0x2A01
#define GL_POLYGON_OFFSET_LINE			0x2A02
#define GL_POLYGON_OFFSET_FILL			0x8037
#define GL_POLYGON_TOKEN			0x0703
#define GL_POLYGON_SMOOTH_HINT			0x0C53
#define GL_POLYGON_BIT				0x00000008
#define GL_POLYGON_STIPPLE_BIT			0x00000010
#define GL_TRIANGLE_FAN				0x0006
#define GL_RGB					0x1907
#define GL_RGBA					0x1908
#define GL_RGBA_MODE				0x0C31
#define GL_RGB4					0x804F
#define GL_RGB5					0x8050
#define GL_RGB8					0x8051
#define GL_RGB10				0x8052
#define GL_RGB12				0x8053
#define GL_RGB16				0x8054
#define GL_RGBA2				0x8055
#define GL_RGBA4				0x8056
#define GL_RGB5_A1				0x8057
#define GL_RGBA8				0x8058
#define GL_RGB10_A2				0x8059
#define GL_RGBA12				0x805A
#define GL_RGBA16				0x805B
#define GL_RGB_SCALE				0x8573
#define GL_FLOAT				0x1406
#define GL_TEXTURE_MIN_FILTER			0x2801
#define GL_LINEAR_ATTENUATION			0x1208
#define GL_LINEAR				0x2601
#define GL_LINEAR_MIPMAP_NEAREST		0x2701
#define GL_LINEAR_MIPMAP_LINEAR			0x2703
#define GL_TEXTURE_MAG_FILTER			0x2800
#define GL_TEXTURE_WRAP_S			0x2802
#define GL_CLAMP_TO_EDGE			0x812F
#define GL_TEXTURE_1D                           0x0DE0
#define GL_TEXTURE_WRAP_T			0x2803
#define GL_POINTS                               0x0000


typedef unsigned int    GLenum;
typedef int             GLint;
typedef unsigned int    GLuint;
typedef float           GLfloat;
typedef double          GLdouble;
typedef int             GLsizei;
typedef unsigned int	GLbitfield;
typedef unsigned char   GLubyte;

void glBegin (GLenum mode);
void glEnd (void);
void glClear (GLbitfield mask);
void glClearColor (float red, float green, float blue, float alpha);
void glBindTexture (GLenum target, GLuint texture);
void glColor3f (GLfloat red, GLfloat green, GLfloat blue);
void glColorMaterial (GLenum face, GLenum mode);
void glDisable (GLenum cap);
void glEnable (GLenum cap);
void glFinish (void);
void glGetDoublev (GLenum pname, GLdouble * params);
void glGetIntegerv (GLenum pname, GLint * params);
void glHint (GLenum target, GLenum mode);
void glLightfv (GLenum light, GLenum pname, const GLfloat *params);
void glGetLightfv (GLenum light, GLenum pname, GLfloat * params);
void glLightModeli (GLenum pname, GLint param);
void glLineWidth (GLfloat width);
void glPointSize (GLfloat size);
void glNormal3d (GLdouble nx, GLdouble ny, GLdouble nz);
void glOrtho (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top,
	      GLdouble nearVal, GLdouble farVal);
void glShadeModel (GLenum mode);
void glTexCoord1d (GLdouble s);
void glTexCoord2f (GLfloat s, GLfloat t);
void glTexImage1D (GLenum target, GLint level, GLint internalFormat,
		   GLsizei width, GLint border, GLenum format, GLenum type,
		   const void * data);
void glTexParameteri (GLenum target, GLenum pname, GLint param);
void glVertex3d (GLdouble x, GLdouble y, GLdouble z);
void glVertex3f (GLfloat x, GLfloat y, GLfloat z);

GLenum glGetError (void);
void glGetFloatv (GLenum pname, GLfloat * params);
void glMultMatrixf (const GLfloat * m);
void glLoadIdentity (void);
void glScalef (GLfloat x, GLfloat y, GLfloat z);
void glTranslatef (GLfloat x, GLfloat y, GLfloat z);
void glRotatef (GLfloat angle, GLfloat x, GLfloat y, GLfloat z);
void glMatrixMode (GLenum mode);
void glPopMatrix (void);
void glPushMatrix (void);
void glLoadMatrixd (const GLdouble * m);
