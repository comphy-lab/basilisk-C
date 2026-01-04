#include <grid/gpu/glad.h>
#include <GLFW/glfw3.h>
#if DEBUG_OPENGL
#include <grid/gpu/debug.h>
#endif

#pragma autolink -L$BASILISK/grid/gpu -lglfw -lgpu -ldl

static struct {
  ///// GPU /////
  GLFWwindow * window;
  GLuint * ssbo;
  bool fragment_shader;
  int current_shader, nssbo;
  size_t max_ssbo_size, current_size;
} GPUContext = {
  .current_shader = -1,
};

static void gpu_check_error (const char * stmt, const char * fname, int line)
{
  GLenum err = glGetError();
  if (err != GL_NO_ERROR) {
    fprintf (stderr, "%s:%d: error: OpenGl %08x for '%s;'\n",
	     fname, line, err, stmt);
    abort();
  }
}

#ifdef NDEBUG
// helper macro that checks for GL errors.
#define GL_C(stmt) do {	stmt; } while (0)
#else
// helper macro that checks for GL errors.
#define GL_C(stmt) do {							\
    stmt;								\
    gpu_check_error (#stmt, __FILE__, LINENO);				\
  } while (0)
#endif

void gpu_free_solver (void)
{
  GL_C (glFinish());
  GL_C (glBindFramebuffer (GL_FRAMEBUFFER, 0));
  glDeleteBuffers (GPUContext.nssbo, GPUContext.ssbo);
  free (GPUContext.ssbo);
  glfwTerminate();
  GPUContext.window = NULL;
}

static char * getShaderLogInfo (GLuint shader)
{
  char * infoLog = NULL;
  GLint len;
  GL_C (glGetShaderiv (shader, GL_INFO_LOG_LENGTH, &len));
  if (len > 0) {
    GLsizei actualLen;
    infoLog = malloc (len);
    GL_C (glGetShaderInfoLog (shader, len, &actualLen, infoLog));
  }
  return infoLog;
}

static char * getProgramLogInfo (GLuint program)
{
  char * infoLog = NULL;
  GLint len;
  GL_C (glGetProgramiv (program, GL_INFO_LOG_LENGTH, &len));
  if (len > 0) {
    GLsizei actualLen;
    infoLog = malloc (len);
    GL_C (glGetProgramInfoLog (program, len, &actualLen, infoLog));
  }
  return infoLog;
}

char * gpu_errors (const char * errors, const char * source, char * fout);

static GLuint createShaderFromString (const char * shaderSource,
				      const GLenum shaderType)
{
  GLuint shader;

  GL_C (shader = glCreateShader(shaderType));
  GL_C (glShaderSource (shader, 1, &shaderSource, NULL));
  GL_C (glCompileShader (shader));

  GLint compileStatus;
  GL_C (glGetShaderiv (shader, GL_COMPILE_STATUS, &compileStatus));
  if (compileStatus != GL_TRUE) {
    char * info = getShaderLogInfo (shader);
#if PRINTSHADERERROR
    fputs (shaderSource, stderr);
    fputs (info, stderr);
#endif
    char * error = gpu_errors (info, shaderSource, NULL);
    fputs (error, stderr);
    sysfree (error);
    free (info);
    glDeleteShader (shader);
    return 0;
  }

  return shader;
}

static GLuint loadNormalShader (const char * vsSource, const char * fsShader)
{
  GLuint vs = 0;
  if (vsSource) {
    vs = createShaderFromString (vsSource, GL_VERTEX_SHADER);
    if (!vs)
      return 0;
  }
  GLuint fs = createShaderFromString (fsShader, vsSource ? GL_FRAGMENT_SHADER : GL_COMPUTE_SHADER);
  if (!fs)
    return 0;
  
  GLuint shader = glCreateProgram();
  if (vs)
    glAttachShader (shader, vs);
  glAttachShader (shader, fs);
  glLinkProgram (shader);

  GLint Result;
  glGetProgramiv (shader, GL_LINK_STATUS, &Result);

  if (Result == GL_FALSE) {
    char * info = getProgramLogInfo (shader);
    fprintf (stderr, "GLSL: could not link shader \n\n%s\n%s\n%s\n",
	     info, vsSource, fsShader);
    free (info);
    glDeleteProgram (shader);
    shader = 0;
  }

  if (shader) {
    if (vs)
      glDetachShader (shader, vs);
    glDetachShader (shader, fs);
  }

  if (vs)
    glDeleteShader (vs);
  glDeleteShader (fs);
  
  return shader;
}

typedef struct {
  char * s;
  int index;
} GLString;

GLString gpu_limits_list[] = {
  {"GL_MAX_DRAW_BUFFERS", GL_MAX_DRAW_BUFFERS},
  {"GL_MAX_VERTEX_UNIFORM_COMPONENTS", GL_MAX_VERTEX_UNIFORM_COMPONENTS},
  {"GL_MAX_VERTEX_UNIFORM_BLOCKS", GL_MAX_VERTEX_UNIFORM_BLOCKS},
  {"GL_MAX_VERTEX_OUTPUT_COMPONENTS", GL_MAX_VERTEX_OUTPUT_COMPONENTS},
  {"GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS", GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_VERTEX_OUTPUT_COMPONENTS", GL_MAX_VERTEX_OUTPUT_COMPONENTS},
#if 1  
  {"GL_MAX_TESS_CONTROL_UNIFORM_COMPONENTS", GL_MAX_TESS_CONTROL_UNIFORM_COMPONENTS},
  {"GL_MAX_TESS_CONTROL_UNIFORM_BLOCKS", GL_MAX_TESS_CONTROL_UNIFORM_BLOCKS},
  {"GL_MAX_TESS_CONTROL_INPUT_COMPONENTS", GL_MAX_TESS_CONTROL_INPUT_COMPONENTS},
  {"GL_MAX_TESS_CONTROL_OUTPUT_COMPONENTS", GL_MAX_TESS_CONTROL_OUTPUT_COMPONENTS},
  {"GL_MAX_TESS_CONTROL_TEXTURE_IMAGE_UNITS", GL_MAX_TESS_CONTROL_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_TESS_CONTROL_OUTPUT_COMPONENTS", GL_MAX_TESS_CONTROL_OUTPUT_COMPONENTS},
  {"GL_MAX_TESS_EVALUATION_UNIFORM_COMPONENTS", GL_MAX_TESS_EVALUATION_UNIFORM_COMPONENTS},
  {"GL_MAX_TESS_EVALUATION_UNIFORM_BLOCKS", GL_MAX_TESS_EVALUATION_UNIFORM_BLOCKS},
  {"GL_MAX_TESS_EVALUATION_INPUT_COMPONENTS", GL_MAX_TESS_EVALUATION_INPUT_COMPONENTS},
  {"GL_MAX_TESS_EVALUATION_OUTPUT_COMPONENTS", GL_MAX_TESS_EVALUATION_OUTPUT_COMPONENTS},
  {"GL_MAX_TESS_EVALUATION_TEXTURE_IMAGE_UNITS", GL_MAX_TESS_EVALUATION_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_TESS_EVALUATION_OUTPUT_COMPONENTS", GL_MAX_TESS_EVALUATION_OUTPUT_COMPONENTS},
#endif
#if 1  
  {"GL_MAX_COMPUTE_UNIFORM_COMPONENTS", GL_MAX_COMPUTE_UNIFORM_COMPONENTS},
  {"GL_MAX_COMPUTE_UNIFORM_BLOCKS", GL_MAX_COMPUTE_UNIFORM_BLOCKS},
  {"GL_MAX_COMPUTE_TEXTURE_IMAGE_UNITS", GL_MAX_COMPUTE_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_SHADER_STORAGE_BLOCK_SIZE", GL_MAX_SHADER_STORAGE_BLOCK_SIZE},
  {"GL_MAX_COMPUTE_SHADER_STORAGE_BLOCKS", GL_MAX_COMPUTE_SHADER_STORAGE_BLOCKS},
  {"GL_MAX_TEXTURE_BUFFER_SIZE", GL_MAX_TEXTURE_BUFFER_SIZE},
#endif
  {"GL_MAX_GEOMETRY_UNIFORM_COMPONENTS", GL_MAX_GEOMETRY_UNIFORM_COMPONENTS},
  {"GL_MAX_GEOMETRY_UNIFORM_BLOCKS", GL_MAX_GEOMETRY_UNIFORM_BLOCKS},
  {"GL_MAX_GEOMETRY_INPUT_COMPONENTS", GL_MAX_GEOMETRY_INPUT_COMPONENTS},
  {"GL_MAX_GEOMETRY_OUTPUT_COMPONENTS", GL_MAX_GEOMETRY_OUTPUT_COMPONENTS},
  {"GL_MAX_GEOMETRY_TEXTURE_IMAGE_UNITS", GL_MAX_GEOMETRY_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_GEOMETRY_OUTPUT_COMPONENTS", GL_MAX_GEOMETRY_OUTPUT_COMPONENTS},
  {"GL_MAX_FRAGMENT_UNIFORM_COMPONENTS", GL_MAX_FRAGMENT_UNIFORM_COMPONENTS},
  {"GL_MAX_FRAGMENT_UNIFORM_BLOCKS", GL_MAX_FRAGMENT_UNIFORM_BLOCKS},
  {"GL_MAX_FRAGMENT_INPUT_COMPONENTS", GL_MAX_FRAGMENT_INPUT_COMPONENTS},
  {"GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS", GL_MAX_COMBINED_TEXTURE_IMAGE_UNITS},
  {"GL_MAX_FRAGMENT_IMAGE_UNIFORMS", GL_MAX_FRAGMENT_IMAGE_UNIFORMS},
  {"GL_MAX_IMAGE_UNITS", GL_MAX_IMAGE_UNITS},
  {NULL}
};

void printWorkGroupsCapabilities()
{
  int workgroup_size[3];
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 0, &workgroup_size[0]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 1, &workgroup_size[1]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 2, &workgroup_size[2]);
  printf ("Maximal workgroup sizes:\n\tx:%u\n\ty:%u\n\tz:%u\n",
	  workgroup_size[0], workgroup_size[1], workgroup_size[2]);

  int workgroup_count[3];
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &workgroup_count[0]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 1, &workgroup_count[1]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 2, &workgroup_count[2]);
  printf ("Maximum number of local invocations:\n\tx:%u\n\ty:%u\n\tz:%u\n",
	  workgroup_count[0], workgroup_count[1], workgroup_count[2]);

  int workgroup_invocations;
  glGetIntegerv (GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS, &workgroup_invocations);
  printf ("Maximum workgroup invocations:\n\t%u\n", workgroup_invocations);
}

#if TRACE == 3
# define tracing_foreach(name, file, line) tracing(name, file, line)
# define end_tracing_foreach(name, file, line) end_tracing(name, file, line)
#else
# define tracing_foreach(name, file, line)
# define end_tracing_foreach(name, file, line)
#endif

macro2 BEGIN_FOREACH()
{
  if (_gpu_done_)
    _gpu_done_ = false;
  else {
    tracing_foreach ("foreach", S__FILE__, S_LINENO);
    {...}
    end_tracing_foreach ("foreach", S__FILE__, S_LINENO);
  }
}

typedef struct {
  coord p, * box, n; // region
  int level; // level
} RegionParameters;

bool gpu_end_stencil (ForeachData * loop, const RegionParameters * region,
		      External * externals, const char * kernel);

macro2 foreach_stencil_generic (char flags, Reduce reductions,
				int _parallel, External * _externals, const char * _kernel)
{
  tracing_foreach ("foreach", S__FILE__, S_LINENO);
  static ForeachData _loop = { .fname = S__FILE__, .line = S_LINENO, .first = 1 };
  _loop.parallel = _parallel;
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
    _attribute[s.i].input = _attribute[s.i].output = false;
    _attribute[s.i].width = 0;
  }
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0}; NOT_UNUSED (point);
  RegionParameters _region = {0};
  
  {...}
  
#if PRINTIO
  if (baseblock) {
    fprintf (stderr, "%s:%d:", _loop.fname, _loop.line);
    for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i)
      if (_attribute[s.i].input || _attribute[s.i].output)
	fprintf (stderr, " %s:%d:%c:%d", _attribute[s.i].name, s.i,
		 _attribute[s.i].input && _attribute[s.i].output ? 'a' :
		 _attribute[s.i].input ? 'r' : 'w',
		 _attribute[s.i].width);
    fprintf (stderr, "\n");
  }
#endif // PRINTIO
  bool _first = _loop.first;
  _loop.first = 0; // to avoid warnings in check_stencil
  check_stencil (&_loop);
  _loop.first = _first;
  _gpu_done_ = gpu_end_stencil (&_loop, &_region, _externals, _kernel);
  _loop.first = 0;
  end_tracing_foreach ("foreach", S__FILE__, S_LINENO);
}

macro2 foreach_stencil (char flags, Reduce reductions,
			int _parallel, External * _externals, const char * _kernel)
{
  foreach_stencil_generic (flags, reductions, _parallel, _externals, _kernel)
    {...}
}

macro2 foreach_level_stencil (int _level, char flags, Reduce reductions,
			      int _parallel, External * _externals, const char * _kernel)
{
  foreach_stencil_generic (flags, reductions, _parallel, _externals, _kernel) {
    _region.level = _level + 1;
    {...}
  }
}

macro2 foreach_point_stencil (double _xp, double _yp, double _zp,
			      char flags, Reduce reductions,
			      int _parallel, External * _externals, const char * _kernel)
{
  foreach_stencil_generic (flags, reductions, _parallel, _externals, _kernel) {
    _region.p = (coord){ _xp, _yp, _zp };
    _region.n = (coord){ 1, 1 };
    {...}
  }
}

macro2 foreach_region_stencil (coord _p, coord _box[2], coord _n,
			       char flags, Reduce reductions,
			       int _parallel, External * _externals, const char * _kernel)
{
  foreach_stencil_generic (flags, reductions, _parallel, _externals, _kernel) {
    _region.p = _p, _region.box = _box, _region.n = _n;
    {...}
  }
}

macro2 foreach_vertex_stencil (char flags, Reduce reductions,
				  int _parallel, External * _externals,
				  const char * _kernel)
{
  foreach_stencil_generic (flags, reductions, _parallel, _externals, _kernel) {
    _loop.vertex = true;
    {...}
  }
}

macro2 foreach_face_stencil (char flags, Reduce reductions, const char * _order,
			     int _parallel, External * _externals,
			     const char * _kernel)
{
  foreach_stencil_generic (flags, reductions, _parallel, _externals, _kernel)
    {...}
}

macro2 foreach_coarse_level_stencil (int _level, char flags, Reduce reductions,
				     int _parallel, External * _externals,
				     const char * _kernel)
{
  foreach_level_stencil (_level, flags, reductions, _parallel, _externals, _kernel)
    {...}
}

macro2 foreach_level_or_leaf_stencil (int _level, char flags, Reduce reductions,
				      int _parallel, External * _externals,
				      const char * _kernel)
{
  foreach_level_stencil (_level, flags, reductions, _parallel, _externals, _kernel)
    {...}
}

@ifndef tracing
  @ def tracing(func, file, line) do {
    if (glFinish) glFinish();
    tracing(func, file, line);
  } while(0) @
  @ def end_tracing(func, file, line) do {
    if (glFinish) glFinish();
    end_tracing(func, file, line);
  } while(0) @
@endif

void realloc_ssbo()
{
  if (!datasize)
    return;
  
  GLint max_ssbo_size;
  GL_C (glGetIntegerv (GL_MAX_SHADER_STORAGE_BLOCK_SIZE, &max_ssbo_size));
#if 1
  GPUContext.max_ssbo_size = 128*(max_ssbo_size/128);
#else // for testing multi SSBOs
  GPUContext.max_ssbo_size = 128*(5*field_size()/2*sizeof(float)/128);
#endif
  
  size_t totalsize = field_size()*datasize;
  int nssbo = totalsize/GPUContext.max_ssbo_size;
  if (nssbo*GPUContext.max_ssbo_size < totalsize)
    nssbo++;

  GLint max_ssbo_blocks;
  GL_C (glGetIntegerv (GL_MAX_COMPUTE_SHADER_STORAGE_BLOCKS, &max_ssbo_blocks));
  if (nssbo > max_ssbo_blocks) {
    fprintf (stderr,
	     "%s:%d: error: cannot allocate %ld bytes\n"
	     "%s:%d: error: maximum allowed is %d x %ld bytes\n",
	     __FILE__, LINENO, totalsize, __FILE__, LINENO,
	     max_ssbo_blocks, GPUContext.max_ssbo_size);
    exit (1);
  }

#if DEBUGALLOC
  fprintf (stderr, "resizing from %ld to %ld nssbo %d pnssbo %d max_ssbo %ld\n", 
	   GPUContext.current_size, totalsize, nssbo, GPUContext.nssbo, GPUContext.max_ssbo_size);
#endif
  assert (totalsize > GPUContext.current_size);
  size_t size = max(GPUContext.nssbo - 1, 0)*GPUContext.max_ssbo_size;
  size_t current_size = GPUContext.current_size - size;
  assert (current_size >= 0 && current_size <= GPUContext.max_ssbo_size);
  GPUContext.current_size = totalsize;
  totalsize -= size;
  assert (totalsize >= 0);

  if (current_size > 0) {
    size_t size = min (totalsize, GPUContext.max_ssbo_size);
    totalsize -= size;
    if (current_size < GPUContext.max_ssbo_size) {
      GLuint ssbo;
      GL_C (glGenBuffers (1, &ssbo));
      GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, ssbo));
      GL_C (glBufferData (GL_SHADER_STORAGE_BUFFER, size, NULL, GL_DYNAMIC_READ));
      GL_C (glBindBuffer (GL_COPY_READ_BUFFER, GPUContext.ssbo[GPUContext.nssbo - 1]));
      GL_C (glBindBuffer (GL_COPY_WRITE_BUFFER, ssbo));
      GL_C (glCopyBufferSubData (GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, current_size));
      GL_C (glDeleteBuffers (1, &GPUContext.ssbo[GPUContext.nssbo - 1]));
      GPUContext.ssbo[GPUContext.nssbo - 1] = ssbo;
    }
#if DEBUGALLOC    
    else
      fprintf (stderr, "  skipping fully allocated %ld\n", current_size);
#endif
  }

#if DEBUGALLOC
  fprintf (stderr, "  need to allocate %ld\n", totalsize);
#endif
  if (nssbo > GPUContext.nssbo) {
    assert (totalsize > 0);
    GPUContext.ssbo = realloc (GPUContext.ssbo, nssbo*sizeof(GLuint));
#if DEBUGALLOC
    fprintf (stderr, "  allocating %d buffers\n", nssbo - GPUContext.nssbo);
#endif
    GL_C (glGenBuffers (nssbo - GPUContext.nssbo, GPUContext.ssbo + GPUContext.nssbo));
    while (GPUContext.nssbo < nssbo) {
      GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, GPUContext.ssbo[GPUContext.nssbo]));
      size_t size = min (totalsize, GPUContext.max_ssbo_size);
#if DEBUGALLOC
      fprintf (stderr, "  allocating buffer %d size %ld\n", GPUContext.nssbo, size);
#endif
      GL_C (glBufferData (GL_SHADER_STORAGE_BUFFER, size, NULL, GL_DYNAMIC_READ));
      totalsize -= size;
      GPUContext.nssbo++;
    }
  }
  GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, 0));
#if DEBUGALLOC
  fprintf (stderr, "done resizing %d %ld %ld\n", GPUContext.nssbo, GPUContext.current_size, totalsize);
#endif  
  assert (totalsize == 0);
}

static void gpu_cpu_sync_scalar (scalar s, char * sep, GLenum mode);

void realloc_scalar_gpu (int size)
{
  realloc_scalar (size);
  realloc_ssbo();
}

void gpu_boundary_level (scalar * list, int l)
{
  scalar * list1 = NULL;
  for (scalar s in list)
    if (s.gpu.stored > 0)
      list1 = list_prepend (list1, s);
  if (list1) {
    void cartesian_boundary_level (scalar * list, int l); 
    cartesian_boundary_level (list1, l);
    free (list1);
  }
}

#define realloc_scalar(size) realloc_scalar_gpu (size)
