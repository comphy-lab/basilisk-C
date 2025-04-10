double cpu_reduction (GLuint * src, size_t offset, size_t nb, const char op)
{
  GL_C (glMemoryBarrier (GL_BUFFER_UPDATE_BARRIER_BIT));
  double result;
  switch (op) {
  case '+': result =   0.;   break;
  case 'M': result = - HUGE; break;
  case 'm': result =   HUGE; break;
  default: assert (false);
  }
  offset *= sizeof(real), nb *= sizeof(real);
  int index = offset/GPUContext.max_ssbo_size;
  offset -= index*GPUContext.max_ssbo_size;
  while (nb) {
    GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, src[index]));
    size_t size = min (nb, GPUContext.max_ssbo_size - offset);
    real * a = glMapBufferRange (GL_SHADER_STORAGE_BUFFER, offset, size, GL_MAP_READ_BIT);
    assert (a);
    switch (op) {
    case '+':
      for (size_t i = 0; i < size/sizeof(real); i++, a++)
	result += *a;
      break;
    case 'M':
      for (size_t i = 0; i < size/sizeof(real); i++, a++)
	result = fmax (result, *a);
      break;
    case 'm':
      for (size_t i = 0; i < size/sizeof(real); i++, a++)
	result = fmin (result, *a);
      break;
    default: assert (false);
    }
    assert (glUnmapBuffer (GL_SHADER_STORAGE_BUFFER));
    nb -= size, offset = 0, index++;
  }
  GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, 0));
  return result;
}

double gpu_reduction (size_t offset, const char op, const RegionParameters * region, size_t nb)
{
  const int stride = 64, nwgr = 64;
  bool is_foreach_point = (region->n.x == 1 && region->n.y == 1);
  if (!is_foreach_point && nb < nwgr*stride)
    return cpu_reduction (GPUContext.ssbo, offset, nb, op);
  
  GLuint * br = gpu_grid->reduct;
  if (!br[0]) {
    GL_C (glGenBuffers (2, br));
    for (int i = 0; i < 2; i++) {
      GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, br[i]));
      size_t size = (sq((size_t)N + 1)/stride + 1)*sizeof(real);
      assert (size < GPUContext.max_ssbo_size); // must fit within a single SSBO
      GL_C (glBufferData (GL_SHADER_STORAGE_BUFFER, size, NULL, GL_DYNAMIC_READ));
    }
    GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, 0));
  }
  
  const char * operation;
  static const char * opsum =
    "  real reduct = 0.;\n"
    "  for (uint j = 0; j < stride; j++, index++) {\n"
    "    val = _data_val(index);\n"
    "    reduct += val;\n"
    "  }\n";
  static const char * opmin =
    "  real reduct = val;\n"
    "  for (uint j = 0; j < stride; j++, index++) {\n"
    "    val = _data_val(index);\n"
    "    if (val < reduct) reduct = val;\n"
    "  }\n";
  static const char * opmax =
    "  real reduct = val;\n"
    "  for (uint j = 0; j < stride; j++, index++) {\n"
    "    val = _data_val(index);\n"
    "    if (val > reduct) reduct = val;\n"
    "  }\n";
  switch (op) {
  case '+': operation = opsum; break;
  case 'M': operation = opmax; break;
  case 'm': operation = opmin; break;
  default: // unknown reduction operation
    assert (false);
  }

  char nwgrs[20], strides[20];
  snprintf (nwgrs, 19, "%d", nwgr);
  snprintf (strides, 19, "%d", stride);
  
  char * fs =
    str_append (NULL,
		"#version 430\n", glsl_preproc,
		"layout (std430, binding = 0) writeonly buffer _reduct_layout {"
		" real _reduct[]; };\n"
		"layout (std430, binding = 1) readonly buffer _data_layout { real f[]; } _data");
  if (GPUContext.nssbo > 1) {
    char a[20], s[30];
    snprintf (a, 19, "%d", GPUContext.nssbo);
    snprintf (s, 29, "%ld", GPUContext.max_ssbo_size/sizeof(real));
    fs = str_append (fs, "[", a, "];\n"
		     "#define _data_val(index) _data[(index)/", s, "].f[(index)%", s, "]\n");
  }
  else
    fs = str_append (fs, ";\n"
		     "#define _data_val(index) _data.f[index]\n");
  fs = str_append (fs,
		   "layout (location = 3) uniform uint offset;\n"
		   "layout (location = 4) uniform uint nb;\n"
		   "layout (location = 5) uniform uint nbr;\n"
		   "layout (local_size_x = ", nwgrs, ") in;\n"
		   "void main() {\n"
		   "if (gl_GlobalInvocationID.x < nb) {\n"
		   "  uint stride = ", strides, ";\n"
		   "  uint index = stride*gl_GlobalInvocationID.x;\n"
		   "  if (index + stride > nbr) stride = nbr - index;\n"
		   "  index += offset;\n"
		   "  real val = _data_val(index);\n",
		   operation,
		   "  _reduct[gl_GlobalInvocationID.x] = reduct;\n"
		   "}}\n");
  Shader * shader;
  Adler32Hash hash;
  a32_hash_init (&hash);
  a32_hash_add (&hash, fs, strlen (fs));
  shader = load_shader (fs, a32_hash (&hash), NULL);
  assert (shader);
  if (shader->id != GPUContext.current_shader) {
    GL_C (glUseProgram (shader->id));
    GPUContext.current_shader = shader->id;
  }
  const GLint loffset = 3, lnb = 4, lnbr = 5;

  int start = offset*sizeof(real)/GPUContext.max_ssbo_size;
  offset -= start*GPUContext.max_ssbo_size/sizeof(real);
  for (int i = start; i < GPUContext.nssbo; i++)
    GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 1 + i - start, GPUContext.ssbo[i]));

  if (is_foreach_point) {
    real result = 0.;
    int i = (region->p.x - X0)/L0*N;
    int j = (region->p.y - Y0)/L0*N;
    if (i >= 0 && i < N && j >= 0 && j < N) {
      offset += i*N + j;
      GL_C (glUniform1ui (loffset, offset));
      GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 0, br[0]));
      GL_C (glUniform1ui (lnbr, 1));
      GL_C (glUniform1ui (lnb, 1));
      GL_C (glMemoryBarrier (GL_SHADER_STORAGE_BARRIER_BIT));
      GL_C (glDispatchCompute (1, 1, 1));
      GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, br[0]));
      GL_C (glGetBufferSubData (GL_SHADER_STORAGE_BUFFER, 0, sizeof (real), &result));
      GL_C (glBindBuffer (GL_SHADER_STORAGE_BUFFER, 0));
    }
    return result;
  }
  
  GL_C (glUniform1ui (loffset, offset));
  int src = 0, dst = 1;
  while (nb >= nwgr*stride) {
    GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 0, br[dst]));
    GL_C (glUniform1ui (lnbr, nb));
    if (nb % stride) {
      nb /= stride;
      nb++;
    }
    else
      nb /= stride;
    int ng = nb/nwgr;
    if (ng*nwgr < nb)
      ng++;
    GL_C (glUniform1ui (lnb, nb));
    GL_C (glMemoryBarrier (GL_SHADER_STORAGE_BARRIER_BIT));
    GL_C (glDispatchCompute (ng, 1, 1));
    swap (int, src, dst);
    if (offset) {
      GL_C (glUniform1ui (loffset, 0));
      offset = 0;
    }
    GL_C (glBindBufferBase (GL_SHADER_STORAGE_BUFFER, 1, br[src]));
  }

  return cpu_reduction (br + src, 0, nb, op);
}
