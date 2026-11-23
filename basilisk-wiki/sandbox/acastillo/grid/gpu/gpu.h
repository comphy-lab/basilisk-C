char * gpu_errors (const char * errors, const char * source, char * fout,
                   const char * lang);

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
  static ForeachData _loop = { .fname = S__FILE__, .func = S__func__, .line = S_LINENO, .first = 1 };
  _loop.parallel = _parallel;
  if (baseblock)
    for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
      _attribute[s.i].stencil.io = 0;
      _attribute[s.i].stencil.width = 0;
    }
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0}; NOT_UNUSED (point);
  RegionParameters _region = {0};
  
  {...}
  
#if PRINTIO
  if (baseblock) {
    fprintf (stderr, "%s:%d:", _loop.fname, _loop.line);
    for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i)
      if ((_attribute[s.i].stencil.io & s_input) || (_attribute[s.i].stencil.io & s_output))
	fprintf (stderr, " %s:%d:%c:%d", _attribute[s.i].name, s.i,
                 (_attribute[s.i].stencil.io & s_input) &&
                 (_attribute[s.i].stencil.io & s_output) ? 'a' :
		 (_attribute[s.i].stencil.io & s_input) ? 'r' : 'w',
		 _attribute[s.i].stencil.width);
    fprintf (stderr, "\n");
  }
#endif // PRINTIO
  bool _first = _loop.first;
  _loop.first = 0; // to avoid warnings in check_stencil
  check_stencil (&_loop);
#if _MPI
  {
    bool _any_f = false;
    foreach_dimension()
      if (_loop.listf.x) _any_f = true;
    if (_any_f || _loop.listc) {
      scalar * _merged = NULL;
      for (scalar _s in baseblock) {
	bool _in = list_lookup (_loop.listc, _s);
	foreach_dimension()
	  if (_loop.listf.x && list_lookup (_loop.listf.x, _s))
	    _in = true;
	if (_in)
	  _merged = list_append (_merged, _s);
      }
      free (_loop.listc);
      _loop.listc = _merged;
      foreach_dimension() {
	free (_loop.listf.x);
	_loop.listf.x = NULL;
      }
    }
  }
#endif
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
    _region.n = (coord){ 1, 1, 1 };
    {...}
  }
}

macro2 foreach_region_stencil (coord _p, coord _box[2], coord _n,
			       char flags, Reduce reductions,
			       int _parallel, External * _externals, const char * _kernel)
{
  foreach_stencil_generic (flags, reductions, _parallel, _externals, _kernel) {
#if _CUDA
    _loop.parallel = 2;
#endif
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

void gpu_boundary_level (scalar * list, int l)
{
  scalar * list_cpu = NULL;
#if _MPI && _CUDA
  scalar * list_gpu = NULL;
#endif
  for (scalar s0 in list)
    if (!is_constant(s0) && s0.i >= 0) {
      scalar s = s0.block < 0 ? (scalar){s0.i + s0.block} : s0;
      if (s.block <= 0) continue;
      if (s.gpu.stored > 0)
	list_cpu = list_prepend (list_cpu, s);
#if _MPI && _CUDA
      else
	list_gpu = list_append (list_gpu, s);
#endif
    }
  if (list_cpu) {
    void cartesian_boundary_level (scalar * list, int l);
    cartesian_boundary_level (list_cpu, l);
    free (list_cpu);
  }
#if _MPI && _CUDA
  if (list_gpu) {
    void mpi_boundary_level_gpu (scalar * list, int l);
    mpi_boundary_level_gpu (list_gpu, l);
    free (list_gpu);
  }
#endif
}

void realloc_scalar_gpu (int size);
void realloc_scalar_cpu (int size) { realloc_scalar(size); }

#define realloc_scalar(size) realloc_scalar_gpu (size)
