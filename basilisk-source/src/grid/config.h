#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#if NVTX
 @include <nvtx3/nvToolsExt.h>
#endif

@define unmap(x,y)
@define trash(x)  // data trashing is disabled by default. Turn it on with
                  // -DTRASH=1

auto macro2 BEGIN_FOREACH() {{...}}

#if _OPENMP
@ include <omp.h>
@ define OMP(x) Pragma(#x)

macro OMP_SERIAL()
{
  @undef OMP
  @define OMP(x)
  ; // necessary so that the preproc above is included
  {...}
  @undef OMP
  @define OMP(x) Pragma(#x)
  ; // necessary so that the preproc above is included
}
   
#elif _MPI

@ define OMP(x)
macro OMP_SERIAL() {{...}}

@ include <mpi.h>
static int mpi_rank, mpi_npe;
@ define tid() mpi_rank
@ define pid() mpi_rank
@ define npe() mpi_npe

#else // not MPI, not OpenMP

@ define OMP(x)
macro OMP_SERIAL() {{...}}

#endif // _MPI

#if _CADNA
# include <cadna.h>
#endif // CADNA

#if __cplusplus
# define delete delete_qcc
# define right right_qcc
# define left left_qcc
# define norm norm_qcc
# define new new_qcc
#endif // _cplusplus

@define _NVARMAX 65536
@define is_constant(v) ((v).i >= _NVARMAX)
@define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

@define systderr  stderr
@define systdout  stdout
#if _MPI
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
@ def not_mpi_compatible()
do {
  if (npe() > 1) {
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);
    exit (1);
  }
} while(0)
@
@ define system(command) (pid() == 0 ? system(command) : 0)
#else
@ define qstderr() stderr
@ define qstdout() stdout
@ define ferr      stderr
@ define fout      stdout
@ define not_mpi_compatible()
#endif

// redirect assert
static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (stderr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}
#undef assert
#ifndef LINENO
# define LINENO __LINE__
#endif
#define assert(a) if (!(a)) qassert (__FILE__, LINENO, #a)

// Memory tracing
@define sysmalloc malloc
@define syscalloc calloc
@define sysrealloc realloc
@define sysfree free
@define systrdup strdup

#if MTRACE
# include "mtrace.h"
#else // !MTRACE
@ define pmalloc(s,func,file,line)    malloc(s)
@ define pcalloc(n,s,func,file,line)  calloc(n,s)
@ define prealloc(p,s,func,file,line) realloc(p,s)
@ define pfree(p,func,file,line)      free(p)
@ define pstrdup(s,func,file,line)    strdup(s)
#endif // !MTRACE

#define qrealloc(p, size, type) p = (type *) realloc (p, (size)*sizeof(type))
#define qmalloc(size, type) ((type *) malloc ((size)*sizeof(type)))
#define qcalloc(size, type) ((type *) calloc (size, sizeof(type)))

#include "array.h"

// Function tracing

#if TRACE == 1 // with Extrae library
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func     = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = strdup (func);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  assert (t->stack.len > 0);
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    free (*func);
  free (t->index.p);
  free (t->stack.p);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}
#if 0
#define TRACE_TYPE(func) (strncmp (func, "mpi_", 4) ?		\
			  &trace_func : &trace_func)
#else
#define TRACE_TYPE(func) &trace_func
#endif
@  define tracing(func, file, line)     trace_push (TRACE_TYPE(func), func)
@  define end_tracing(func, file, line) trace_pop (TRACE_TYPE(func), func)

#elif TRACE // built-in function tracing

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if _MPI
  double min, max;
#endif // _MPI
} TraceIndex;
				      
struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
		       double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {strdup(func), strdup(file), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));
#if 0
  fprintf (stderr, "trace %s:%s:%d t: %g sum: %g\n",
	   func, file, line, t[0], t[1]);
#endif
#if NVTX
  nvtxRangePush (func);
#endif
}

static void end_tracing (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  assert (Trace.stack.len >= 2*sizeof(double));
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];
#if 0
  fprintf (stderr, "end trace %s:%s:%d ts: %g te: %g dt: %g sum: %g\n",
	   func, file, line, t[0], te, dt, t[1]);
#endif
  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
#if NVTX
  nvtxRangePop();
#endif
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if _MPI
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if _MPI
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self,  min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self,  max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
	      self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
	      tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif // _MPI
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
	       t->calls, t->total, t->self, t->self*100./total);
#if _MPI
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    free (t->func), free (t->file);

  free (Trace.index.p);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;
  
  free (Trace.stack.p);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else // disable tracing
@  define tracing(...)
@  define end_tracing(...)
#endif

// OpenMP / MPI
  
#if _OPENMP

@define tid() omp_get_thread_num()
@define pid() 0
@define npe() omp_get_num_threads()
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)

#elif _MPI

static bool in_prof = false;
static double prof_start, _prof;
@def prof_start(name)
  assert (!in_prof); in_prof = true;
  prof_start = MPI_Wtime();
@
@def prof_stop()
  assert (in_prof); in_prof = false;
  _prof = MPI_Wtime();
  mpi_time += _prof - prof_start;
@

#if FAKE_MPI
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)
#else // !FAKE_MPI
trace
int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
		     MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  return MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);
}

@def mpi_all_reduce(v,type,op) {
  prof_start ("mpi_all_reduce");
  union { int a; float b; double c;} global;
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);
  memcpy (&(v), &global, sizeof (v));
  prof_stop();
}
@

trace
void mpi_all_reduce_array (void * v, MPI_Datatype datatype, MPI_Op op, int elem)
{
  size_t size;
  if (datatype == MPI_DOUBLE) size = sizeof (double);
  else if (datatype == MPI_INT) size = sizeof (int);
  else if (datatype == MPI_LONG) size = sizeof (long);
  else if (datatype == MPI_C_BOOL) size = sizeof (bool);
  else if (datatype == MPI_UNSIGNED_CHAR) size = sizeof (unsigned char);
  else {
    fprintf (stderr, "unknown reduction type\n");
    fflush (stderr);
    abort();
  }
  void * global = malloc (elem*size), * tmp = malloc (elem*size);
  memcpy (tmp, v, elem*size);
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);
  memcpy (v, global, elem*size);
  free (global), free (tmp);
}
#endif // !FAKE_MPI

@define QFILE FILE // a dirty trick to avoid qcc 'static FILE *' rule

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else // not MPI, not OpenMP

@define tid() 0
@define pid() 0
@define npe() 1
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)

#endif // not MPI, not OpenMP

macro2 OMP_PARALLEL() {{...}}
@define OMP_PARALLEL(...) OMP(omp parallel S__VA_ARGS__)

@define NOT_UNUSED(x) (void)(x)

macro2 VARIABLES() { _CATCH; }
@define _index(a,m)    (a.i)
@define val(a,k,l,m)   data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;

/* undefined value */
/* Initialises unused memory with "signaling NaNs".  
 * This is probably not very portable, tested with
 * gcc (Debian 4.4.5-8) 4.4.5 on Linux 2.6.32-5-amd64.
 * This blog was useful:
 *   http://codingcastles.blogspot.co.nz/2008/12/nans-in-c.html 
 */
@if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
@ if __APPLE__
@   include <stdint.h>
@   include "fp_osx.h"
@ endif
@if _GPU
@  define enable_fpe(flags)
@else
@  define enable_fpe(flags)  feenableexcept (flags)
@endif
@  define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  assert (sizeof (int64_t) == sizeof (double));
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
@else // !((_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA && !_GPU)
@  define undefined ((double) DBL_MAX)
@  define enable_fpe(flags)
@  define disable_fpe(flags)
static void set_fpe (void) {}
@endif

// Pipes

static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qrealloc (qpopen_pipes, n + 2, FILE *);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  free (qpopen_pipes);
  qpopen_pipes = NULL;
}

#define popen  qpopen
#define pclose qpclose

// files with pid

FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}

#include "../ast/symbols.h"

enum typedef_kind_t {
  sym_SCALAR = sym_root + 1,
  sym_VECTOR,
  sym_TENSOR,
  sym_COORD,
  sym__COORD,
  sym_VEC4,
  sym_IVEC
};

@define attroffset(x) (offsetof(_Attributes,x))

/**
These are placeholders for internally-defined macros. */

typedef int Reduce;

macro2 foreach_face (char flags = 0, Reduce reductions = None,
			const char * order = "xyz")
{{...}}
macro2 einstein_sum() {{...}}
macro2 diagonalize (int a) {{...}}

/**
Macros overloaded by the interpreter. */

@define dimensional(...)
#define show_dimension(...) show_dimension_internal (__VA_ARGS__ + 10293847566574839201.)
@define show_dimension_internal(...)
@define display_value(...)
@define interpreter_verbosity(...)
