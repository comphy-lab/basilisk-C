#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

@if _OPENMP
@ include <omp.h>
@ define OMP(x) Pragma(#x)
@elif _MPI

@ define OMP(x)

@ include <mpi.h>
static int mpi_rank, mpi_npe;
@ define tid() mpi_rank
@ define pid() mpi_rank
@ define npe() mpi_npe

@else // not MPI, not OpenMP

@ define OMP(x)

@endif // _MPI

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

#define pi 3.14159265358979
#undef HUGE
#define HUGE 1e30
#define nodata HUGE
@define _NVARMAX 65536
@define is_constant(v) ((v).i >= _NVARMAX)
@define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

@define max(a,b) ((a) > (b) ? (a) : (b))
@define min(a,b) ((a) < (b) ? (a) : (b))
@define sq(x) ((x)*(x))
@define cube(x) ((x)*(x)*(x))
@define sign(x) ((x) > 0 ? 1 : -1)
@define sign2(x) ((x) > 0 ? 1 : (x) < 0 ? -1 : 0)
@define noise() (1. - 2.*rand()/(double)RAND_MAX)
@define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
#define swap(type,a,b) do { type __tmp = a; a = b; b = __tmp; } while(0)
@define unmap(x,y)

@define trash(x)  // data trashing is disabled by default. Turn it on with
                  // -DTRASH=1

@define systderr  stderr
@define systdout  stdout

@if _MPI
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
@else
@ define qstderr() stderr
@ define qstdout() stdout
@ define ferr      stderr
@ define fout      stdout
@ define not_mpi_compatible()
@endif

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

@if MTRACE

struct {
  FILE * fp;                     // trace file
  size_t total, max;             // current and maximum allocated memory
  size_t overhead, maxoverhead;  // current and maximum profiling overhead
  size_t nr;                     // current number of records
  size_t startrss, maxrss;       // starting and maximum system ressource usage
  char * fname;                  // trace file name
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
	     c, f->id, pmtrace.nr, pmtrace.total, f->total);
@if _GNU_SOURCE
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
@endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
			    const char * func, const char * file, int line,
			    char c)
{
  assert (d != NULL);
  OMP (omp critical)
  {
    d->id = pmfunc_index(func, file, line);
    d->size = size;
    pmfunc * f = &pmfuncs[d->id - 1];
    f->total += size;
    if (f->total > f->max)
      f->max = f->total;
    pmtrace.total += size;
    pmtrace.overhead += sizeof(pmdata);
    if (pmtrace.total > pmtrace.max) {
      pmtrace.max = pmtrace.total;
      pmtrace.maxoverhead = pmtrace.overhead;
    }
    pmfunc_trace (f, c);
  }
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", stderr);
    if (d->size == 0)
      fputs (", possible double free()", stderr);
    else
      fputs (", not traced?", stderr);
    fputs (", aborting...\n", stderr);
    abort();
    return ptr;
  }
  else
  OMP (omp critical)
  {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (stderr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
	       f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (stderr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
	       pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
  }
  return d;
}

static void * pmalloc (size_t size,
		       const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
		       size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
		       const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
			const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
					      sizeof(pmdata) + size),
		       size, func, file, line, '>');
}

static void pfree (void * ptr,
		   const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
		       const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

@if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
@endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
@if MTRACE < 3
  fprintf (stderr,
	   "*** MTRACE: max resident  set size: %10ld bytes\n"
	   "*** MTRACE: max traced memory size: %10ld bytes"
	   " (tracing overhead %.1g%%)\n"
	   "%10s    %20s   %s\n",
	   pmtrace.maxrss*1024,
	   pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
	   "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (stderr, "%10ld    %20s   %s:%d\n",
	     p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
	     "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
	     "total(\"%s\") w l t 'total'",
	     fname,
	     pmtrace.startrss*1024.,
	     pmtrace.startrss*1024.,
	     fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
	       ",func(\"%s\",%d) w l t '%s'",
	       fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (stderr,
	     "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
	     fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
@endif // MTRACE < 3

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (stderr, "%s:%d: error: %ld bytes leaked here\n",
	       p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
@if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", stderr);
@endif
    pmfuncs_free();
  }
}

@else // !MTRACE
@ define pmalloc(s,func,file,line)    malloc(s)
@ define pcalloc(n,s,func,file,line)  calloc(n,s)
@ define prealloc(p,s,func,file,line) realloc(p,s)
@ define pfree(p,func,file,line)      free(p)
@ define pstrdup(s,func,file,line)    strdup(s)
@endif // !MTRACE

#define qrealloc(p, size, type) p = (type *) realloc (p, (size)*sizeof(type))
#define qmalloc(size, type) ((type *) malloc ((size)*sizeof(type)))
#define qcalloc(size, type) ((type *) calloc (size, sizeof(type)))

// Arrays

typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = qmalloc (1, Array);
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  free (a->p);
  free (a);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = realloc (a->p, a->max);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = realloc (a->p, a->len);
  free (a);
  return p;
}

// Function tracing

@if TRACE == 1 // with Extrae library
@include <extrae_user_events.h>

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

@elif TRACE // built-in function tracing

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
@if _MPI
  double min, max;
@endif // _MPI
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
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

@if _MPI
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
@endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
@if _MPI
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
@endif // _MPI
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
	       t->calls, t->total, t->self, t->self*100./total);
@if _MPI
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
@endif
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

@else // disable tracing
@  define tracing(...)
@  define end_tracing(...)
@endif

// OpenMP / MPI
  
@if _OPENMP

@define tid() omp_get_thread_num()
@define pid() 0
@define npe() omp_get_num_threads()
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)

@elif _MPI

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

@if FAKE_MPI
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)
@else // !FAKE_MPI
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
@def mpi_all_reduce_array(v,type,op,elem) {
  prof_start ("mpi_all_reduce");
  type * global = malloc ((elem)*sizeof(type)), * tmp = malloc ((elem)*sizeof(type));
  for (int i = 0; i < elem; i++)
    tmp[i] = (v)[i];
  MPI_Datatype datatype;
  if (!strcmp(#type, "double")) datatype = MPI_DOUBLE;
  else if (!strcmp(#type, "int")) datatype = MPI_INT;
  else if (!strcmp(#type, "long")) datatype = MPI_LONG;
  else if (!strcmp(#type, "bool")) datatype = MPI_C_BOOL;
  else if (!strcmp(#type, "unsigned char")) datatype = MPI_UNSIGNED_CHAR;
  else {
    fprintf (stderr, "unknown reduction type '%s'\n", #type);
    fflush (stderr);
    abort();
  }
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);
  for (int i = 0; i < elem; i++)
    (v)[i] = global[i];
  free (global), free (tmp);
  prof_stop();
}
@

@endif // !FAKE_MPI

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
@if MTRACE == 1
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
@endif
  }
}

@else // not MPI, not OpenMP

@define tid() 0
@define pid() 0
@define npe() 1
@define mpi_all_reduce(v,type,op)
@define mpi_all_reduce_array(v,type,op,elem)

@endif // not MPI, not OpenMP

@define OMP_PARALLEL() OMP(omp parallel)

@define NOT_UNUSED(x) (void)(x)

@define VARIABLES      _CATCH;
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
@  define enable_fpe(flags)  feenableexcept (flags)
@  define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  assert (sizeof (int64_t) == sizeof (double));
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
@else // !((_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA)
@  define undefined ((double) DBL_MAX)
@  define enable_fpe(flags)
@  define disable_fpe(flags)
static void set_fpe (void) {}
@endif

// the grid
typedef struct {
  long n;       // number of (leaf) cells for this process
  long tn;      // number of (leaf) cells for all processes
  int depth;    // the depth for this process
  int maxdepth; // the maximum depth for all processes
} Grid;
Grid * grid = NULL;
// coordinates of the lower-left corner of the box
double X0 = 0., Y0 = 0., Z0 = 0.;
// size of the box
double L0 = 1. [1];
// number of grid points
#if dimension <= 2
int N = 64;
#else
int N = 16;
#endif

typedef struct { int i; } scalar;

typedef struct {
  scalar x;
#if dimension > 1
  scalar y;
#endif
#if dimension > 2
  scalar z;
#endif
} vector;

typedef struct {
  scalar * x;
#if dimension > 1
  scalar * y;
#endif
#if dimension > 2
  scalar * z;
#endif
} vectorl;

typedef struct {
  vector x;
#if dimension > 1
  vector y;
#endif
#if dimension > 2
  vector z;
#endif
} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
			   omp_out.x += omp_in.x,
			   omp_out.y += omp_in.y,
			   omp_out.z += omp_in.z))

#if dimension == 1
# define norm(v) fabs(v.x[])
# define dv() (Delta*cm[])
#elif dimension == 2
# define norm(v) (sqrt(sq(v.x[]) + sq(v.y[])))
# define dv() (sq(Delta)*cm[])
#else // dimension == 3
# define norm(v) (sqrt(sq(v.x[]) + sq(v.y[]) + sq(v.z[])))
# define dv() (cube(Delta)*cm[])
#endif

void normalize (coord * n)
{
  double norm = 0.;
  foreach_dimension()
    norm += sq(n->x);
  norm = sqrt(norm);
  foreach_dimension()
    n->x /= norm;
}

void origin (double x = 0., double y = 0., double z = 0.) {
  X0 = x; Y0 = y; Z0 = z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }

// boundary conditions for each direction/variable

#if dimension == 1
  enum { right, left };
#elif dimension == 2
  enum { right, left, top, bottom };
#else
  enum { right, left, top, bottom, front, back };
#endif
int nboundary = 2*dimension;

#define none -1

@define _dirichlet(expr, ...)             (2.*(expr) - val(_s,0,0,0))
@define _dirichlet_homogeneous(...)       (- val(_s,0,0,0))
@define _dirichlet_face(expr,...)         (expr)
@define _dirichlet_face_homogeneous(...)  (0.)
@define _neumann(expr,...)                (Delta*(expr) + val(_s,0,0,0))
@define _neumann_homogeneous(...)         (val(_s,0,0,0))

double  * _constant = NULL;
size_t datasize = 0;
typedef struct _Point Point;

#include "grid/boundaries.h"

// attributes for each scalar

typedef struct {
  double (** boundary)             (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient)              (double, double, double);
  void   (* delete)                (scalar);
  char * name;
  struct {
    int x;
#if dimension > 1
    int y;
#endif
#if dimension > 2
    int z;
#endif
  } d; // staggering
  vector v;
  int face;
  bool   nodump, freed;
  int    block;
  scalar * depends; // boundary conditions depend on other fields
} _Attributes;

static _Attributes * _attribute = NULL;

#define foreach_block() // treatment of block data is disabled by default
#define foreach_blockf(s)

// lists

int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  for (scalar s in list) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  qrealloc (list, len + 2, scalar);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  qrealloc (list, len + 2, scalar);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  for (scalar t in list)
    if (t.i == s.i)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    for (scalar s1 in l)
      if (s1.i == s.i)
	return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    for (scalar s in l)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  for (scalar s in l2)
    l3 = list_append (l3, s);
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  for (scalar s in l)
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", s.name);
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  for (vector v in list) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  qrealloc (list, len + 2, vector);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  for (vector w in list) {
    bool id = true;
    foreach_dimension()
      if (w.x.i != v.x.i)
	id = false;
    if (id)
      return list;
  }
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    for (vector v in l)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
    foreach_dimension() {
      assert (s->i >= 0);
      v.x = *s++;
    }
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  for (tensor t in list) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  qrealloc (list, len + 2, tensor);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
    foreach_dimension() {
      assert (v->x.i >= 0);
      t.x = *v++;
    }
    list = tensors_append (list, t);
  }
  return list;
}

static inline bool is_vertex_scalar (scalar s)
{
  foreach_dimension()
    if (s.d.x != -1)
      return false;
  return true;
}

scalar * all = NULL; // all the fields
scalar * baseblock = NULL; // base block fields

// basic methods

scalar (* init_scalar)        (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector)        (vector, const char *);
vector (* init_face_vector)   (vector, const char *);
tensor (* init_tensor)        (tensor, const char *);
void   (* scalar_clone)       (scalar, scalar);

#define vector(x) (*((vector *)&(x)))

// events 

typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL; // all events

int iter = 0, inext = 0; // current step and step of next event
double t = 0, tnext = 0; // current time and time of next event
void init_events (void);
void event_register (Event event);
static void _init_solver (void);

void init_solver()
{
  Events = malloc (sizeof (Event));
  Events[0].last = 1;
  _attribute = calloc (datasize/sizeof(double), sizeof (_Attributes));
  int n = datasize/sizeof(double);
  all = (scalar *) malloc (sizeof (scalar)*(n + 1));
  baseblock = (scalar *) malloc (sizeof (scalar)*(n + 1));
  for (int i = 0; i < n; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[n].i = all[n].i = -1;
@if _CADNA
  cadna_init (-1);
@endif
@if _MPI
  mpi_init();
@elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
@endif
}

// timers

@if _MPI
static double mpi_time = 0.;
@endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
@if _MPI
  t.tm = mpi_time;
@endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) + 
	  (tvend.tv_usec - t.tv.tv_usec)/1e6);
}

// Constant fields // fixme: should be removed

const face vector zerof[] = {0.,0.,0.};
const face vector unityf[] = {1.,1.,1.};
const scalar unity[] = 1.;
const scalar zeroc[] = 0.;

// Metric

const face vector unityf0[] = {1.[0],1.[0],1.[0]};
const scalar unity0[] = 1.[0];
(const) face vector fm = unityf0; // fixme: should just be fm = {1,1,1};
(const) scalar cm = unity0; // fixme: should just be cm = 1.;

// Embedded boundaries
// these macros are overloaded in embed.h

#define SEPS 0.
#define face_gradient_x(a,i) ((a[i] - a[i-1])/Delta)
#define face_gradient_y(a,i) ((a[0,i] - a[0,i-1])/Delta)
#define face_gradient_z(a,i) ((a[0,0,i] - a[0,0,i-1])/Delta)
#define face_value(a,i)      ((a[i] + a[i-1])/2.)
#define center_gradient(a)   ((a[1] - a[-1])/(2.*Delta))

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

// matrices

void * matrix_new (int n, int p, size_t size)
{
  void ** m = qmalloc (n, void *);
  char * a = qmalloc (n*p*size, char);
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
	for (k = 0; k < n; k++) {
	  if (ipiv[k] == -1) {
	    if (fabs (m[j][k]) >= big) {
	      big = fabs (m[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }
	}
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++) 
	swap (double, m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
	dum = m[ll][icol];
	m[ll][icol] = 0.0;
	for (l = 0; l < n; l++)
	  m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
	swap (double, m[k][indxr[l]], m[k][indxc[l]]);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  free (((void **) m)[0]);
  free (m);
}

// Solver cleanup

typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}

// Default objects to display

static char * display_defaults = NULL;

static void free_display_defaults() {
  free (display_defaults);
}

void display (const char * commands, bool overwrite = false)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (overwrite) {
    free (display_defaults);
    display_defaults = malloc (strlen(commands) + 2);
    strcpy (display_defaults, "@");
    strcat (display_defaults, commands);
  }
  else {
    if (!display_defaults)
      display_defaults = strdup ("@");
    display_defaults =
      realloc (display_defaults,
	       strlen(display_defaults) + strlen(commands) + 1);
    strcat (display_defaults, commands);
  }
}

#define display_control(val, ...)

typedef struct {
  double x;
#if dimension > 1
  double y;
#endif
#if dimension > 2
  double z;
#endif
} _coord;

#if LAYERS
# include "grid/layers.h"
#endif

#include "grid/stencils.h"

#define dirichlet(expr)                 _dirichlet(expr, point, neighbor, _s, data)
#define dirichlet_face(expr)            _dirichlet_face(expr, point, neighbor, _s, data)
#define neumann(expr)                   _neumann(expr, point, neighbor, _s, data)

typedef struct {
  coord x, y, z;
} mat3;

OMP(omp declare reduction (+ : mat3 :
			   omp_out.x.x += omp_in.x.x,
			   omp_out.x.y += omp_in.x.y,
			   omp_out.x.z += omp_in.x.z,
			   omp_out.y.x += omp_in.y.x,
			   omp_out.y.y += omp_in.y.y,
			   omp_out.y.z += omp_in.y.z,
			   omp_out.z.x += omp_in.z.x,
			   omp_out.z.y += omp_in.z.y,
			   omp_out.z.z += omp_in.z.z
			   ))
