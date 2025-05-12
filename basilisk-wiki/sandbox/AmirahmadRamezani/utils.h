// Default parameters
// maximum timestep
double DT = 1e10;
// CFL number
double CFL = 0.5;

void runge_kutta (int stages,
		  double t, double dt,
		  int nv, scalar f[nv], scalar df[stages][nv], 
		  void (* advance) (double t, scalar f[nv], scalar df[nv]),
		  void (* update)  (double t, scalar f[nv]))
{
  switch (stages) {
  case 1:
    (* advance) (t, f, df[0]);
    foreach()
      for (int v = 0; v < nv; v++)
	f[v][] += df[0][v][]*dt;
    (* update) (t + dt, f);
    break;

  case 2:
    (* advance) (t, f, df[0]);
    foreach()
      for (int v = 0; v < nv; v++)
	df[0][v][] = f[v][] + df[0][v][]*dt/2.;
    (* update) (t + dt/2., df[0]);

    (* advance) (t + dt/2., df[0], df[1]);
    foreach()
      for (int v = 0; v < nv; v++)
	f[v][] += df[1][v][]*dt;
    (* update) (t + dt, f);
    break;

  default:
    /* not implemented yet */
    assert(false);
  }
}

double change (scalar v, scalar vn)
{
  double max = 0.;
  foreach(reduction(max:max)) {
    double dv = fabs (v[] - vn[]);
    if (dv > max)
      max = dv;
    vn[] = v[];
  }
  return max;
}

typedef struct {
  double cpu;   // CPU time (sec)
  double real;  // Wall-clock time (sec)
  double speed; // Speed (points.steps/sec)
  double min;   // Minimum MPI time (sec)
  double avg;   // Average MPI time (sec)
  double max;   // Maximum MPI time (sec)
  size_t tnc;   // Number of grid points
  long   mem;   // Maximum resident memory (kB)
} timing;

timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
@if _MPI
  s.avg = mpi_time - t.tm;
@endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    foreach(reduction(+:tnc)) tnc++;
    s.tnc = tnc;
    tnc *= i;
  }
  else
    s.tnc = tnc;
@if (_GNU_SOURCE || _DARWIN_C_SOURCE)
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
@else
  s.mem = 0;
@endif
@if _MPI
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
@else
  s.min = s.max = s.avg = 0.;
@endif
  s.speed = s.real > 0. ? tnc/s.real : -1;
  return s;
}

void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
	   "\n# " GRIDNAME 
	   ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
	   i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
@if _MPI
  fprintf (fout,
	   "# %d procs, MPI: min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   s.min, 100.*s.min/s.real,
	   s.avg, 100.*s.avg/s.real,
	   s.max, 100.*s.max/s.real);
@endif
}

typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
  foreach(reduction(max:max) reduction(+:avg) 
	  reduction(+:rms) reduction(+:volume) @) 
    if (f[] != nodata) {
      double v = fabs(f[]);
      if (v > max) max = v;
      volume += dv();
      avg    += dv()*v;
      rms    += dv()*sq(v);
    }
  norm n;
  n.avg = avg/volume;
  n.rms = sqrt(rms/volume);
  n.max = max;
  n.volume = volume;
  return n;
}

typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
	  reduction(max:max) reduction(min:min) @) 
    if (f[] != nodata) {
      volume += dv();
      sum    += dv()*f[];
      sum2   += dv()*sq(f[]);
      if (f[] > max) max = f[];
      if (f[] < min) min = f[];
    }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}

static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}

double theta = 1.3;

@if _QCCACC
/* Make parameter theta available on device */
ACC(acc declare copyin(theta))
/* Compile minmod2 for device */
ACC(acc routine seq)
@endif
double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}

@if _QCCACC
/*
   The OpenACC version fixes "gradient" to minmod2,
   function pointers cannot by used (yet) in PGI's
   OpenACC implementation.
*/
void gradients (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
  scalar s; vector v;
  for (s,v in f,g) {
    if (s.gradient)
      foreach(@) {
	foreach_dimension()
	  v.x[] = minmod2 (s[-1], s[], s[1])/Delta;
      }
    else // centered
      foreach(@) {
	foreach_dimension()
	  v.x[] = (s[1] - s[-1])/(2.*Delta);
      }
  }
  boundary ((scalar *) g);
}
@else
void gradients (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
  foreach() {
    scalar s; vector v;
    for (s,v in f,g) {
      if (s.gradient)
	foreach_dimension()
	  v.x[] = s.gradient (s[-1], s[], s[1])/Delta;
      else // centered
	foreach_dimension()
	  v.x[] = (s[1] - s[-1])/(2.*Delta);
    }
  }
  boundary ((scalar *) g);
}
@endif

/**
## Vorticity

For convenience, we also define a function which, given a velocity
field $\mathbf{u}$, fills a scalar field $\omega$ with the vorticity
field
$$
\omega = \partial_x u_y - \partial_y u_x
$$ */

void vorticity (const vector u, scalar omega)
{
  #if dimension > 1
    foreach()
      omega[] = (u.y[1,0] - u.y[-1,0] + u.x[0,-1] - u.x[0,1])/(2.*Delta);
    boundary ({omega});
  #endif
}

/**
Performance statistics are stored in this structure. */

struct {
  // total number of (leaf) cells advanced for this process
  long nc;
  // total number of (leaf) cells advanced for all processes
  long tnc;
  // real time elapsed since the start
  double t;
  // average computational speed (leaves/sec)
  double speed;
  // global timer
  timer gt;
} perf;

void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}

#include "output.h"
