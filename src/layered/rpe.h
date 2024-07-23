/**
# Resting Potential Energy

See [Ilicak et al., 2012](#ilicak2012), appendix A and [Petersen et
al., 2015](#petersen2015), section 2.3. This should obviously be much
better documented than this.

Note that other more computationally-efficient options exist, based
for example on sampling of histograms of density.

~~~bib
@article{ilicak2012,
title = {Spurious dianeutral mixing and the role of momentum closure},
journal = {Ocean Modelling},
volume = {45-46},
pages = {37-58},
year = {2012},
issn = {1463-5003},
doi = {https://doi.org/10.1016/j.ocemod.2011.10.003},
url = {https://www.sciencedirect.com/science/article/pii/S1463500311001685},
author = {Mehmet Ilicak and Alistair J. Adcroft and Stephen M. Griffies and 
Robert W. Hallberg}
}

@article{petersen2015,
title = {Evaluation of the arbitrary {L}agrangian–{E}ulerian vertical coordinate 
         method in the {MPAS}-Ocean model},
journal = {Ocean Modelling},
volume = {86},
pages = {93-113},
year = {2015},
issn = {1463-5003},
doi = {https://doi.org/10.1016/j.ocemod.2014.12.004},
url = {https://www.sciencedirect.com/science/article/pii/S1463500314001796},
author = {Mark R. Petersen and Douglas W. Jacobsen and Todd D. Ringler and 
Matthew W. Hecht and Mathew E. Maltrud},
}
~~~
*/

#include "utils.h"
#if dimension == 2
# include "fractions.h"
#endif

struct _Zarea {
  scalar zb;  // compulsory
  double * A, * V; // compulsory
  int n;      // compulsory
  double max; // optional (default 0.)
  double min; // optional (default zb.min)
};
typedef struct _Zarea Zarea;

#define A(i) (a->A[i])

double area_integral (Zarea * a, double z1, double z2, double * Az)
{
  assert (z2 >= z1);
  double dz = (a->max - a->min)/(a->n - 1);

  int i2 = (z2 - a->min)/(a->max - a->min)*(a->n - 1);
  double f1, f2, z;
  if (i2 < 0)
    z = z2 = 0.;
  else if (i2 < a->n - 1) {
    f1 = A(i2), z = a->min + dz*i2;
    f2 = A(i2) + (z2 - z)/dz*(A(i2 + 1) - f1);
  }
  else {
    z = a->max, f1 = f2 = A(a->n - 1);
    i2 = a->n - 1;
  }
  double v2, zv2;
  if (z < z2) {
    v2 = (z2 - z)*(f2 + f1)/2.;
    zv2 = ((sq(z2) - sq(z))/2.*(f1*z2 - f2*z) +
	   (cube(z2) - cube(z))/3.*(f2 - f1))/(z2 - z);
  }
  else
    v2 = 0., zv2 = 0.;
  
  int i1;
  double v1, zv1;
  if (z1 <= a->min)
    v1 = 0., zv1 = 0., i1 = -1;
  else {
    i1 = (z1 - a->min)/(a->max - a->min)*(a->n - 1);
    double f1, f2, z;
    if (i1 >= a->n - 1) {
      f1 = f2 = A(a->n - 1);
      z = a->max;
    }
    else {
      f1 = A(i1) + (z1 - a->min - i1*dz)/dz*(A(i1 + 1) - A(i1));
      if (i2 > i1)
	f2 = A(i1 + 1), z = a->min + dz*(i1 + 1);
      else
	f2 = A(i1), z = a->min + dz*i1;
    }
    if (z != z1) {
      v1 = (z - z1)*(f2 + f1)/2.;
      zv1 = ((sq(z) - sq(z1))/2.*(f1*z - f2*z1) +
	     (cube(z) - cube(z1))/3.*(f2 - f1))/(z - z1);
    }
    else
      v1 = 0., zv1 = 0.;
  }
    
  double V = v1 + v2;
  *Az = zv1 + zv2;
  for (int i = i1 + 1; i < i2; i++) {
    z = a->min + i*dz, z2 = z + dz;
    f1 = A(i), f2 = A(i+1);
    V += (z2 - z)*(f2 + f1)/2.;
    *Az += ((sq(z2) - sq(z))/2.*(f1*z2 - f2*z) +
	    (cube(z2) - cube(z))/3.*(f2 - f1))/(z2 - z);
  }
  return V;
}

double area_z2 (Zarea * a, double z1, double dv, double * Az)
{
  double zm = z1, z2 = z1 + 10.*(a->max - a->min);
  double vol = area_integral (a, z1, z2, Az);
  assert (vol > dv);
  while (fabs (z2 - zm) > dry) {
    double z = (zm + z2)/2.;
    vol = area_integral (a, z1, z, Az);
    if (vol > dv) z2 = z;
    else zm = z;
    //    fprintf (stderr, "^^^ dv = %g: zm: %g z2: %g vol: %g\n", dv, zm, z2, vol);
  }
  //  assert (fabs(dv - vol) < 1e-6*dv);
  return (zm + z2)/2.;
}

Zarea zarea (scalar zb, double * A, double * V, int n,
	     double max = 0., double min = 0.)
{
  if (min == 0.)
    min = statsf (zb).min;
  for (int i = 0; i < n; i++) {
    double val = min + i*(max - min)/(n - 1);
#if dimension == 2    
    scalar f[];
    fractions (zb, f, val = val);
    stats sf = statsf(f);
    A[i] = sf.volume - sf.sum;
#else
    double area = 0.;
    foreach(reduction(+:area)) {
      double a = (zb[] + zb[-1])/2. - val, b = (zb[] + zb[1])/2. - val;
      if (a*b < 0.) {
	double x = a/(a - b);
	area += dv()*(a < 0. ? x : 1. - x);
      }
      else if (a < 0.)
	area += dv();
    }
    A[i] = area;
#endif
  }
  return (Zarea){ zb, A, V, n, max, min };
}

Zarea zvolume (scalar zb, double * A, double * V, int n,
	       double max = 0., double min = 0.)
{
  Zarea s = zarea (zb, A, V, n, max, min);
  double volume = 0.;
  for (int i = 0; i < n; i++) {
    double dz = (s.max - s.min)/(n - 1);
    volume += dz*A[i];
    V[i] = volume;
  }
  return s;
}

double volumez (Zarea s, double volume)
{
  if (volume <= 0.) return s.min;
  if (volume >= s.V[s.n-1]) return s.max;
  int a = 0, b = s.n - 1;
  while (b > a + 1) {
    int i = (a + b)/2;
    if (s.V[i] > volume) b = i;
    else a = i;
  }
  double c = (s.V[a] - volume)/(s.V[a] - s.V[b]);
  double dz = (s.max - s.min)/(s.n - 1);
  return (s.min + a*dz)*(1. - c) + (s.min + b*dz)*c;
}

double zareaval (Zarea s, double z)
{
  double dz = (s.max - s.min)/(s.n - 1);
  int j = (z - s.min)/dz;
  assert (j >= 0 && j < s.n - 1);
#if 0  
  printf ("A %d %g %g %g\n", j, s.A[j], s.A[j+1],
	  s.A[j] + (z - s.min)/dz*(s.A[j+1] - s.A[j]));
#endif
  return s.A[j] + (z - s.min)/dz*(s.A[j+1] - s.A[j]);
}

double barycenter (Zarea s, double z1, double z2)
{
  if (z1 == z2)
    return z1;
  assert (z1 >= s.min && z1 < z2 && z2 <= s.max);
  int n = 1;
  double dz = (z2 - z1)/n, sz = 0., sa = 0., z = z1 + dz/2.;
  for (int i = 0; i < n; i++) {
    double a = zareaval (s, z);
    sa += a, sz += z*a, z += dz;
  }
  //  printf ("# %g %g %g\n", z1, z2, sz/sa);
  return sz/sa;
}

typedef struct {
  double rho, dv, z;
} Parcel;

int heavier_than (const void * a, const void * b)
{
  Parcel * p1 = (Parcel *)a, * p2 = (Parcel *)b;
  return p1->rho > p2->rho ? -1 : 1;
}

double energy (double * PE = NULL, double * KE = NULL)
{
  double vPE = 0., vKE = 0.;
  foreach(reduction(+:vPE) reduction(+:vKE)) {
    double z = zb[];
    foreach_layer() {
      z += h[]/2;
      double mass = h[]*dv()*(1. + drho(T[]));
      vPE += mass*z;
#if NH
      vKE += mass*sq(w[])/2.;
#endif
      foreach_dimension()
	vKE += mass*sq(u.x[])/2.;
      z += h[]/2;
    }
  }
  vPE *= G;
  if (PE) *PE = vPE;
  if (KE) *KE = vKE;
  return vPE + vKE;
}

// Parallel sort

size_t psort (void ** base, size_t nmemb, size_t size,
	      int (* compar)(const void *, const void *))
{
#if _MPI // fixme: sorting is not done in parallel
  // Number of MPI processes and current rank
  int nproc, rank;
  MPI_Comm_size (MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  int counts[nproc];
  // Each process tells the root how many elements it holds
  nmemb *= size;
  MPI_Gather (&nmemb, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Displacements in the receive buffer for MPI_GATHERV
  int disps[nproc];
  // Displacement for the first chunk of data - 0
  for (int i = 0; i < nproc; i++)
    disps[i] = i > 0 ? (disps[i - 1] + counts[i - 1]) : 0;

  if (rank == 0) {
    size_t nmemb1 = disps[nproc-1] + counts[nproc-1];
    void * alldata = malloc (nmemb1);
    MPI_Gatherv (*base, nmemb, MPI_BYTE,
		 alldata, counts, disps, MPI_BYTE, 0, MPI_COMM_WORLD);
    free (*base);
    *base = alldata;
    nmemb = nmemb1/size;
    qsort (*base, nmemb, size, compar);
  }
  else
    MPI_Gatherv (*base, nmemb, MPI_BYTE,
		 NULL, counts, disps, MPI_BYTE, 0, MPI_COMM_WORLD);
#else
  qsort (*base, nmemb, size, compar);
#endif
  return nmemb;
}

// Resting Potential Energy: see e.g. Ilicak et al, 2012, Appendix A

trace
double RPE (int n = 100)
{
  double A[n], V[n];
  Zarea vol = zarea (zb, A, V, n);
#if 0
  for (int i = 0; i < vol.n; i++) {
    double dz = (vol.max - vol.min)/(vol.n - 1);
    fprintf (stderr, "%g %g %g\n", vol.min + i*dz, V[i], - statsf (zb).sum);
  }

  for (double volume = 0.; volume <= V[vol.n-1]; volume += V[vol.n-1]/35.)
    printf ("%g %g\n", volumez (vol, volume), volume);
  exit (0);
#endif
  
  long nb = 0;
  foreach_leaf() nb++;
  nb *= nl;
  Parcel * p = malloc (sizeof(Parcel)*nb);
  long i = 0;
  foreach_leaf()
    foreach_layer() {
      assert (i < nb);
      p[i++] = (Parcel){ 1. + drho(T[]), h[]*dv() };
    }
  
  nb = psort ((void **) &p, nb, sizeof (Parcel), heavier_than);

  double RPE = HUGE;
  if (pid() == 0) {
    RPE = 0.;
    double volume = 0.;
    double z = vol.min;
    for (int i = 0; i < nb; i++)
      if (p[i].dv > 0.) {
	volume += p[i].dv;
#if 1
	double Az, zt = area_z2 (&vol, z, p[i].dv, &Az);
	p[i].z = Az/p[i].dv;
	RPE += p[i].rho*Az; // see e.g. (1) in Ilicak et al, 2012
#else
	double zt = volumez (vol, volume);
	p[i].z = barycenter (vol, z, zt);
	
	//	fprintf (stderr, "&& %g %g %g %g\n", zt, zt1, p[i].z, Az/p[i].dv);

	p[i].z = Az/p[i].dv;
	RPE += p[i].rho*p[i].z*p[i].dv; // see e.g. (1) in Ilicak et al, 2012
#endif
	
	//	fprintf (stdout, "%g %g %g %g\n", p[i].rho, p[i].dv, p[i].z, zt);
	z = zt;	
      }
  }
  //  fprintf (stdout, "****\n");
  mpi_all_reduce (RPE, MPI_DOUBLE, MPI_MIN);
  free (p);
  
  return RPE*G;
}
