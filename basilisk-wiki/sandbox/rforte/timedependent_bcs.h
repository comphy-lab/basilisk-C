/**
# Time-dependent inlet boundary conditions

This file contains:

- the `InletData` structure, 
- binary input loading,
- space-time interpolation functions,
- inlet velocity boundary conditions.
*/

typedef struct {
  int nt, ny, nz;
  double * times;
  double * y;
  double * z;
  double * Ux;
  double * Uy;
  double * Uz;
} InletData;

/**
Maps `(it, iy, iz)` to the corresponding index in the flattened array.
*/

#define IDX3(D,it,iy,iz) ((((size_t)(it))*(D)->ny + (iy))*(D)->nz + (iz))     // function like macro     

/**
Linearly interpolates between `a` and `b`.
The parameter `s` is the normalized interpolation coordinate:

- `s = 0` returns `a`,
- `s = 1` returns `b`.

*/

inline double lerp (double a, double b, double s) { return a + s*(b - a); }

inline double bilerp (double q00, double q10,  double q01, double q11,  double sy, double sz)
{ return (1. - sy) * (1. - sz) * q00 + sy * (1. - sz)*q10 + (1. - sy) * sz * q01 +  sy * sz * q11; }

/**
Returns the lower index `i` such that `a[i] <= q <= a[i+1]`.

The array `a` must be sorted in ascending order.
If `q` is outside the range, the first or last valid interval is returned.
*/

inline int lower_index (const double * a, int n, double q)
{
// a: ordered array of points
// q: point to localize
// n: dimension of a
// OUTPUT: index i s.t. a[i] <= q <= a[i+1] 

  if (q <= a[0])       return 0;        // return first interval
  if (q >= a[n - 1])   return n - 2;    // return last interval

  int lo = 0, hi = n - 1;
  while (hi - lo > 1) {                 // restric interval 
    int mid = (lo + hi)/2;
    if (a[mid] <= q) lo = mid;          // if q is on the right, move lower border
    else             hi = mid;          // if q is on the left move upper border
  }
  return lo;
}

/**
Loads the inlet dataset from a binary file.
The file is assumed to contain:

- the dimensions `nt`, `ny`, `nz`,
- the coordinate arrays `times`, `y`, and `z`,
- the velocity components `Ux`, `Uy`, and `Uz`.

Memory is allocated dynamically for all arrays stored in `InletData`.
The function also prints the dataset dimensions and coordinate ranges
for diagnostic purposes.
*/
int load_inlet_data (const char * fname, InletData * D)
{
  // Read binary file
  FILE * fp = fopen(fname, "rb");

  // Store integers
  fread(&D->nt, sizeof(int), 1, fp);
  fread(&D->ny, sizeof(int), 1, fp);
  fread(&D->nz, sizeof(int), 1, fp);
  fprintf(stderr, "Loading inlet data '%s': nt=%d ny=%d nz=%d\n", fname, D->nt, D->ny, D->nz);

  // Store arrays
  D->times = (double *) malloc(D->nt*sizeof(double));   
  D->y     = (double *) malloc(D->ny*sizeof(double));
  D->z     = (double *) malloc(D->nz*sizeof(double));

  // Total lenght
  size_t N = (size_t)D->nt*D->ny*D->nz;

  D->Ux = (double *) malloc(N*sizeof(double));
  D->Uy = (double *) malloc(N*sizeof(double));
  D->Uz = (double *) malloc(N*sizeof(double));

  fread(D->times, sizeof(double), D->nt, fp);
  fread(D->y,     sizeof(double), D->ny, fp);
  fread(D->z,     sizeof(double), D->nz, fp);
  fread(D->Ux, sizeof(double), N, fp);
  fread(D->Uy, sizeof(double), N, fp);
  fread(D->Uz, sizeof(double), N, fp);
  fclose(fp);

  // Print y,z,t ranges
  fprintf(stderr, "  y: [%g, %g]\n", D->y[0], D->y[D->ny - 1]);
  fprintf(stderr, "  z: [%g, %g]\n", D->z[0], D->z[D->nz - 1]);
  fprintf(stderr, "  t: [%g, %g]\n", D->times[0], D->times[D->nt - 1]);

  return 0;
}

/**
Samples one inlet velocity component at the query point `(tq, yq, zq)`.
The procedure is:

1. clamp the query coordinates to the dataset bounds,
2. locate the surrounding indices in time and space,
3. perform bilinear interpolation in the `(y,z)` plane at two
   consecutive time levels,
4. linearly interpolate the result in time.

This provides a smooth space-time reconstruction of the inlet data
from the discrete values stored in the binary file.
*/
double sample_inlet_component (const InletData * D, const double * U, double tq, double yq, double zq)
{
  // Limit values between max and min ranges of the Data
  tq = clamp(tq, D->times[0], D->times[D->nt - 1]);
  yq = clamp(yq, D->y[0], D->y[D->ny - 1]);
  zq = clamp(zq, D->z[0], D->z[D->nz - 1]);

  // Find lower indexes
  int it = lower_index(D->times, D->nt, tq);
  int iy = lower_index(D->y,     D->ny, yq);
  int iz = lower_index(D->z,     D->nz, zq);
  
  // Intervals widths
  double dtloc = D->times[it + 1] - D->times[it];
  double dyloc = D->y[iy + 1]     - D->y[iy];
  double dzloc = D->z[iz + 1]     - D->z[iz];
  
  // Compute local normalized coordinates
  double st = (dtloc > 0.) ? (tq - D->times[it]) / dtloc : 0.;
  double sy = (dyloc > 0.) ? (yq - D->y[iy])     / dyloc : 0.;
  double sz = (dzloc > 0.) ? (zq - D->z[iz])     / dzloc : 0.;

  // Extract points at t(i), at all 4 cell centers to do bilinear interpolation
  double u000 = U[IDX3(D,it,   iy,   iz)];
  double u100 = U[IDX3(D,it,   iy+1, iz)];
  double u010 = U[IDX3(D,it,   iy,   iz+1)];
  double u110 = U[IDX3(D,it,   iy+1, iz+1)];

  // Extract points at t(i + 1) to do bilinear interpolation
  double u001 = U[IDX3(D,it+1, iy,   iz)];
  double u101 = U[IDX3(D,it+1, iy+1, iz)];
  double u011 = U[IDX3(D,it+1, iy,   iz+1)];
  double u111 = U[IDX3(D,it+1, iy+1, iz+1)];

  // Bilinear interpolation in space
  double u_t0 = bilerp(u000, u100, u010, u110, sy, sz);         // time i
  double u_t1 = bilerp(u001, u101, u011, u111, sy, sz);         // time i+1

  // Liner interpolation in time
  return lerp(u_t0, u_t1, st);      
}
