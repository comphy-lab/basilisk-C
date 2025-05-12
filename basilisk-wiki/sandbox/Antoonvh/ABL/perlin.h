/**
# Generate Perlin's noise.

Here we generate a scalar map via Perlin's noise procedure. 
 */
double *gradp; //Perlin's gradients

#define DELTAP_x (L0/(double)(nx))
#define DELTAP_z (L0/(double)(nz))
#ifndef SMOOTHSTP
#define SMOOTHSTP(x) (6*cube(x)*sq(x) - 15*sq(sq(x)) + 10*cube(x)) //Perlin's function
#endif
#define interps(a0, a1, w) (a0*(1. - w) + a1*w)
#define IND(j,k) ((2*(j*nz + k)))

void init_perlin (int nx, int nz) {
  gradp = (double*)malloc((nx)*(nz)*sizeof(double)*2);
  for (int j = 0; j < nx; j++) {
    for (int k = 0; k < nz; k++) {
      double ang = noise()*pi;
      gradp[IND(j,k)]     = sin(ang);
      gradp[IND(j,k) + 1] = cos(ang);
    }
  }
}

double dotprdt (double xp, double xg, double zp, double zg, int j, int k, int nx, int nz) {
  double l = sqrt(sq(DELTAP_x) + sq(DELTAP_z));
  return (((xp - xg)*gradp[IND(j,k)] + (zp - zg)*gradp[IND(j,k) + 1])/l);
}

double perlin (double xp, double zp, int nx, int nz) {
  
  int j = (int)(((xp - X0)/DELTAP_x) + 0.5);
  int k = (int)(((zp - Z0)/DELTAP_z) + 0.5);
  double xl = X0 + ((double)j - 0.5)*DELTAP_x;
  double xr = xl + DELTAP_x;
  double zb = Z0 + ((double)k - 0.5)*DELTAP_z;
  double zt = zb + DELTAP_z;
  j--; k--;
  if (j < 0) 
    j += nx;
  if (k < 0)
    k += nz;
  double xw = (xp - xl)/DELTAP_x;
  double zw = (zp - zb)/DELTAP_z;
  double n0 = dotprdt (xp, xl, zp, zb, j, k, nx, nz);
  j++;
  if (j >= nx)
    j -= nx;
  double n1 = dotprdt (xp, xr, zp, zb, j, k, nx, nz);
  double ix0 = interps (n0, n1, SMOOTHSTP(xw));
  k++;
  if (k >= nz)
    k -= nz;
  n1 = dotprdt (xp, xr, zp, zt, j, k, nx, nz);
  j--;
  if (j < 0)
    j += nx;
  n0 = dotprdt (xp, xl, zp, zt, j, k, nx, nz);
  double ix1 = interps (n0, n1, SMOOTHSTP(xw));
  return interps (ix0, ix1, SMOOTHSTP(zw));
}
/**
## Test  

* [Generate Topographies](tp.c)
*/
