#ifndef INTERPOLATE_H
#define INTERPOLATE_H

double PI;
void cart2pol(double x, double y, double *rho, double *phi);
void pol2cart(double rho, double phi, double *x, double *y);
double interp_lin(double x[], double y[], int Nx, int Ny, double xi, double yi, double F[]);

#endif
