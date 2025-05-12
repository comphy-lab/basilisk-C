/**
# Example: Planetary geostrophic model

This is the driver file for pg.h.  We define the grid, the topography,
the initial conditions and the output routines.

uncomment #define H5FILE_NAME "fields.h5" to output fields in hdf5 format
*/

//#define H5FILE_NAME "fields.h5"
#include "grid/multigrid.h"
#include "pg.h"

/**
 define topography: h is the topography, hx and hy the x and y
 derivatives */

double d = 0.1;
double h0 = 0.05;
#define g(x)    ( 1-exp(-sq(x)/(2*sq(d))) + h0)
#define gp(x)   ( (x)/sq(d)*exp(-sq(x)/(2*sq(d))))

double h (double x, double y){
  return ( g(x-0.0)*g(1.0-x)*g(y-ys)*g(1.0-y+ys) );
}

double hx (double x, double y){
  return ( (gp(x-0.0)*g(1.0-x)-g(x-0.0)*gp(1.0-x))*g(y-ys)*g(1.0-y+ys));
}

double hy (double x, double y){
  return ( g(x-0.0)*g(1.0-x)*(gp(y-ys)*g(1.0-y+ys)-g(y-ys)*gp(1.0-y+ys)));
}

/* double h  (double x, double y){ return ( 1.0);} */
/* double hx (double x, double y){ return ( 0.0);} */
/* double hy (double x, double y){ return ( 0.0);} */


/**
 spatially varying non dimensional diffusivity coef. in dimensional units
$$
\kappa^* = \kappa \frac{N^2H^4}{\beta L^3}
$$
 */

double k (double x, double y, double s) { return (1e-2);}

/**
   wind stress and wind stress derivative
*/

double tau0 = 2e-2;

double taux   (double x, double y){ return (tau0*y*sin(2*(y-ys)*pi));}
double taux_y (double x, double y){ return (2*pi*tau0*y*cos(2*(y-ys)*pi) + tau0*sin(2*(y-ys)*pi));}
/* double taux   (double x, double y){ return (0.);} */
/* double taux_y (double x, double y){ return (0.);} */
double tauy   (double x, double y){ return (0.);}
double tauy_x (double x, double y){ return (0.);}


int main() {
  N = 64;
  nl = 20;

  CFL = 0.02;
  DT = 1.e-2;
  tend = 0.2;

  /**
     physical parameters : $a$ is the aspect ratio of the basin $H/L$.
     $r$ is the non dimensional friction coefficient. The dimensional
     friction is $r^* = r \beta L $. tau_surf is the surface
     relaxation time scale. */

  a = 0.7; 
  r = 0.1; // idealy, should be ~0.02 but BTsolver doesn't like it.
  tau_surf = 5e-2;

  /**
     We use a pseudo SOR impletation because for small r, the
     ellicptic system is not diagonally dominant. For omega = 1, we
     recover the default basilisk implementation. omega<1 slows down
     the convergence rate, but at least it converges. */

  omega = 0.2;

  ys = 0.2;
  origin (0.0, ys);

//  outdir = "outdir/";

  run(); 
}


/**
   Initial conditions and surface forcing (keep the name of this event
   'init' so that BC are set automatically: see pg.h init event)*/
event init (t = 0) {
  for (int l = 0; l < nl; l++) {
    scalar b = bl[l+1];
    foreach() 
      b[] = 0.0;
  }
  boundary (bl);

  foreach() 
    b_surf[] = 2*cos(pi*(y-ys));

/**
   write topo in file
 */

  char name[80];
  sprintf (name, "%s%s", outdir, "topo.dat");
  FILE * fp = fopen (name, "w");
  output_field ({hc}, fp);
  fclose(fp);

}

// stdout
event writestdout (i++) {
  /* omega = 1.1*omega; */
  /* if (omega > 0.2) omega = 0.2; */
  /* printf (" omega = %g\n", omega); */
  printf ("i = %i, dt = %g, t = %g\n", i, dt, t);
}


/* event adapt (i++) { */
/*  astats s = adapt_wavelet (bl, (double []){1.0}, maxlevel = 10); */
/*  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc); */
/* } */


/**
   write output
 */
event writestate (t = 0; t <= tend+1e-10;  t +=0.2) {
  printf ("i = %d t = %g\n", i, t);
  
  char name[80];
  sprintf (name,"%sstate%09d.dat", outdir, i);
  FILE * fp = fopen (name, "w");
  
  scalar * lout = list_copy(bl);
  scalar * lout1 = list_concat(lout,(scalar *) ul); free(lout);
  lout = list_concat(lout1, wl); free(lout1);
  lout = list_append(lout,psibt);
  output_field (lout, fp);
  free (lout);
  fclose(fp);

#ifdef H5FILE_NAME
  backup_fields (bl, ul, i)
#endif

}

/**
   We plot the buoyancy and velocity fields on a sigma surface.
   
   ~~~pythonplot Buoyancy and velocity at $\sigma = -1$
   import numpy as np
   import matplotlib.pyplot as plt
   import glob

   plt.ion()
   dir0 = "./"

   file1 = 'topo.dat'
   file2 = 'state*'

   allfiles2 = sorted(glob.glob(dir0 + file2));
   nb_files  = len(allfiles2);

   data = np.loadtxt(allfiles2[-1])
   N1, nl1 = data.shape
   N1 = int(np.sqrt(N1))
   nl = int((nl1 - 6)/4) # 5 extra rows: x,y,b(+2),w(+1),psibt
   N = N1 - 1
   
   # build vertical grid
   topo = np.loadtxt(dir0 + file1)
   topo = topo[:,2].reshape((N1,N1,1))
   ds = 1./nl
   sc = np.linspace(-1+0.5*ds,-0.5*ds,nl)
   z = topo*sc.reshape((1,1,nl))

   x = data[:,0].reshape((N1,N1))
   y = data[:,1].reshape((N1,N1))
   b = data[:,3:nl+3].reshape((N1,N1,nl))
   uv = data[:,nl+4:3*nl+4].reshape((N1,N1,2*nl))
   u = uv[:,:,::2]
   v = uv[:,:,1::2]
   w = data[:,3*nl+4:4*nl+5].reshape((N1,N1,nl+1))
   psibt = data[:,4*nl+5].reshape((N1,N1))
   
   # choose vertical level (0=bottom, nl-1=top)
   l = nl-1
   
   #plot velocity vector every nsk pts
   nsk = 5

   plt.figure()
   plt.pcolormesh(x,y,b[:,:,l])
   plt.colorbar()
   plt.contour(x,y,topo[:,:,0],colors='w',linewidths=0.5)
   plt.contour(x,y,psibt,colors='r',linewidths=0.5)
   plt.quiver(x[::nsk,::nsk],y[::nsk,::nsk],u[::nsk,::nsk,l],v[::nsk,::nsk,l])
   plt.xlabel('x')
   plt.ylabel('y')
   plt.title('Buoyancy (color), velocity (arrows), topography (contours), and BT stream function (red)')
   plt.savefig('xyplot.png')
   ~~~

   and on a vertical section

   ~~~pythonplot Buoyancy and velocity at $y = L/2$
   plt.figure()
   xsec = np.tile(x[:,0],(nl,1)).T
   ysec = np.tile(y[0,:],(nl,1)).T
   iy = int(N/2)-1
   psi = v[:,iy,:]
   vmax = np.max(np.abs(psi))
   psi[0,0] = vmax
   psi[1,0] = -vmax
   plt.contour(xsec,z[:,iy,:],b[:,iy,:],colors="k",linewidths=1.0)
   plt.contourf(xsec,z[:,iy,:],psi,30,cmap=plt.cm.bwr)
   plt.plot(x[:,iy],z[:,iy,0],"k",linewidth=1.5)
   plt.fill_between(x[:,iy], z[:,iy,0],np.min(z[:,iy,0]) + 0*z[:,iy,0],facecolor="0.7")
   plt.colorbar()
   plt.xlabel("x")
   plt.ylabel("z")
   plt.title('V vel (color), buoyancy (contours)')
   plt.savefig('xzplot.png')
   ~~~
*/
