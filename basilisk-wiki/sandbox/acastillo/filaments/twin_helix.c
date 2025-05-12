/**
# Simulated Navier-Stokes closely-spaced helical vortex pair.

![The movie shows a $\lambda_2=0$ isosurface on top of a slice of $u_z$](twin_helix/twin_helix2.mp4)

![A view from the side showing the rotation of the vortex pair. Not much is happening yet.](twin_helix/twin_helix3.mp4)

*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "PointTriangle.h"

int n_seg;
int minlevel = 4;
int maxlevel = 10;
double as;
double Hs1, Rs1, n_turns1;
double Hs2, Rs2, n_turns2;

int main() {
  L0 = 42;
  X0 = Y0 = Z0 = -L0/2;
  DT = 0.01;
  N = 1<<(minlevel+maxlevel)/2;
  periodic(top);
  periodic(left);
  periodic(back);
  run();
}

event adapt (i++){
 adapt_wavelet ((scalar*){u}, (double[]){1e-4, 1e-4, 1e-4}, maxlevel, minlevel);
}

/**
## Initialisation
For the initialisation process, we proceed as follows.

1. The position of the vortex segments is given by a space-curve $c$
parametrized as function of t0.

2. We require the first and second derivatives
of $c$, as to define a local Serret-Frenet frame $(\vec{t},\vec{n},\vec{b})$.
$$
\vec{t} = \frac{d\vec{c}}{dt_0}/\|\frac{d\vec{c}}{dt_0}\|
$$
$$
\vec{n} = \frac{d^2\vec{c}}{dt_0^2}/\|\frac{d^2\vec{c}}{dt_0^2}\|
$$
$$
\vec{b} = \vec{t} \times \vec{n}
$$

3. Each position $\vec{p}$ is projected into the local Frenet-Serret frame
to obtain a set of local coordinates, such that:
$$
(\vec{p} - \vec{c}(t_0)) \cdot \vec{t} = 0
$$
$$
(\vec{p} - \vec{c}(t_0)) \cdot \vec{n} = x_n
$$
$$
(\vec{p} - \vec{c}(t_0)) \cdot \vec{b} = x_b
$$
which is done through a minization process.

4. We use the local coordinates $(x_n, x_b)$ to define a radial coordinate
$\rho$ required to compute the vorticity of a Lamb-Oseen vortex as
$$
\vec{\omega} = \Gamma/(\pi a^2) \exp(-\rho^2/a^2) \cdot \vec{t}
$$
where $\Gamma$ is the circulation and $a$ the core size.

5. We consider the flow is expressed as a vector potential,
$$
\vec{u} = \nabla \times \vec{\psi}
$$
such that
$$
\nabla^2 \vec{\psi} = -\vec{\omega}
$$
which requires solving a Poisson problem for each vorticity component.

6. We use the vector potential to evaluate the velocity field.
*/
#include "view.h"
#include "filaments.h"
#include "filaments.c"
#include "lambda2.h"
#include "../output_fields/output_vtu_foreach.h"

#define RAD (sqrt(sq(x) + sq(y)))
event init (t = 0) {
  FILE * fp ;
  int restoring = 1;
  if (restoring){
    printf( "Restoring from previous run \n" );
    restore("dump_init_t25");
  }
  else {
    printf( "Initial conditions from filamentary solutions \n" );
    // Step 1a and 2a: Load external vortex
    fp = fopen("hext.asc", "r");
    fscanf(fp, "%d",  &n_seg);
    fscanf(fp, "%lf", &Rs1);
    fscanf(fp, "%lf", &as);
    fscanf(fp, "%lf", &Hs1);

    n_turns1 = L0/Hs1;
    double t0[n_seg];
    coord c1[n_seg], dc1[n_seg], d2c1[n_seg];

    for (int i = 0; i < n_seg; i++)
      fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t0[i], &c1[i].x, &c1[i].y, &c1[i].z, &dc1[i].x, &dc1[i].y, &dc1[i].z, &d2c1[i].x, &d2c1[i].y, &d2c1[i].z);
    fclose (fp);

    // Step 1b and 2b: Load internal vortex
    fp = fopen("hint.asc", "r");
    fscanf(fp, "%d",  &n_seg);
    fscanf(fp, "%lf", &Rs2);
    fscanf(fp, "%lf", &as);
    fscanf(fp, "%lf", &Hs2);

    n_turns2 = L0/Hs2;
    coord c2[n_seg], dc2[n_seg], d2c2[n_seg];

    for (int i = 0; i < n_seg; i++)
      fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t0[i], &c2[i].x, &c2[i].y, &c2[i].z, &dc2[i].x, &dc2[i].y, &dc2[i].z, &d2c2[i].x, &d2c2[i].y, &d2c2[i].z);
    fclose (fp);

    // Refine aroud the region of interest
    refine  (RAD < Rs1 + 16*as && RAD > Rs2 - 16*as && level < maxlevel - 2);
    refine  (RAD < Rs1 +  8*as && RAD > Rs2 -  8*as && level < maxlevel - 1);
    refine  (RAD < Rs1 +  4*as && RAD > Rs2 -  4*as && level < maxlevel);

    cells();
    save ("helical_pair0.png");

    view (theta = pi/6, phi = pi/6, fov=30);
    box();
    cells();
    draw_space_curve(n_seg, c1);
    draw_space_curve(n_seg, c2);
    save ("helical_pair1.png");

    // Compute the local Frenet-Serret frame for each curve
    coord tvec1[n_seg], nvec1[n_seg], bvec1[n_seg];
    coord tvec2[n_seg], nvec2[n_seg], bvec2[n_seg];
    for (int i = 0; i < n_seg; i++){
      foreach_dimension(){
        tvec1[i].x =  dc1[i].x/sqrt(vecdot( dc1[i],  dc1[i]));
        nvec1[i].x = d2c1[i].x/sqrt(vecdot(d2c1[i], d2c1[i]));

        tvec2[i].x =  dc2[i].x/sqrt(vecdot( dc2[i],  dc2[i]));
        nvec2[i].x = d2c2[i].x/sqrt(vecdot(d2c2[i], d2c2[i]));
      }
      bvec1[i] = vecdotproduct(tvec1[i], nvec1[i]);
      bvec2[i] = vecdotproduct(tvec2[i], nvec2[i]);
    }

    // Steps 3 and 4: Local projection and evaluate vorticity
    vector omega[];
    foreach(){
      coord pcar = {x,y,z};
      coord val_omega1 = get_vorticity_filament(pcar, n_seg, as, t0, c1, tvec1, nvec1, bvec1, 0);
      coord val_omega2 = get_vorticity_filament(pcar, n_seg, as, t0, c2, tvec2, nvec2, bvec2, 0);
      foreach_dimension()
        omega.x[] = -val_omega1.x - 0.95*val_omega2.x;
    }

    // Step 5: Solve a Poisson problem for each vorticity component
    vector psi[];
    foreach()
      foreach_dimension()
        psi.x[] = 0.;
    boundary ((scalar*){psi,omega});
    poisson (psi.x, omega.x);
    poisson (psi.y, omega.y);
    poisson (psi.z, omega.z);

    // Step 6: Compute the velocity from the vector potential field
    foreach(){
      u.x[] = ((psi.z[0,1,0] - psi.z[0,-1,0]) - (psi.y[0,0,1] - psi.y[0,0,-1]))/(2.*Delta);
      u.y[] = ((psi.x[0,0,1] - psi.x[0,0,-1]) - (psi.z[1,0,0] - psi.z[-1,0,0]))/(2.*Delta);
      u.z[] = ((psi.y[1,0,0] - psi.y[-1,0,0]) - (psi.x[0,1,0] - psi.x[0,-1,0]))/(2.*Delta);
    }
    boundary ((scalar *){u});

    adapt_wavelet ((scalar*){u}, (double[]){1e-4, 1e-4, 1e-4}, maxlevel, minlevel);
  }
  printf( "Initial conditions ... Done! \n" );

  double reynolds= 16100;
  const face vector muc[] = {1./reynolds,1./reynolds,1./reynolds};
  mu = muc;

  scalar l2[];
  lambda2 (u, l2);

  view (theta = pi/6, phi = pi/6, fov=30);
  box();
  isosurface ("l2", 0);
  save ("helical_pair2.png");

  fp = fopen("vortex_x0.asc", "w");
  fputs ("[1]t\t [2]tag\t [3]Gamma\t [4]mu_x\t [5]mu_y\t [6]mu_z\t [7]M20\t [8]M02\t [9]M11\t [10]a\t [11]b\t [12]c\t [13]e \n", fp);
  fclose(fp);

  fp = fopen("vortex_y0.asc", "w");
  fputs ("[1]t\t [2]tag\t [3]Gamma\t [4]mu_x\t [5]mu_y\t [6]mu_z\t [7]M20\t [8]M02\t [9]M11\t [10]a\t [11]b\t [12]c\t [13]e \n", fp);
  fclose(fp);

  fp = fopen("twin_helix.asc", "w");
  fputs ("[1]t\t [2]ekin\t [3]enstrophy\t [4]helicity \n", fp);
  fclose (fp);
}

/**
## Ellipticity of the vortex cores
The goal is to evaluate several quantities that are pertinent to describe the
flow field such as circulation, sizes, and ellipticity. When a vorticity field
possesses an elliptical shape, these quantities are sufficient to define its
geometry. In this case, we measure the ellipticity from the intersection of
each vortex with the plane $z=0$.
We also track the evolution of the kinetic energy, the viscous
enery dissipation rate, as well as the total helicity.
*/

#include "3d/ellipticity.h"
event logfile (t+=0.02) {

  vector omega[];
  vorticity3d(u, omega);

  scalar omega_mag[];
  foreach()
    omega_mag[] = sqrt(sq(omega.x[]) + sq(omega.y[]) + sq(omega.z[]));

  scalar l2[];
  lambda2 (u, l2);

  scalar m1[], m2[];
  foreach(){
    m1[] = ((l2[] < -0.001) & (abs(z-4) < 5.25/4.0)) & (y > 0);
    m2[] = ((l2[] < -0.001) & (abs(z-4) < 5.25/4.0)) & (x > 0);
  }

  int n1 = tag(m1);
  int n2 = tag(m2);

  FILE * fp = fopen("vortex_x0.asc", "a");
  vorticity_moments_plane(omega_mag, m1, fp, (coord){1,0,0}, 0.);
  fclose(fp);

  fp = fopen("vortex_y0.asc", "a");
  vorticity_moments_plane(omega_mag, m2, fp, (coord){0,1,0}, 0.);
  fclose(fp);

  double kin=0., ens=0., hel=0.;
  foreach(reduction(+:kin), reduction(+:ens), reduction(+:hel)){
    foreach_dimension(){
      kin += dv() * sq(u.x[]);
      ens += dv() * sq(omega.x[]);
      hel += dv() * u.x[] * omega.x[];
    }
  }

  fp = fopen("twin_helix.asc", "a");
  fprintf (fp, "%f %.9g %.9g %.9g \n", t, kin, ens, hel);
  fclose (fp);
}

/**
We also create a movie with isocontours of $\lambda_2=0$ superimposed to
different quantities
*/

event movie (t += 0.02) {
  scalar l2[];
  lambda2 (u, l2);
  foreach()
    l2[] = l2[] > -0.01 ? nodata : l2[];
  stats f = statsf (l2);
  foreach()
    l2[] = l2[] == nodata ? 0 : l2[];
  boundary ({l2});

  vector omega[];
  vorticity3d(u, omega);

  view (theta = pi/6, phi = pi/6, fov=30);
  isosurface ("l2", -f.stddev);
  save ("twin_helix1.mp4");

  isosurface ("l2", -f.stddev);
  squares ("u.z", linear = false, alpha=-L0/2);
  save ("twin_helix2.mp4");

  view (camera="left", fov=20);
  squares ("omega.x", linear = false, n= {1,0,0});
  save ("twin_helix3.mp4");
}

event snapshot (t = 30){
  dump();
}

/**
We extract a slice of the plane X=0 every 0.5 t.u., and export to Paraview.
These should serve for post-processing and validation. If `vtk` is available we
can plot these slices directly in `Python` (.vtu files are not uploaded to  the
wiki). There is still work to do for the tracking of the vortices and follow the
relaxation process, and spatial resolution could be refined, but it's a starting
point.

~~~pythonplot Slice of the normal velocity and vorticity components at the plane $z=0$ for $t=5$.
import numpy as np
import os
from vtk import vtkXMLUnstructuredGridReader
from vtk.util import numpy_support
from vtk.numpy_interface import dataset_adapter
import matplotlib.pyplot as plt

params = {'backend': 'ps',
          'axes.labelsize': 12,
          'font.size': 12,
          'legend.fontsize': 12,
          'font.size': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'font.family': 'serif',
          'figure.figsize': (3.75*3, 3.75*6/5),
          'contour.negative_linestyle':'dashed'}

plt.rcParams.update(params)

path = "./vtu/twin_helix_z0_t5.vtu"
reader = vtkXMLUnstructuredGridReader()
reader.SetFileName(path)
reader.Update()
data = reader.GetOutput()
points = numpy_support.vtk_to_numpy(data.GetPoints().GetData())
cells = numpy_support.vtk_to_numpy(data.GetCells().GetData())
ntri = cells.size//5
tri = np.take(cells,[n for n in range(cells.size) if n%5 != 0]).reshape(ntri,4)


cell_centers = np.zeros([ntri,3], float)
for i in range(ntri):
  cell_centers[i,0] = np.sum(points[tri[i,0:4],0])/4
  cell_centers[i,1] = np.sum(points[tri[i,0:4],1])/4
  cell_centers[i,2] = np.sum(points[tri[i,0:4],2])/4

ds = dataset_adapter.WrapDataObject(reader.GetOutput())
lambda2   = numpy_support.vtk_to_numpy(ds.CellData['l2'])
vorticity = numpy_support.vtk_to_numpy(ds.CellData['omega.x'])
velocity  = numpy_support.vtk_to_numpy(ds.CellData['u.x'])

f, axes = plt.subplots(1,2, sharex=True, sharey=False)
ax= np.ravel(axes)
cntr1 = ax[0].tricontourf(cell_centers[:,0], cell_centers[:,1], velocity[:,2],  levels=51, cmap=plt.cm.jet, vmin=-1.0, vmax=1.0)
cntr2 = ax[1].tricontourf(cell_centers[:,0], cell_centers[:,1], vorticity[:,2], levels=51, cmap=plt.cm.PuOr_r, vmin=-2.5, vmax=2.5)

clb1 = f.colorbar(cntr1, ax=ax[0], ticks=np.arange(-6,8,0.25))
clb2 = f.colorbar(cntr2, ax=ax[1], ticks=np.arange(-6,8,0.5))

for i in range(2):
  ax[i].set_yticks(np.arange(-20,20,5))
  ax[i].set_xlim([-10,10])
  ax[i].set_ylim([-10,10])
  ax[i].set_xlabel('$x$')
  ax[i].set_ylabel('$y$')
  ax[i].set_aspect('equal', 'box')

for c in cntr1.collections:
  c.set_edgecolor("face")

for c in cntr2.collections:
  c.set_edgecolor("face")


plt.tight_layout()
plt.savefig('./slice_z0.png', dpi=100)
~~~


~~~pythonplot Slice of the normal velocity and vorticity components at the plane $x=0$ for $t=5$.
import numpy as np
import os
from vtk import vtkXMLUnstructuredGridReader
from vtk.util import numpy_support
from vtk.numpy_interface import dataset_adapter
import matplotlib.pyplot as plt

params = {'backend': 'ps',
          'axes.labelsize': 12,
          'font.size': 12,
          'legend.fontsize': 12,
          'font.size': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'font.family': 'serif',
          'figure.figsize': (3.75*3, 3.75*6/5),
          'contour.negative_linestyle':'dashed'}

plt.rcParams.update(params)

path = "./vtu/twin_helix_x0_t5.vtu"
reader = vtkXMLUnstructuredGridReader()
reader.SetFileName(path)
reader.Update()
data = reader.GetOutput()
points = numpy_support.vtk_to_numpy(data.GetPoints().GetData())
cells = numpy_support.vtk_to_numpy(data.GetCells().GetData())
ntri = cells.size//5
tri = np.take(cells,[n for n in range(cells.size) if n%5 != 0]).reshape(ntri,4)

cell_centers = np.zeros([ntri,3], float)
for i in range(ntri):
  cell_centers[i,0] = np.sum(points[tri[i,0:4],0])/4
  cell_centers[i,1] = np.sum(points[tri[i,0:4],1])/4
  cell_centers[i,2] = np.sum(points[tri[i,0:4],2])/4

ds = dataset_adapter.WrapDataObject(reader.GetOutput())
lambda2   = numpy_support.vtk_to_numpy(ds.CellData['l2'])
vorticity = numpy_support.vtk_to_numpy(ds.CellData['omega.x'])
velocity  = numpy_support.vtk_to_numpy(ds.CellData['u.x'])

f, axes = plt.subplots(1,2, sharex=True, sharey=False)
ax= np.ravel(axes)
cntr1 = ax[0].tricontourf(cell_centers[:,2], cell_centers[:,1], velocity[:,0],  levels=51, cmap=plt.cm.jet, vmin=-0.15, vmax=0.15)
cntr2 = ax[1].tricontourf(cell_centers[:,2], cell_centers[:,1], vorticity[:,0], levels=51, cmap=plt.cm.PuOr_r, vmin=-10, vmax=10)

for c in cntr1.collections:
  c.set_edgecolor("face")

for c in cntr2.collections:
  c.set_edgecolor("face")

for i in range(2):
  ax[i].set_xlim([0,10])
  ax[i].set_ylim([2.0,12.0])
  ax[i].set_xlabel('$z$')
  ax[i].set_ylabel('$y$')

clb1 = f.colorbar(cntr1, ax=ax[0], ticks=np.arange(-0.2,0.2,0.05))
clb2 = f.colorbar(cntr2, ax=ax[1], ticks=np.arange(-15,15,5))

plt.tight_layout()
plt.savefig('./slice_x0.png', dpi=100)
~~~


~~~pythonplot Tracking the vortex centers inside a window between $z=-2$ and $z=2$.
import numpy as np
import os
from vtk import vtkXMLUnstructuredGridReader
from vtk.util import numpy_support
from vtk.numpy_interface import dataset_adapter
import matplotlib.pyplot as plt

def get_from_vtk(path):
	reader = vtkXMLUnstructuredGridReader()
	reader.SetFileName(path)
	reader.Update()
	data = reader.GetOutput()
	points = numpy_support.vtk_to_numpy(data.GetPoints().GetData())
	cells = numpy_support.vtk_to_numpy(data.GetCells().GetData())
	ntri = cells.size//5
	tri = np.take(cells,[n for n in range(cells.size) if n%5 != 0]).reshape(ntri,4)

	cell_centers = np.zeros([ntri,3], float)
	for i in range(ntri):
	  cell_centers[i,0] = np.sum(points[tri[i,0:4],0])/4
	  cell_centers[i,1] = np.sum(points[tri[i,0:4],1])/4
	  cell_centers[i,2] = np.sum(points[tri[i,0:4],2])/4

	ds = dataset_adapter.WrapDataObject(reader.GetOutput())
	lambda2   = numpy_support.vtk_to_numpy(ds.CellData['l2'])
	vorticity = numpy_support.vtk_to_numpy(ds.CellData['omega.x'])
	velocity  = numpy_support.vtk_to_numpy(ds.CellData['u.x'])
	return cell_centers, lambda2, vorticity, velocity



params = {'backend': 'ps',
          'axes.labelsize': 12,
          'font.size': 12,
          'legend.fontsize': 12,
          'font.size': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'font.family': 'serif',
          'figure.figsize': (3.75*15/5, 3.75*7.5/5),
          'contour.negative_linestyle':'dashed'}

plt.rcParams.update(params)

f, axes = plt.subplots(1,2, sharex=False, sharey=False)
ax= np.ravel(axes)



for j,i in enumerate(range(0,25,2)):
	cell_centers, lambda2, vorticity, velocity = get_from_vtk("./vtu/twin_helix_y0_t{0:g}.vtu".format(i))
	ax[0].tricontour(cell_centers[:,2], cell_centers[:,0], lambda2[:], levels=[-10,-1,-0.1,-0.01,0], colors='k', linewidths=0.5, alpha=0.05*j+0.1)

# Read the time series and separate each vortex using the minimal distance between
# consecutive time steps
file1 = 'vortex_y0.asc'
data = np.loadtxt(file1, skiprows=1)
data = data[~np.isnan(data).any(axis=1)]
dt = 0.02
t = 0.02
from_index = np.zeros([np.shape(data)[0],], int)
to_index = np.zeros([np.shape(data)[0],], int)

for it, t in enumerate(np.arange(t, t + 25, dt)):
    index1 = np.where( np.abs(data[:,0] - t) < 1e-8 )[0]
    for i in index1:
        mindist = 1e30
        index2 = np.where( np.abs(data[:,0] - (t+dt)) < 1e-8 )[0]
        for j in index2:
            dist = np.sqrt((data[i,3]-data[j,3])**2 + (data[i,4]-data[j,4])**2 + (data[i,5]-data[j,5])**2)

            if (dist < mindist):
                mindist = dist
                from_index[i] = i
                to_index[i]   = j

index = 0
ix1 = 0
while 1 :
    current = to_index[index]
    index = np.where(from_index[:] == current)[0]
    if index.size == 0:
        break
    ix1 = np.append(ix1, index[0])

index = 1
ix2 = 1
while 1 :
    current = to_index[index]
    index = np.where(from_index[:] == current)[0]
    if index.size == 0:
        break
    ix2 = np.append(ix2, index[0])

vortex1 = data[ix1,:]
vortex2 = data[ix2,:]

ax[0].set_xlim([-0.5,4.5])
ax[0].set_ylim([5.0,9.0])

ax[0].scatter(vortex1[::5,5], vortex1[::5,3], color='C0', s=3)
ax[0].scatter(vortex2[::5,5], vortex2[::5,3], color='C1', s=3)

ax[0].set_xlabel('$z$')
ax[0].set_ylabel('$y$')

ax[1].scatter(vortex1[:,0], vortex1[:,9], color='C0', s=3)
ax[1].scatter(vortex2[:,0], vortex2[:,9], color='C1', s=3)

ax[1].set_xlabel('$t$')
ax[1].set_ylabel('$a$')
ax[1].set_xlim([0,25])
ax[1].set_ylim([0,0.20])

plt.tight_layout()
plt.savefig('tracking.png')
~~~

*/

event export_to_paraview (t += 1.0){
  printf( "Exporting slices to VTU format \n" );
  fprintf (stdout, "%d \n", i);
  scalar l2[];
  lambda2 (u, l2);

  vector omega[];
  vorticity3d(u, omega);

  scalar m[];
  foreach()
    m[] = ((l2[] < -0.001) & (abs(z)< 5.25/4)) & (x > 0);

  int n = tag (m);

  char name[80];
  FILE * fp ;
  sprintf(name, "twin_helix_x0_t%g.vtu", t);
  fp = fopen(name, "w"); output_vtu_plane ((scalar *) {l2, m}, (vector *) {u, omega}, fp, (coord){1,0,0}, 0.0); fclose (fp);

  sprintf(name, "twin_helix_y0_t%g.vtu", t);
  fp = fopen(name, "w"); output_vtu_plane ((scalar *) {l2, m}, (vector *) {u, omega}, fp, (coord){0,1,0}, 0.0); fclose (fp);

  sprintf(name, "twin_helix_z0_t%g.vtu", t);
  fp = fopen(name, "w"); output_vtu_plane ((scalar *) {l2, m}, (vector *) {u, omega}, fp, (coord){0,0,1}, 0.0); fclose (fp);
}
