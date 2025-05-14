#include "../output_fields/output_vtu_foreach.h"

event snapshots2d (t += tsample3) {
  printf( "Exporting slices to VTU format \n" );

  static int nsnap2d = 0;
  char name[80];
  sprintf(name, "snapshots2d_x0_%3.3d", nsnap2d);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){1,0,0}, 0.0);

  sprintf(name, "snapshots2d_y0_%3.3d", nsnap2d);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){0,1,0}, 0.0);

  sprintf(name, "snapshots2d_z0_%3.3d", nsnap2d);
  output_vtu_plane ((scalar *) {l2}, (vector *) {u, omega}, name, (coord){0,0,1}, 0.0);

  nsnap2d++;
}

event snapshots3d (t += tsample4) {
  printf( "Exporting snapshots to VTU format \n" );

  static int nsnap3d = 0;
  char name[80];
  sprintf(name, "snapshots3d_%3.3d", nsnap3d);
  output_vtu ((scalar *) {l2}, (vector *) {u, omega}, name);
  nsnap3d++;

  sprintf (name, "dump-%g", t);
  dump (name);
}

event movie (t += tsample1) {

  view(camera="iso", fov=4*L0);
  isosurface ("l2", 0);
  squares ("u.z", linear = false, alpha=-L0/2);

  box();
  save ("lambda2.mp4");
}
