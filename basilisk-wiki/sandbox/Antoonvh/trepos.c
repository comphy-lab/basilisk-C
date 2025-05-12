/**
# Test for the reposition function

![Before repositioning](trepos/before.png)

Not all particles are repositioned:

![After repositioning](trepos/after.png)

![After repositioning 2](trepos/after2.png)

![After repositioning 3](trepos/after3.png)

*/
#include "vof-tracer-particles.h"
#include "view.h"
#include "scatter2.h"

vector uf[];
scalar f[];

int main() {
  L0 = 3;
  X0 = Y0 = -L0/2;
  init_grid (8);
  int pnr = 100;
  Particles P = new_vof_tracer_particles (sq(pnr), 1);
  place_in_square (P, (struct Init_P){pnr});
  fraction(f, sq(1) - sq(x-0.0432) - sq(y + 0.1423));
  boundary ({f});
  scatter (P, s = 10);
  draw_vof ("f", lc = {1,0,1}, lw = 3);
  save ("before.png");

  reposition (f, P); 
  cells(lc = {0,1,1});
  scatter (P, s = 10);
  draw_vof ("f", lc = {1,0,1}, lw = 3);
  save ("after.png");

  Particles P2 = new_vof_tracer_particles (sq(pnr), 0);
  place_in_square (P2, (struct Init_P){pnr});
  reposition (f, P2); 
  cells (lc = {0,1,1});
  scatter (P2, s = 10);
  draw_vof ("f", lc = {1,0,1}, lw = 3);
  save ("after2.png");
  
  Particles P3 = new_vof_tracer_particles (sq(pnr), 3);
  place_in_square (P3, (struct Init_P){pnr});
  reposition (f, P3); 
  cells (lc = {0,1,1});
  scatter (P3, s = 10);
  draw_vof ("f", lc = {1,0,1}, lw = 3);
  save ("after3.png");
  
  run();
}
