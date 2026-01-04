#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

#define r0 6e-5
double tmax=5e-5;
double dtuser=1e-6;


// u.n[right]  = dirichlet(0.0);
// p[right]    = neumann(0.);
// pf[right]   = neumann(0.);

u.n[left]  = neumann(0.0);
p[left]    = dirichlet(0.);
pf[left]   = dirichlet(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[top] = neumann(0.);
p[top]   = dirichlet(0.);
pf[top]  = dirichlet(0.);


u.n[bottom] = neumann(0.);
p[bottom]   = dirichlet(0.);
pf[bottom]  = dirichlet(0.);


u.n[embed] = dirichlet( x > 0 ? 1e-6 : -1e-6 );
u.t[embed] = dirichlet( y > 0 ? -1e-6 : 1e-6 );
//u.n[embed] = dirichlet(x/r0*1e-6);


face vector muv[];
scalar rhog[];


int main(void){
  L0 = 15.*r0; //m
  origin(-L0/2.,-L0/2.);
  DT=dtuser;
  mu = muv;
  rho = rhog;
  run();
}

event init(i=0){
  init_grid(64);
  solid(cs,fs,sq(x)+sq(y)-sq(r0)); 
  fractions_cleanup(cs,fs,1e-4);
}

event properties(i++){
  foreach(){
    rhog[] = (cs[] > 0 ? 1. : nodata);
  }
  foreach_face(){
    muv.x[]=0.2*fm.x[];
  }
}


event view(t=0; t<tmax; t+=2*dt){
  clear();
  scalar U[];
  foreach(){
    U[] = sqrt(sq(u.x[])+sq(u.y[]));
  }
  view (fov = 25, tx = -0.2, width = 1100, height = 800);
  squares ("U", spread = -1,
     cbar = true, border = true, pos = {0.5,  0.0},
     label = "U", mid = true, format = "%9.5f", levels = 10);
  squares("U", linear=false, map=jet);
  draw_vof ("cs",filled=0, fc = {0.,0.8,0.8}, lw = 2);
  box();
  vectors("u",scale = 0.1);
  save ("u.mp4");

}

event end(t=tmax);

/**
![Results](embed_stefanflow/u.mp4)(width="80%")
*/