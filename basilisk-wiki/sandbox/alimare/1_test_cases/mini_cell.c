/**
#Double diffusion of a scalar with a moving embed boundary.

We simulate the diffusion of two tracers separated by an embedded boundary. The
interface moves using the Stefan relation :

$$
  \mathbf{v}_{pc} = \frac{1}{L_H}(\lambda_1 \nabla T_L - \lambda_1 \nabla T_S)
  $$
where $L_H$ is the latent heat of water, $T_L$ and $T_S$ are the temperature
fields on both sides of the interface.

The full algorithm is done on two iterations. It is the following :

  1. Define the diffusion coefficient
  1. Diffuse one tracer
  1. Exchange the inside and the outside of the embed boundary
  1. Define the diffusion coefficient
  1. Diffuse the other tracer
  1. Exchange the inside and the outside of the embed boundary
  1. Move the interface using gradients calculated on both sides of the
  interface
  1. Reinit the level_set function

![Animation of cs*u.x + (1-cs)*u2.x.](mini_cell/visu.mp4)(loop)



~~~gnuplot Phase change velocity
set term pngcairo size 800,800
f(x)  = a  - x/b
f2(x) = a2 - x/b2  
f3(x) = a3 - x/b3  
fit f(x)  'log1' u ($1):(log($7)) via a ,b
fit f2(x) 'log2' u ($1):(log($7)) via a2,b2
fit f3(x) 'log3' u ($1):(log($7)) via a3,b3
ftitle(a,b) = sprintf("%.3fexp^{t/{%4.4f}}",a,b)
latent(a) = sprintf("L_H  %4.4f",a)
set grid
set xlabel "time"
set ylabel "log(v_{pc})"
#set logscale y
aa = 100.
#plot 'log1' u 1:(log($7)) w l  t latent(aa/1), \
#      f(x) w l lw 3 t ftitle(a,b)
plot 'log1' u 1:(log($7)) w l  t latent(aa/1),\
         f(x) w l dt 2 t ftitle(a,b), \
    'log2' u 1:(log($7)) w l pointtype 2  t latent(aa/2), \
        f2(x) dt 2 t ftitle(a2,b2), \
    'log3'  u 1:(log($7)) w l pointtype 2 lc rgb "blue" t latent(aa/3), \
        f3(x) dt 2 t ftitle(a3,b3)
~~~

~~~gnuplot Temperature in the cell located at 
set title 'Temperature in one cell' font 'Helvetica,20'
set key left
plot 'log1' u 1:2 w l t 'Liquid Temperature',  \
     'log1' u 1:4 w l  t 'Solid Temperature',  \
     'log1' u 1:5 w l lw 2 dt 2 lt 8 t 'Approximated Temperature'
~~~


*/

#define DOUBLE_EMBED 1
#define LevelSet     1

#include "embed.h"
#include "advection.h"
#include "diffusion.h"
#include "../level_set.h"
#include "view.h"

#define T_eq         0.
#define TL_inf       1.
#define TS_inf       -1.

int MIN_LEVEL, MAXLEVEL; 
double H0;
double latent_heat;
char filename [100];
FILE * fp1;

#define DT_MAX  1.

#define T_eq         0.
#define plane(x, y, H) (y - H)

scalar TL[], TS[], dist[];
scalar * tracers = {TL, TS};
scalar * level_set = {dist};
face vector v_pc[];
face vector muv[];
mgstats mgT;
scalar grad1[], grad2[];

double lambda[2];

int     nb_cell_NB =  1 << 3 ;  // number of cells for the NB
double  NB_width ;              // length of the NB
  
mgstats mg1,mg2;

TL[embed]  = dirichlet(T_eq);
TL[top]    = dirichlet(TL_inf); 

TS[embed]  = dirichlet(T_eq);
TS[bottom] = dirichlet(TS_inf); 
int j;
int k_loop = 0;
/**
The domain is 4 units long, centered vertically. */

int main() {
  L0 = 4.;
  CFL = 0.5;
  origin (-L0/2., -L0/2.);
  
  j = 1;
  for (j=1;j<=3;j++){

/**
Here we set up the parameters of our simulation. The latent heat $L_H$, the
initial position of the interface $h_0$ and the resolution of the grid.
*/
    latent_heat  = 100./j;
    MAXLEVEL = MIN_LEVEL = 5;

    H0 = -0.3*L0; 
    // H0 = -1.501*L0 / (1 << MAXLEVEL); 
    N = 1 << MAXLEVEL;
    snprintf(filename, 100,  "log%d", j);
    fp1 = fopen (filename,"w");

    init_grid (N);
    run();
    // fclose(fp1); 
  }
}

event init(t=0){


  DT = L0/(1 << MAXLEVEL);

  NB_width = nb_cell_NB * L0 / (1 << MAXLEVEL);

  lambda[0] = 1.;
  lambda[1] = 1.;

  foreach(){
      dist[] = plane(x,y,H0);
  }
  boundary ({dist});
  restriction({dist});

  fractions (dist, cs, fs);
  boundary({cs,fs});
  restriction({cs,fs});

  foreach() {
    TL[] = T_eq;
    TS[] = T_eq;
  }

  foreach_face(){
    v_pc.x[] = 0.;
  }

  boundary({TL,TS});
  restriction({TL,TS});

 
}


event properties(i++){
  

  foreach_face()
    muv.x[] = lambda[i%2]*fs.x[];
  boundary((scalar *) {muv});
}

event tracer_diffusion(i++){
  int kk;
  for (kk=1;kk<=(40*k_loop+1);kk++){
    if(i%2==0){
      mg1 = diffusion(TL, L0/(1 << MAXLEVEL), D = muv);
    }
    else{
      mg2 = diffusion(TS, L0/(1 << MAXLEVEL), D = muv );
    }
  }
}


event LS_advection(i++,last){
  if(i%2 == 1){
    double L_H       = latent_heat;  

    scalar cs0[];


    foreach(){
      cs0[]   = cs[];
      cs[]    = 1.-cs[];
    }
    foreach_face(){
      fs.x[]  = 1.-fs.x[];
    }

    boundary({cs,fs,cs0});
    restriction({cs,fs});

    phase_change_velocity_LS_embed (cs, fs ,TL, TS, v_pc, dist, L_H, 1.05*NB_width,
      nb_cell_NB,lambda);

    advection_LS (level_set, v_pc, dt);
    boundary ({dist});
    restriction({dist});

    fractions (dist, cs, fs);
    boundary({cs,fs});
    restriction({cs,fs});

    foreach(){
      cs[]      = 1.-cs[];
    }
    foreach_face(){
      fs.x[]      = 1.-fs.x[];
    }

    boundary({cs,fs});
    restriction({cs,fs});
    // event ("properties");

    k_loop = 0;
    foreach(){
      if(cs0[] != 1. && cs[] ==1.)k_loop = 1;
    }
  }
}



event LS_reinitialization(i++,last){
  if(i>0 && i%2==1){
    LS_reinit2(dist,L0/(1 << MAXLEVEL), 1.4*NB_width,
      1.4*nb_cell_NB);
  }
}

#if DOUBLE_EMBED
event double_calculation(i++,last){
// modify the fs , cs, copy the outer fields in the partially covered cells

  foreach(){
    cs[]      = 1.-cs[];
  }
  foreach_face(){
    fs.x[]      = 1.-fs.x[];
  }

  boundary({cs,fs});
  restriction({cs,fs});
}
#endif

event movies ( i++,last;t<300.)
{
  if(i%2 == 1 && t > 50) {
    if(i%20==1 && j ==1){
    boundary({TL,TS});
    scalar visu[];
    foreach(){
      visu[] = (cs[])*TL[]+(1.-cs[])*TS[] ;
    }
    boundary({visu});
    
    draw_vof("cs");
    squares("visu", min =-1., max = 1.);
    save ("visu.mp4");
    
    boundary((scalar *){v_pc});
    }

    stats s2    = statsf(v_pc.y);
    Point p     = locate(-1.51*L0/(1<<MIN_LEVEL),-3.51*L0/(1<<MIN_LEVEL));

    double cap  = capteur(p, TL);
    double cap3 = capteur(p, TS);
    double cap2 = capteur(p, cs);
    double T_point = cap2*cap + (1.-cap2)*cap3;
    fprintf (fp1, "%.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",
      t, cap, cap2, cap3, T_point, s2.min, s2.max);
    fprintf(fp1, "## %g %g %d\n", mg1.resa, mg2.resa, k_loop);
  }
}

// event dumps(t=25.){
//   dump();
// }
