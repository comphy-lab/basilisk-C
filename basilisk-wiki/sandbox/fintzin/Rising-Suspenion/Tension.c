#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "output.h"

#include "parameters.h"

#include "algebra.h"
#include "MISC.h"
FILE * fDATA;

event my_debug(i++){
  if(i==0){
    fDATA = fopen("DATA.csv","w");
    fprintf(fDATA,"t,A_x,A_y,G_x_x,G_y_y,G_x_x_x,G_y_y_y,P_x,P_y,M_x,M_y,As_x,As_y,Ms_x,Ms_y,MsV_x,MsV_y,vol,D,n\n");
  }
  double vol = 0;
  coord pos = {0,0};
  coord asym = {0,0};
  coord G = {0,0};
  coord A={0,0};
  coord As={0,0};
  coord M = {0};
  coord Ms = {0};
  coord MsV = {0};
  foreach(){
    coord poss= POS;
    vol += f[] * dv();
    foreach_dimension()
      pos.x += poss.x * f[] * dv();
  }

  pos.x /= vol;
  pos.y /= vol;
  double rad =sqrt(vol/M_PI);
  int nc = 1<<(LEVEL);
  double cell_per_diam = rad*2 / Ls * nc;

  foreach(){
    coord X = {x,y};
    foreach_dimension(){
      G.x    += sq(X.x - pos.x) * f[] * dv();
      asym.x += cube(X.x - pos.x) * f[] * dv();
    }
  }

  foreach_face(){
    coord X = {x,y};
    A.x += a.x[]*dv();
    M.x += (X.x - pos.x) * a.x[]*((f[]+f[-1])/2*(rho1-rho2)+rho2) * dv();
  }

  foreach_face(){
    coord X = {x,y};
    if (f[] != f[-1] && fm.x[] > 0.){
      As.x += a.x[]*dv();
      Ms.x += (X.x - pos.x) * a.x[]*((f[]+f[-1])/2*(rho1-rho2)+rho2) * dv();
    }
  }

  foreach(){
    coord X = {x,y};
    foreach_dimension()
    if (f[1] != f[-1] && fm.x[] > 0.){
      MsV.x += (X.x - pos.x) * (a.x[1]+a.x[])/2 * (f[]*(rho1-rho2)+rho2) * dv();
    }
  }

  fprintf(fDATA,"%g,%g,%g,%g,%g,%g,%g",t,A.x,A.y,G.x,G.y,asym.x,asym.y);
  fprintf(fDATA,",%g,%g",pos.x,pos.y);
  fprintf(fDATA,",%g,%g",M.x,M.y);
  fprintf(fDATA,",%g,%g",As.x,As.y);
  fprintf(fDATA,",%g,%g",Ms.x,Ms.y);
  fprintf(fDATA,",%g,%g",MsV.x,MsV.y);
  fprintf(fDATA,",%g,%g,%g\n",vol,rad,cell_per_diam);
}


event init (t=0){
  double Rm = 0.5;
  double a = 0.7;
  double b = sq(Rm)/a;
  double k = 1.5;
  foreach(){
    if(!strcmp(Shape, "Sph"))
      f[] = sq(x) + sq(y) < sq(D/2); //sphere 
    if(!strcmp(Shape, "Elli"))
      f[] = sq(x/a) + sq(y/b) < 1.; // ellipse
    if(!strcmp(Shape, "Asym"))
      f[] = sq(x) + sq(y) * (1+k*x) < sq(D/2); // ovoid 
    u.x[]=u.y[]=u.z[]=0.;
  }
}


event movies(t=0; t<=1; t+=0.1){
  box();
  draw_vof("f", lw = 2.);
  squares("u.y");
  save("uy.mp4");
  clear();

  box();
  draw_vof("f", lw = 2.);
  squares("p");
  save("P.mp4");
  clear();
  
  box();
  squares("a.y");
  save("a.mp4");
  clear();

  box();
  squares("f");
  save("f.mp4");
  clear();
}

event stop(t=TMAX){
  return 1;
}

/** The flow parameters are calculated based on the previously defined non-dimensional parameters. 
 * The tolerance is reduced to 1e-4 and the boundaries are set to be periodic*/
int main(){
  size(Ls);
  init_grid(1<<(LEVEL));
  origin(-Ls/2.,-Ls/2.,-Ls/2.);
  rho1=r1;
  mu1=mu_d;//sqrt(fabs(1-rho_r)*rho_r*g*(D*D*D))*rho1/(Ga*mu_r);
  rho2=rho_f;//rho1*rho_r;
  mu2=mu_f;//mu1*mu_r;
  f.sigma=sig;//fabs(1-rho_r)*rho1*g*D*D/Bo;
  TOLERANCE=1e-4;
  foreach_dimension()
    periodic(right);
  run();
}