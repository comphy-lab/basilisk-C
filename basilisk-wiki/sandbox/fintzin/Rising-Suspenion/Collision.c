#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "view.h"
#include "output.h"

#include "parameters.h"

#include "algebra.h"
scalar f1[], f2[], * interfaces = {f1,f2};
/**Function that return true if the current process is of rank 0 */
int main_process(){
  #if _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  return rank == 0;
  #else
  return 1;
  #endif
}
FILE * fDATA;
event computation(t=0;t+=0.05){
  if(i==0){
    fDATA = fopen("DATA.csv","w");
    if(main_process()){
      fprintf(fDATA,"t,A1_x,A1_y,G1_xx,G1_yy,G1_xxx,G1_yyy,P1_x,P1_y,M1_xx,M1_yy,U1_x,U1_y,vol1,D1,n1,");
      fprintf(fDATA,"A2_x,A2_y,G2_xx,G2_yy,G2_xxx,G2_yyy,P2_x,P2_y,M2_xx,M2_yy,U2_x,U2_y,vol2,D2,n2,");
      fprintf(fDATA,"A_x,A_y,Rhou_x,Rhou_y\n");
    }
  }else{
    fDATA = fopen("DATA.csv","a");
  }
  if(main_process())
    fprintf(fDATA,"%g",t);

  for(scalar f in interfaces){

    double vol = 0;
    coord pos = {0,0};
    coord asym = {0,0};
    coord G = {0,0};
    coord As={0,0};
    coord U={0,0};
    coord Ms = {0,0};
    foreach(reduction(+:vol) reduction(+:pos) reduction(+:U)){
      coord poss = {x,y};
      vol += f[] * dv();
      foreach_dimension(){
        pos.x += poss.x * f[] * dv();
        U.x += u.x[] * f[] *dv();
      }
    }

    double rad,cell_per_diam=0.;
    if(vol != 0.){
      foreach_dimension(){
        pos.x /= vol;
        U.x /= vol;
      }
      rad =sqrt(vol/M_PI);
      cell_per_diam = rad*2 / Ls * (1<<(LEVEL));
    }

    foreach(reduction(+:G) reduction(+:asym)){
      coord X = {x,y};
      foreach_dimension(){
        G.x    += sq(X.x - pos.x) * f[] * dv();
        asym.x += cube(X.x - pos.x) * f[] * dv();
      }
    }

    foreach_face(reduction(+:As) reduction(+:Ms)){
      coord X = {x,y};
      if (f[] != f[-1] && fm.x[] > 0.){
        As.x += a.x[]*dv();
        Ms.x += (X.x - pos.x) * a.x[]*dv();
      }
    }


    
    if(main_process()){
      fprintf(fDATA,",%g,%g,%g,%g,%g,%g",As.x,As.y,G.x,G.y,asym.x,asym.y);
      fprintf(fDATA,",%g,%g",pos.x,pos.y);
      fprintf(fDATA,",%g,%g",Ms.x,Ms.y);
      fprintf(fDATA,",%g,%g",U.x,U.y);
      fprintf(fDATA,",%g,%g,%g",vol,rad,cell_per_diam);
    }
  }
  coord RhoU={0,0};
  coord A={0,0};
  foreach(reduction(+:RhoU))
    foreach_dimension()
      RhoU.x += u.x[];//average momentum
  foreach_face(reduction(+:A))
    A.x += a.x[]*dv();
  

  if(main_process()){
    fprintf(fDATA,",%g,%g",A.x,A.y);
    fprintf(fDATA,",%g,%g",RhoU.x,RhoU.y);
    fprintf(fDATA,"\n");
  }
  fclose(fDATA);
}


event init (t=0){
  fraction(f1, - (sq(x + 1.) + sq(y) - sq(D/2)));
  fraction(f2, - (sq(x - 1.) + sq(y) - sq(D/2)));
  foreach()
    u.x[] = f1[] - f2[];
}


event movies(t=0; t<=TMAX; t+=0.1){
  box();
  draw_vof("f1", lw = 2.);
  draw_vof("f2", lw = 2.);
  squares("u.x");
  save("ux.mp4");
  clear();

  box();
  draw_vof("f1", lw = 2.);
  draw_vof("f2", lw = 2.);
  squares("p");
  save("P.mp4");
  clear();
  
  box();
  squares("a.x");
  save("a.mp4");
  clear();
}

event stop(t=TMAX){
  return 1;
}

int main(){
  size(Ls);
  init_grid(1<<(LEVEL));
  origin(-Ls/2.,-Ls/2.,-Ls/2.);
  // rho1=r1;
  // mu1=mu_d;//sqrt(fabs(1-rho_r)*rho_r*g*(D*D*D))*rho1/(Ga*mu_r);
  // rho2=rho_f;//rho1*rho_r;
  // mu2=mu_f;//mu1*mu_r;
  const face vector muc[] = {mu_f,mu_f};
  mu = muc;
  f1.sigma=f2.sigma=sig;//fabs(1-rho_r)*rho1*g*D*D/Bo;
  TOLERANCE=1e-4;
  foreach_dimension()
    periodic(right);
  run();
}