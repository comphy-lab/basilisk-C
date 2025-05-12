#define BGHOSTS 2
#define growthTimestep 0.1
#include "run.h"
#include "timestep-vdw.h"
#include "difference.h"

scalar normGradRho[];
scalar rho[],rhop[];
vector q[],u[];
scalar mu[];
double P0( double x);



(const) double lambda;
(const) double miu;

event defaults (i = 0)
{
  CFL = 0.1;
}

double dtmax;

event init (i = 0)
{
  boundary (all);
  event ("properties");
  dtmax = DT;
  event ("stability");
}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  double mydtmax = HUGE;
  foreach(reduction(min:mydtmax)){
    if ( miu != 0.) {

      double dt = 0.1*sq(Delta)/miu;

      if (dt < mydtmax) mydtmax = dt;
    }
  }
  dtmax = min(mydtmax,dtmax);
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  dt = dtnext (timestep (u, q, dtmax));
}
    scalar pp[];
    scalar uu[];
	  scalar vv[];

    vector adv[];
    vector surf1[];
    vector surf2[];
    vector vis[];


event predictor_corrector(i++,last)
{
  scalar laprho[],rhop[];
  vector qp[];

  trash({rhop});
    foreach()
     {

      rhop[] = rho[] - dt*(backward_x(point,q.x)+backward_y(point,q.y));
      foreach_dimension(){
	      u.x[] = q.x[]/rho[];}
     }
    boundary({rhop});
    boundary(( scalar *) {u});
    
    foreach(){
      vv[]=mu[]*(center_y(point,u.x)+center_x(point,u.y));
      pp[]=P0(rho[])-lambda*rho[]*laplace(point,rho);    
      uu[]=rho[]*u.x[]*u.y[];  

      vis.x[]=mu[]*(4.0/3.0*center_x(point,u.x)-2.0/3.0*center_y(point,u.y));
      vis.y[]=mu[]*(4.0/3.0*center_y(point,u.y)-2.0/3.0*center_x(point,u.x));

      adv.x[]=rho[]*sq(u.x[]);
      adv.y[]=rho[]*sq(u.y[]);

      surf1.x[]=0.5*lambda*(sq(center_y(point,rho))-sq(forward_x(point,rho)));
      surf1.y[]=0.5*lambda*(sq(center_x(point,rho))-sq(forward_y(point,rho)));
      surf2.x[]=-lambda*forward_y(point,rho)*center_x(point,rho);
      surf2.y[]=-lambda*forward_x(point,rho)*center_y(point,rho);

    }
    trash({qp});
    foreach(){
        qp.x[] = q.x[] - dt*         (backward_x(point,adv.x) +     backward_y(point,uu) + backward_x(point,pp)       
                       -             (backward_x(point,surf1.x ) +     backward_y(point,surf2.x)) 
                       -             (center_x(point,vis.x)  +     center_y(point,vv)));

        qp.y[] = q.y[] - dt*         (backward_y(point,adv.y) +     backward_x(point,uu) +backward_y(point,pp)
                       -             (backward_y(point,surf1.y ) +     backward_x(point,surf2.y))
                       -             (center_y(point,vis.y)  +     center_x(point,vv)));
    }                                    

  boundary(( scalar *) {qp});
  foreach()
  {

    rho[] = 0.5*(rhop[]+rho[]) - 0.5*dt*(forward_x(point,qp.x) + forward_y(point,qp.y));
    foreach_dimension()
	    u.x[] = qp.x[]/rhop[];
  }
  boundary({rho});
  boundary(( scalar *) {u});
  foreach(){

      uu[]=rhop[]*u.x[]*u.y[];
      vv[]=mu[]*(center_y(point,u.x)+center_x(point,u.y));
      pp[]=P0(rhop[])-lambda*rhop[]*laplace(point,rhop);

      vis.x[]=mu[]*(4.0/3.0*center_x(point,u.x)-2.0/3.0*center_y(point,u.y));
      vis.y[]=mu[]*(4.0/3.0*center_y(point,u.y)-2.0/3.0*center_x(point,u.x));
      
      adv.x[]=rhop[]*sq(u.x[]);
      adv.y[]=rhop[]*sq(u.y[]);

      surf1.x[]=0.5*lambda*(sq(center_y(point,rhop))-sq(backward_x(point,rhop)));
      surf1.y[]=0.5*lambda*(sq(center_x(point,rhop))-sq(backward_y(point,rhop)));
      surf2.x[]=-lambda*backward_y(point,rhop)*center_x(point,rhop);
      surf2.y[]=-lambda*backward_x(point,rhop)*center_y(point,rhop);

  }

  foreach(){
  q.x[] = 0.5*(qp.x[] + q.x[]) - 0.5*dt*    (forward_x(point,adv.x) +    forward_y(point,uu) +forward_x(point,pp)
                               -            (forward_x(point,surf1.x ) +    forward_y(point,surf2.x))
                               -            (center_x(point,vis.x) +    center_y(point,vv)));

  q.y[] = 0.5*(qp.y[] + q.y[]) - 0.5*dt*    (forward_y(point,adv.y) +    forward_x(point,uu) +forward_y(point,pp)
                               -            (forward_y(point,surf1.y ) +    forward_x(point,surf2.y))
                               -            (center_y(point,vis.y) +    center_x(point,vv)));
}
  boundary(( scalar *) {q});
  restriction((scalar *){q});

}
	
