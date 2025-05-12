#include "poisson.h"


//set up magnetic fluid parameters
#define MU0 1.25663706e-6
double chi1 = 0.0, chi2 = 0.0;


//used scalars
scalar phi_mag[];
vector H[];

//default boundary conditions on phi. approximately set H field normal to boundary to be symmetric (neuman(0))
phi_mag[left] = 2*phi_mag[] - phi_mag[1,0,0];
phi_mag[right] = 2*phi_mag[] - phi_mag[-1,0,0];
phi_mag[top] = 2*phi_mag[] - phi_mag[0,-1,0];
phi_mag[bottom] = 2*phi_mag[] - phi_mag[0,1,0];
phi_mag[front] = 2*phi_mag[] - phi_mag[0,0,-1];
phi_mag[back] = 2*phi_mag[] - phi_mag[0,0,1];

//H default boundary conditions are default Basilisk boundary conditons


//harmonic averaging of magnetic permeability
double compute_mu(double f){
  return 1.0/(f/(MU0*(chi1+1)) + (1-f)/(MU0*(chi2+1)));
}


//arithmetic averaging of magnetic permeability
//double compute_mu(double f){
//  return f*MU0*(chi1+1) + (1-f)*(MU0*(chi2+1));
//}


//compute phi field with Poisson solver
void compute_phi(){
  boundary({f});
  face vector mu[];
  foreach_face()
    mu.x[] = compute_mu((f[]+f[-1])/2);
  scalar s[];
  foreach()
   s[] = 0.0;

  mgstats q = poisson(a = phi_mag, b = s, alpha = mu);
  fprintf(fout, "%d %d\n", q.i, q.nrelax);

}

face vector Fmm[];


//formulations for magnetic body force density

//mu H formulation
/*
event acceleration(i++){
  boundary({f});
  compute_phi();
  face vector av = a;

  foreach_dimension()
    foreach()
      H.x[] = -(phi_mag[1]-phi_mag[-1])/2/Delta;
  boundary((scalar *){H});

  foreach_face(){
    double Hx = face_value(H.x,0);
    double Hy = face_value(H.y,0);
    double Hz = face_value(H.z,0);
    double dHxdx = face_gradient_x(H.x, 0);
    double dHydx = face_gradient_x(H.y, 0);
    double dHzdx = face_gradient_x(H.z, 0);


    Fmm.x[] = compute_mu((f[]+f[-1])/2)*(Hx*dHxdx + Hy*dHydx + Hz*dHzdx);
    av.x[] += alpha.x[]/fm.x[]*compute_mu((f[]+f[-1])/2)*(Hx*dHxdx + Hy*dHydx + Hz*dHzdx);

  }
}*/


//H^2 grad mu formulation


event acceleration(i++){
  boundary({f});
  compute_phi();
  face vector av = a;

  foreach_dimension()
    foreach()
      H.x[] = -(phi_mag[1]-phi_mag[-1])/2/Delta;
  boundary((scalar *){H});

  foreach_face(){
    double Hx = face_value(H.x,0);
    double Hy = face_value(H.y,0);
    double Hz = face_value(H.z,0);
    Fmm.x[] = (compute_mu(f[])-compute_mu(f[-1]))/Delta*(Hx*Hx + Hy*Hy + Hz*Hz);
    av.x[] -= alpha.x[]/fm.x[]*(compute_mu(f[])-compute_mu(f[-1]))/Delta*(Hx*Hx + Hy*Hy + Hz*Hz);
  }
}
