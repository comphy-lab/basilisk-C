/**
# Runge--Kutta time integrators

The *runge_kutta()* function implements the classical first- (Euler),
second- and fourth-order Runge--Kutta time integrators for evolution
equations of the form
$$
\frac{\partial\mathbf{u}}{\partial t} = L(\mathbf{u}, t)
$$
with $\mathbf{u}$ a vector (i.e. list) of evolving fields and $L()$ a
generic, user-defined operator.

Given $\mathbf{u}$, the initial time *t*, a timestep *dt* and the
function $L()$ which should fill *kl* with the right-hand-side of the
evolution equation, the function below will return $\mathbf{u}$ at
time $t + dt$ using the Runge--Kutta scheme specified by *order*.

We also have an implementation for a Strong-stability-preserving Runge-Kutta-3
scheme in this header-file.
*/


static double updaterk (scalar * ul, scalar * kl, double t, double dt,
		      void (* Lu) (scalar * ul, double t, scalar * kl),
		      scalar * dul, double w)
{

  scalar * u1l = list_clone (ul);
  vector * ulf = NULL;
  scalar * ulc = NULL;
  vector * u1lf = NULL;
  scalar * u1lc = NULL;
  vector * klf = NULL;
  scalar * klc = NULL;
  vector * dulf = NULL;
  scalar * dulc = NULL;
   
  scalar f,g,h,i;
  for(f,g,h,i in ul,kl,dul,u1l){
     if(f.face){
       ulf  = vectors_append(ulf,  f.v);
       klf  = vectors_append(klf,  g.v);
       dulf = vectors_append(dulf, h.v); 
       u1lf = vectors_append(u1lf, i.v);
     }
     else{
       ulc  = list_append(ulc,  f);
       klc  = list_append(klc,  g);
       dulc = list_append(dulc, h);
       u1lc = list_append(u1lc, i); 
     }
  }

  if(ulc!=NULL){
    foreach(){
      scalar u1,u,k;
      for(u1,u,k in u1lc,ulc,klc)
         u1[] = u[] + dt*k[];
    }
  }
  if(ulf!=NULL){
    foreach_face(){
      vector u1,u,k;
      for(u1,u,k in u1lf,ulf,klf)
         u1.x[] = u.x[] + dt*k.x[];
    }
  }
  boundary (u1l);
  
  Lu (u1l, t + dt, kl);
 

  if(ulc!=NULL){
    foreach(){
      scalar du,k;
      for(du,k in dulc,klc)
         du[] += w*k[];
    }
  }
  if(ulf!=NULL){
    foreach_face(){
      vector du,k;
      for(du,k in dulf,klf)
         du.x[] += w*k.x[];
    }
  }
  boundary(dul);
 
  free(ulf);
  free(ulc);
  free(u1lf);
  free(u1lc);
  free(klf);
  free(klc);
  free(dulf);
  free(dulc); 

  delete (u1l), free (u1l);
  return w;
}

void runge_kutta (scalar * ul, double t, double dt,
		  void (* Lu) (scalar * ul, double t, scalar * kl),
		  int order)
{
  scalar * dul = list_clone (ul);
  scalar * kl = list_clone (ul);
  vector * ulf  = NULL;
  scalar * ulc  = NULL;
  vector * dulf = NULL;
  scalar * dulc = NULL;
  vector * klf  = NULL;
  scalar * klc  = NULL;

  scalar f,g,h;
  for(f,g,h in ul,kl,dul){
     if(f.face){
       ulf  = vectors_append(ulf,  f.v);
       klf  = vectors_append(klf,  g.v);
       dulf = vectors_append(dulf, h.v); 
     }
     else{
       ulc  = list_append (ulc,  f);
       klc  = list_append (klc,  g);
       dulc = list_append (dulc, h);
     }
  }

  Lu (ul, t, kl);
  if(ulc!=NULL){
    foreach(){
      scalar du,k;
      for(du,k in dulc,klc)
         du[] = k[];
    }
  }
  if(ulf!=NULL){
    foreach_face(){
      vector du,k;
      for(du,k in dulf,klf)
         du.x[] = k.x[];
    }
  }
  boundary(dul);

  double w = 1.;
  switch (order) {
  case 1: // Euler
    break;
  case 2:
    w += updaterk (ul, kl, t, dt, Lu, dul, 1.);
    break;
  case 4:
    w += updaterk (ul, kl, t, dt/2., Lu, dul, 2.);
    w += updaterk (ul, kl, t, dt/2., Lu, dul, 2.);
    w += updaterk (ul, kl, t, dt,    Lu, dul, 1.);
    break;
  default:
    assert (false); // not implemented
  }
 
  if(ulc!=NULL){ 
    foreach(){
      scalar u,du;
      for(u,du in ulc,dulc)
         u[] += dt/w*du[];
    }
  }
  if(ulf!=NULL){
    foreach_face(){
      vector u,du;
      for(u,du in ulf,dulf)
         u.x[] += dt/w*du.x[];
    }
  }
  boundary(ul);

  free(ulf);
  free(ulc);
  free(klf);
  free(klc);
  free(dulf);
  free(dulc);
  delete (dul), free (dul);
  delete (kl), free (kl);
}

void SSP_Runge_Kutta3 (scalar * ul, double t, double dt,
		  void (* Lu) (scalar * ul, double t, scalar * kl))
{

  scalar * u1l = list_clone(ul);
  scalar * kln = list_clone(ul);   // L(un)
  scalar * kl1 = list_clone(ul);   // L(u1)
  scalar * kl2 = list_clone(ul);   // L(u2)
  
  vector * ulf   = NULL;
  scalar * ulc   = NULL;
  vector * u1lf  = NULL;
  scalar * u1lc  = NULL;
  vector * klnf  = NULL;
  scalar * klnc  = NULL;
  vector * kl1f  = NULL;
  scalar * kl1c  = NULL;
  vector * kl2f  = NULL;
  scalar * kl2c  = NULL;

  scalar f,g,h,i,j;
  for(f,g,h,i,j in ul,u1l,kln,kl1,kl2){
     if(f.face){
       ulf  = vectors_add (ulf , f.v);
       u1lf = vectors_add (u1lf, g.v);
       klnf = vectors_add (klnf, h.v);
       kl1f = vectors_add (kl1f, i.v);
       kl2f = vectors_add (kl2f, j.v);  
     }
     else{
       ulc   = list_append (ulc , f);
       u1lc  = list_append (u1lc, g);
       klnc  = list_append (klnc, h);
       kl1c  = list_append (kl1c, i);
       kl2c  = list_append (kl2c, j); 
     }
  }

  Lu (ul, t, kln);

  if(ulc!=NULL){
     foreach(){
         scalar f,g,h; 
         for (f,g,h in u1lc,ulc,klnc)
              f[] = g[] + dt*h[];
     }
  }
  if(ulf!=NULL){
     foreach_face(){
         vector fv,gv,hv; 
         for (fv,gv,hv in u1lf,ulf,klnf)
              fv.x[] = gv.x[] + dt*hv.x[];
     }
  }
  boundary(u1l);
  Lu (u1l, t+dt, kl1);

  if(ulc!=NULL){
     foreach(){
         scalar f,g,h,i;
         for (f,g,h,i in u1lc,ulc,kl1c,klnc)
             f[] = g[] + dt*(h[]+i[])/4.;
     }
  }
  if(ulf!=NULL){
     foreach_face(){
         vector fv,gv,hv,iv; 
         for (fv,gv,hv,iv in u1lf,ulf,kl1f,klnf)
            fv.x[] = gv.x[] + dt*(hv.x[]+iv.x[])/4.;
     }  
  }
  boundary(u1l);
  Lu (u1l, t+dt/2., kl2);
  
  if(ulc!=NULL){
      foreach(){
          scalar f,g,h,i;
          for (f,g,h,i in ulc,klnc,kl1c,kl2c)
              f[] += dt*( g[] + h[] + 4.*i[] )/6.;
      }
  } 
  if(ulf!=NULL){
      foreach_face(){
          vector fv,gv,hv,iv; 
          for (fv,gv,hv,iv in ulf,klnf,kl1f,kl2f)
             fv.x[] += dt*( gv.x[] + hv.x[] + 4.*iv.x[] )/6.;
      }
  }
  boundary(ul);

  free(ulf);
  free(ulc);
  free(u1lf);
  free(u1lc);
  free(klnf);
  free(klnc);
  free(kl1f);
  free(kl1c);
  free(kl2f);
  free(kl2c);
  
  delete (kln), free (kln);
  delete (kl1), free (kl1);
  delete (kl2), free (kl2);
  delete (u1l), free (u1l); 
           
}


void runge_kutta_D (scalar * ul, double t, double dt,
		  void (* Lu) (scalar * ul, double t, scalar * kl),
		  int order)
{
   if(order==2 || order==4)
        runge_kutta (tracers, t , dt, tracer_fluxes, order);
   else
        SSP_Runge_Kutta3 (tracers, t , dt, tracer_fluxes);
   

}
