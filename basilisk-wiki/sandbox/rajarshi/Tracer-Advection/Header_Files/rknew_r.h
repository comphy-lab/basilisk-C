/**
# Runge--Kutta time integrators
*/

static double update (scalar * ul, scalar * kl, double t, double dt,
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
       ulf  = vectors_add(ulf,  f.v);
       klf  = vectors_add(klf,  g.v);
       dulf = vectors_add(dulf, h.v); 
       u1lf = vectors_add(u1lf, i.v);
     }
     else{
       ulc  = list_add(ulc,  f);
       klc  = list_add(klc,  g);
       dulc = list_add(dulc, h);
       u1lc = list_add(u1lc, i); 
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

/**
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
time $t + dt$ using the Runge--Kutta scheme specified by *order*. */


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
       ulf  = vectors_add(ulf,  f.v);
       klf  = vectors_add(klf,  g.v);
       dulf = vectors_add(dulf, h.v); 
     }
     else{
       ulc  = list_add(ulc,  f);
       klc  = list_add(klc,  g);
       dulc = list_add(dulc, h);
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
    w += update (ul, kl, t, dt, Lu, dul, 1.);
    break;
  case 4:
    w += update (ul, kl, t, dt/2., Lu, dul, 2.);
    w += update (ul, kl, t, dt/2., Lu, dul, 2.);
    w += update (ul, kl, t, dt,    Lu, dul, 1.);
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