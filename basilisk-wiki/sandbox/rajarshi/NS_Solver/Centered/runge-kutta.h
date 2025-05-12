/**
#Header file containing list of all Runge kutta methods.

We wish to solve the following system of equations numerically \
$$ \frac{\partial u}{\partial t} = L(u,t) $$ \
To achieve this we have this header file containing the standard RK4, RK3, RK2 methods
as well as SSP-RK3 method. \


#Standard fourth order Runge-Kutta method

   $$ u^{(n+1)} = u^{(n)} + \frac{dt}{6} (k_1 + 2k_2 + 2k_3 + k_4) $$ \
   $$ k_1 = L(t_{n},u^{(n)}) $$ \
   $$ k_2 = L(t_{n}+\frac{dt}{2},u^{(n)}+\frac{dt}{2}k_1) $$ \
   $$ k_3 = L(t_{n}+\frac{dt}{2},u^{(n)}+\frac{dt}{2}k_2) $$ \
   $$ k_4 = L(t_{n}+dt,u^{(n)}+ dtk_3) $$ \
*/

void RungeKutta4 () {
  
  vector us[],gs[],Lu1[],Lg1[];
  scalar ps[],Lp1[];
  face vector ufs[],Luf1[];

  foreach(){
    ps[] = p[];
    foreach_dimension(){
       us.x[] = u.x[];
       gs.x[] = g.x[];
    }
  }
  foreach_face()
    ufs.x[] = uf.x[];
  boundary({ps});
  boundary((scalar *){us,gs,ufs});

  Update (us,ps,ufs,gs,dt/2.,Lu1,Lg1,Luf1,Lp1); 
   
  vector Lu2[],Lg2[];
  scalar Lp2[];
  face vector Luf2[];

  Update (us,ps,ufs,gs,dt/2.,Lu2,Lg2,Luf2,Lp2);

  foreach(){
    ps[] = p[] + (dt/2.)*Lp2[];
    foreach_dimension(){
       us.x[] = u.x[] + (dt/2.)*Lu2.x[];
       gs.x[] = g.x[] + (dt/2.)*Lg2.x[];
    }
  }
  foreach_face()
    ufs.x[] = uf.x[] + (dt/2.)*Luf2.x[];
  boundary({ps});
  boundary((scalar *){us,gs,ufs});
  
  vector Lu3[],Lg3[];
  scalar Lp3[];
  face vector Luf3[];

  Update (us,ps,ufs,gs,dt/2.,Lu3,Lg3,Luf3,Lp3);

  foreach(){
    ps[] = p[] + dt*Lp3[];
    foreach_dimension(){
       us.x[] = u.x[] + dt*Lu3.x[];
       gs.x[] = g.x[] + dt*Lg3.x[];
    }
  }
  foreach_face()
    ufs.x[] = uf.x[] + dt*Luf3.x[];
  boundary({ps});
  boundary((scalar *){us,gs,ufs});
  
  vector Lu4[],Lg4[];
  scalar Lp4[];
  face vector Luf4[];

  Update (us,ps,ufs,gs,dt/2.,Lu4,Lg4,Luf4,Lp4);

  foreach(){
    p[] += (dt/6.)*(Lp1[] + 2.*Lp2[] + 2.*Lp3[] + Lp4[])/6.;
    foreach_dimension(){
       u.x[] += (dt/6.)*(Lu1.x[] + 2.*Lu2.x[] + 2.*Lu3.x[] + Lu4.x[])/6.;
       g.x[] += (dt/6.)*(Lg1.x[] + 2.*Lg2.x[] + 2.*Lg3.x[] + Lg4.x[])/6.;
    }
  }
  foreach_face()
    uf.x[] += (dt/6.)*(Luf1.x[] + 2.*Luf2.x[] + 2.*Luf3.x[] + Luf4.x[])/6.;
  boundary({p});
  boundary((scalar *){u,g,uf}); 
  
}

/**
#Standard third order Runge-Kutta method

   $$ u^{n+1} = u^{n} + \frac{dt}{6} (k_1 + 4k_2 + k_3) $$ \
   $$ k_1 = L(t_{n},u^{n}) $$ \
   $$ k_2 = L(t_{n}+\frac{dt}{2},u^{n}+\frac{dt}{2}k_1) $$ \
   $$ k_3 = L(t_{n}+dt,u^{n}+dt(-k_1 + 2k_2)) $$ \
 
*/

void RungeKutta3 () {

  vector us[],gs[],Lu1[],Lg1[];
  scalar ps[],Lp1[];
  face vector ufs[],Luf1[];

  foreach(){
    ps[] = p[];
    foreach_dimension(){
       us.x[] = u.x[];
       gs.x[] = g.x[];
    }
  }
  foreach_face()
    ufs.x[] = uf.x[];
  boundary({ps});
  boundary((scalar *){us,gs,ufs});

  Update (us,ps,ufs,gs,dt/2.,Lu1,Lg1,Luf1,Lp1); 

  vector Lu2[],Lg2[];
  scalar Lp2[];
  face vector Luf2[];

  Update (us,ps,ufs,gs,dt/2.,Lu2,Lg2,Luf2,Lp2);

  foreach(){
    ps[] = p[] + dt*(2.*Lp2[] - Lp1[]);
    foreach_dimension(){
       us.x[] = u.x[] + dt*(2.*Lu2.x[] - Lu1.x[]);
       gs.x[] = g.x[] + dt*(2.*Lg2.x[] - Lg1.x[]);
    }
  }
  foreach_face()
    ufs.x[] = uf.x[] + dt*(2.*Luf2.x[] - Luf1.x[]);
  boundary({ps});
  boundary((scalar *){us,gs,ufs});
  
  vector Lu3[],Lg3[];
  scalar Lp3[];
  face vector Luf3[];

  Update (us,ps,ufs,gs,dt/2.,Lu3,Lg3,Luf3,Lp3);

  foreach(){
    p[] += dt*(Lp1[] + 4.*Lp2[] + Lp3[])/6.;
    foreach_dimension(){
       u.x[] += dt*(Lu1.x[] + 4.*Lu2.x[] + Lu3.x[])/6.;
       g.x[] += dt*(Lg1.x[] + 4.*Lg2.x[] + Lg3.x[])/6.;
    }
  }
  foreach_face()
    uf.x[] += dt*(Luf1.x[] + 4.*Luf2.x[] + Luf3.x[])/6.;
  boundary({p});
  boundary((scalar *){u,g,uf});   
  
}

/**
#Standard second order Runge-Kutta method

   $$ u^{n+1} = u^{n} + \frac{dt}{2} (k_1 + k_2) $$ \
   $$ k_1 = L(t_{n},u^{n}) $$ \
   $$ k_2 = L(t_{n}+dt,u^{n}+dtk_1) $$ \
 
*/

void RungeKutta2 (){
  
  vector us[],gs[],Lu1[],Lg1[];
  scalar ps[],Lp1[];
  face vector ufs[],Luf1[];

  foreach(){
    ps[] = p[];
    foreach_dimension(){
       us.x[] = u.x[];
       gs.x[] = g.x[];
    }
  }
  foreach_face()
    ufs.x[] = uf.x[];
  boundary({ps});
  boundary((scalar *){us,gs,ufs});
  
  Update (us,ps,ufs,gs,dt,Lu1,Lg1,Luf1,Lp1);
  
  vector Lu2[],Lg2[];
  scalar Lp2[];
  face vector Luf2[];
  
  Update (us,ps,ufs,gs,dt,Lu2,Lg2,Luf2,Lp2);
  
  foreach(){
    p[] += dt*(Lp1[] + Lp2[])/2.;
    foreach_dimension(){
       u.x[] += dt*(Lu1.x[] + Lu2.x[])/2.;
       g.x[] += dt*(Lg1.x[] + Lg2.x[])/2.;
    }
  }
  foreach_face()
    uf.x[] = dt*(Luf1.x[] + Luf2.x[])/2.;
 
  boundary({p});
  boundary((scalar *){u,uf,g});
}

/**
#Strong stability preserving third order Runge-Kutta method

   $$ u^{(1)} = u^{(n)} + dtL(u^{(n)}) $$ \
   $$ u^{(2)} = \frac{3}{4}u^{(n)} + \frac{1}{4}u^{(1)} + \frac{1}{4}L(u^{(1)}) $$ \
   $$ u^{(n+1)} = \frac{1}{3}u^{(n)} + \frac{2}{3}u^{(2)} + \frac{2}{3}L(u^{(2)}) $$ \
 
*/


void SSP_RungeKutta3 (){
  
  vector u1[],g1[],Lu[],Lg[];
  scalar p1[],Lp[];
  face vector uf1[],Luf[];

  foreach(){
    p1[] = p[];
    foreach_dimension(){
       u1.x[] = u.x[];
       g1.x[] = g.x[];
    }
  }
  foreach_face()
    uf1.x[] = uf.x[];
  boundary({p1});
  boundary((scalar *){u1,g1,uf1});
  
  Update (u1,p1,uf1,g1,dt,Lu,Lg,Luf,Lp);
  
  vector u2[],g2[];
  scalar p2[];
  face vector uf2[];

  foreach(){
    p2[] = p1[];
    foreach_dimension(){
       u2.x[] = u1.x[];
       g2.x[] = g1.x[];
    }
  }
  foreach_face()
    uf2.x[] = uf1.x[];
  boundary({p2});
  boundary((scalar *){u2,g2,uf2});
  
  Update (u2,p2,uf2,g2,dt/4.,Lu,Lg,Luf,Lp);
  
  foreach(){
    p2[] += (3./4.)*(p[] - p1[]);
    foreach_dimension(){
       u2.x[] += (3./4.)*(u.x[] - u1.x[]);
       g2.x[] += (3./4.)*(g.x[] - g1.x[]);
    }
  }
  foreach_face()
    uf2.x[] += (3./4.)*(uf.x[] - uf1.x[])/2.;
   
  boundary({p2});
  boundary((scalar *){u2,uf2,g2});

  foreach(){
    p[] = (p[]-p2[])/3.;
    foreach_dimension(){
       u.x[] = (u.x[]-u2.x[])/3.;
       g.x[] = (g.x[]-g2.x[])/3.;
    }
  }
  foreach_face()
    uf.x[] = (uf.x[]-uf2.x[])/3.;

  Update (u2,p2,uf2,g2,2.*dt/3.,Lu,Lg,Luf,Lp);

  foreach(){
    p[] += p2[];
    foreach_dimension(){
       u.x[] += u2.x[];
       g.x[] += g2.x[];
    }
  }
  foreach_face()
    uf.x[] += uf2.x[];

  boundary({p});
  boundary((scalar *){u,g,uf});
  
}