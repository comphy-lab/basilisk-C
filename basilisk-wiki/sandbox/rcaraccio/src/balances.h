/**
# Mass balances
This file contains functions to compute and log mass balances in two-phase flow simulations with phase change in porous media.
*/

#ifndef BALANCES
# define BALANCES
#endif

extern int maxlevel;
extern scalar f, porosity;

struct MassBalances 
{
  //total mass variables
  double tot_sol_mass_start;  // Total initial solid mass
  double tot_gas_mass_start;  // Total initial gas mass
  double tot_sol_mass;        // Total current solid mass
  double tot_gas_mass;        // Total current gas mass
  double tot_mass_start;      // Total initial mass 
  double tot_mass;            // Total current mass
  double tot_gas_mass_bdnow;  // Total gas mass that has crossed the boundary
  double tot_gas_mass_bd;     // Total gas mass that has crossed the boundary in the last iteration

  //gas mass variables
  double* gas_mass_bdnow;
  double* gas_mass;
  double* gas_mass_bd;
  double* gas_mass_start;

  //solid mass variables
  double* sol_mass;
  double* sol_mass_start;

  //file variables
  FILE* fb;
  char name[80];
  bool boundaries;
  int print_iter;
};

struct MassBalances mb = {0};

static double face_gradient_bid (Point point, scalar Y, int bid) {
  double grad = 0.;
  switch (bid) {
    case 0: grad = (Y[1,0]  - Y[])/Delta;  break;  // right
    case 1: grad = (Y[-1,0] - Y[])/Delta;  break;  // left
    case 2: grad = (Y[0,1]  - Y[])/Delta;  break;  // top
    case 3: grad = (Y[0,-1] - Y[])/Delta;  break;  // bottom
  }
  return grad;
}

static double face_value_bid (Point point, scalar Y, int bid) {
  double val = 0.;
  switch (bid) {
    case 0: val = 0.5*(Y[1,0]  + Y[]); break;  // right
    case 1: val = 0.5*(Y[-1,0] + Y[]); break;  // left
    case 2: val = 0.5*(Y[0,1]  + Y[]); break;  // top
    case 3: val = 0.5*(Y[0,-1] + Y[]); break;  // bottom
  }
  return val;
}

static double face_flux_bid (Point point, scalar Y, int bid, bool unity=false) {
  double flux = 0., faceval = (unity) ? 1. : face_value_bid (point, Y, bid);
  switch (bid) {
    case 0: flux = +faceval*uf.x[]*Delta; break;  // right
    case 1: flux = -faceval*uf.x[]*Delta; break;  // left
    case 2: flux = +faceval*uf.y[]*Delta; break;  // top
    case 3: flux = -faceval*uf.y[]*Delta; break;  // bottom
  }
  return flux;
}

static void diffusion_boundary (Point point, int bid) {
#ifdef MULTICOMPONENT
  double ff = face_value_bid (point, f, bid);
  double jG[NGS], jGtot = 0.;
  for (int jj = 0; jj < NGS; jj++) {
    if (ff < 1.-F_ERR) {
      scalar YG = YGList_G[jj];
      double gradYG = face_gradient_bid (point, YG, bid);
      scalar Dmix2  = DmixGList_G[jj];
      double Dmix2v = Dmix2[];
      
      double rhoGh;
  #ifdef VARPROP
      rhoGh = face_value_bid(point, rhoGv_G, bid);
  #else
      rhoGh = rhoG;
  #endif

      jG[jj] = rhoGh*Dmix2v*gradYG;
      #if FICK_CORRECTED
      jGtot += jG[jj];
      #endif
    }
  }

  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_G[jj];
    jG[jj] -= jGtot*face_value_bid(point, YG, bid);
  #if AXI
      mb.gas_mass_bdnow[jj] -= jG[jj]*Delta*dt*y;
  #else
      mb.gas_mass_bdnow[jj] -= jG[jj]*Delta*dt;
  #endif
  }

  //reset the aux arrays
  for (int jj = 0; jj < NGS; jj++) {
    jG[jj] = 0.;
  }
  jGtot = 0.;

  for (int jj = 0; jj < NGS; jj++) {
    if (ff > F_ERR) {
      scalar YG = YGList_S[jj];
      double gradYG = face_gradient_bid (point, YG, bid);
      scalar Dmix2  = DmixGList_S[jj];
      double Dmix2v = Dmix2[];

      double rhoGh;
  #ifdef VARPROP
      rhoGh = face_value_bid(point, rhoGv_S, bid);
  #else
      rhoGh = rhoG;
  #endif
      jG[jj] = rhoGh*Dmix2v*gradYG;
  #if FICK_CORRECTED
      jGtot += jG[jj];
  #endif
    }
  }

  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_S[jj];
    jG[jj] -= jGtot*face_value_bid(point, YG, bid);
  #if AXI
      mb.gas_mass_bdnow[jj] -= jG[jj]*Delta*dt*y*ff;
  #else
      mb.gas_mass_bdnow[jj] -= jG[jj]*Delta*dt*ff;
  #endif
  }
#endif
}

scalar U[];
U[right] = dirichlet(1.);
U[top] = dirichlet(1.);
U[left] = dirichlet(1.);
U[bottom] = dirichlet(1.);

static void advection_boundary (Point point, int bid) {
#ifdef MULTICOMPONENT
  foreach_elem (YGList_G, jj) {
    scalar YG = YGList_G[jj];
    double fluxYG = face_flux_bid (point, YG, bid);

    double rhoGh;
#ifdef VARPROP
    rhoGh = face_value_bid(point, rhoGv_G, bid);
#else
    rhoGh = rhoG;
#endif

    mb.gas_mass_bdnow[jj] += rhoGh*fluxYG*dt;
  }
  foreach_elem (YGList_S, jj) {
    scalar YG = YGList_S[jj];
    double fluxYG = face_flux_bid (point, YG, bid);
    double rhoGh;
#ifdef VARPROP
    rhoGh = face_value_bid(point, rhoGv_S, bid);
#else
    rhoGh = rhoG;
#endif

    mb.gas_mass_bdnow[jj] += rhoGh*fluxYG*dt;
  }
#endif

  double rhoGh;
#ifdef VARPROP
  rhoGh = face_value_bid(point, rhoGv_G, bid);
#else
  rhoGh = rhoG;
#endif
  double ff = face_value_bid (point, f, bid);
  mb.tot_gas_mass_bdnow += rhoGh*face_flux_bid (point, U, bid, unity=true)*(1.-ff)*dt; //careful about the flux in the pseudophase
}

static void write_balances(void) {
  //print total mass
  fprintf(mb.fb, "%-18.12f %-18.12f %-18.12f %-18.12f ", t,
                                        mb.tot_sol_mass,
                                        mb.tot_gas_mass,
                                        mb.tot_mass);

#ifdef MULTICOMPONENT
  //print individual gas species mass
  for (int jj=0; jj<NGS; jj++)
    fprintf(mb.fb, "%-18.12f", mb.gas_mass[jj]);
  
  //print individual solid species mass
  for (int jj=0; jj<NSS; jj++)
    fprintf(mb.fb, "%-18.12f", mb.sol_mass[jj]);
#endif

  fprintf(mb.fb, "\n");
  fflush(mb.fb);
}

static void compute_initial_state (void) {
  if (!mb.tot_gas_mass_start && !mb.tot_sol_mass_start) {
    foreach (serial) {
      mb.tot_sol_mass_start += (f[] - porosity[])*rhoS*dv();

      //Internal gas mass
      if (f[] > F_ERR) {
        double rhoGh;
        #ifdef VARPROP
        rhoGh = rhoGv_S[];
        #else
        rhoGh = rhoG;
        #endif
        mb.tot_gas_mass_start += porosity[]*rhoGh*dv(); //ef
      }
      //External gas mass
      if (f[] < 1.-F_ERR) {
        double rhoGh;
        #ifdef VARPROP
        rhoGh = rhoGv_G[];
        #else
        rhoGh = rhoG;
        #endif
        mb.tot_gas_mass_start += (1-f[])*rhoGh*dv(); //(1-f)
      }
    }

    mb.tot_mass_start = mb.tot_sol_mass_start + mb.tot_gas_mass_start;
    
    /**
    We have the option to use the analytical expressions for the initial 
    masses in the case of a spherical particle.
    */
    #ifdef BALANCES_SPHERE
    mb.tot_sol_mass_start = 3.14159265358979323846*sq(0.5)/4*(1-eps0)*rhoS;
    mb.tot_gas_mass_start = 3.14159265358979323846*sq(0.5)/4*(eps0-1)*rhoG + L0*L0*rhoG;
    mb.tot_mass_start = mb.tot_sol_mass_start + mb.tot_gas_mass_start;
    #endif

  #ifdef MULTICOMPONENT
    for (int jj=0; jj<NGS; jj++)
      mb.gas_mass_start[jj] = mb.tot_gas_mass_start*gas_start[jj];
    for (int jj=0; jj<NSS; jj++)
      mb.sol_mass_start[jj] = mb.tot_sol_mass_start*sol_start[jj];
  #endif
  } else {
    fprintf(stderr,"WARNING: Initial state already computed\n");
  }
}

static void compute_balances(void) {
  mb.tot_sol_mass = 0.;
  mb.tot_gas_mass = 0.;
  mb.tot_mass = 0.;
  mb.tot_gas_mass_bdnow = 0.;
  
#ifdef MULTICOMPONENT
  for (int jj=0; jj<NGS; jj++) {
    mb.gas_mass[jj] = 0.;
    mb.gas_mass_bdnow[jj] = 0.;
  }

  for (int jj=0; jj<NSS; jj++)
    mb.sol_mass[jj] = 0.;
#endif

  foreach (serial) {
    //compute total mass
    mb.tot_sol_mass += (f[] - porosity[])*rhoS*dv(); //(f-ef)
    
    if (f[] > F_ERR) {
      double rhoGh;
      #ifdef VARPROP
      rhoGh = rhoGv_S[];
      #else
      rhoGh = rhoG;
      #endif
      mb.tot_gas_mass += porosity[]*rhoGh*dv(); //ef
    }
    //External gas mass
    if (f[] < 1. - F_ERR) {
      double rhoGh;
      #ifdef VARPROP
      rhoGh = rhoGv_G[];
      #else
      rhoGh = rhoG;
      #endif
      mb.tot_gas_mass += (1. - f[])*rhoGh*dv(); //(1-f)
    }

#ifdef MULTICOMPONENT
    for (int jj = 0; jj < NGS; jj++) {
      scalar YG_G = YGList_G[jj];
      scalar YG_S = YGList_S[jj];

      if (f[] > F_ERR) {
        double rhoGh;
        #ifdef VARPROP
        rhoGh = rhoGv_S[];
        #else
        rhoGh = rhoG;
        #endif
        mb.gas_mass[jj] += porosity[]*rhoGh*(YG_S[]/f[])*dv(); //ef/f*(YG*f) = ef*YG
      }
     
      if (f[] < 1. - F_ERR) {
        double rhoGh;
        #ifdef VARPROP
        rhoGh = rhoGv_G[];
        #else
        rhoGh = rhoG;
        #endif
        mb.gas_mass[jj] += (1. - f[])*rhoGh*(YG_G[]/(1. - f[]))*dv(); //(1-ef)/(1-f)*(YG*(1-f)) = (1-ef)*YG
      }
    }

    //compute individual solid species mass
    if (f[] > F_ERR) {
      for (int jj = 0; jj < NSS; jj++) {
        scalar YS = YSList[jj];
        mb.sol_mass[jj] += (f[] - porosity[])*rhoS*(YS[]/f[])*dv(); // (f-ef)*(YS*f/f) = f*(1-e)*YS
      }
    }
#endif
  }

  if (mb.boundaries) {
#ifdef MULTICOMPONENT
    boundary(YGList_G); //do not remove
    boundary(YGList_S); //do not remove
# ifdef VARPROP
    boundary({rhoGv_G, rhoGv_S});
# endif
#endif
    boundary({U});
    for (int b=0; b<nboundary; b++) {
      foreach_boundary (b) {
        diffusion_boundary (point, b);
        advection_boundary (point, b);
      }
    }

#ifdef MULTICOMPONENT
  for (int jj=0; jj<NGS; jj++)
    mb.gas_mass_bd[jj] += mb.gas_mass_bdnow[jj];
#endif

    mb.tot_gas_mass_bd += mb.tot_gas_mass_bdnow;
  }

#ifdef MULTICOMPONENT
  for (int jj=0; jj<NGS; jj++)
    mb.gas_mass[jj] = ((mb.gas_mass[jj] + mb.gas_mass_bd[jj]) - mb.gas_mass_start[jj])/mb.tot_sol_mass_start;

  for (int jj=0; jj<NSS; jj++)
    mb.sol_mass[jj] = mb.sol_mass[jj] / mb.tot_sol_mass_start;
#endif

  mb.tot_gas_mass = ((mb.tot_gas_mass + mb.tot_gas_mass_bd) -  mb.tot_gas_mass_start)/mb.tot_sol_mass_start;
  mb.tot_sol_mass = mb.tot_sol_mass/mb.tot_sol_mass_start;
  mb.tot_mass = mb.tot_gas_mass + mb.tot_sol_mass;
}

event defaults (i = 0) {
  mb.tot_gas_mass_bdnow = 0.;
  mb.tot_gas_mass = 0.;
  mb.tot_gas_mass_bd = 0.;
  mb.tot_mass = 0.;
  mb.tot_sol_mass = 0.;
  mb.tot_gas_mass_start = 0.;
  mb.tot_mass_start = 0.;
  mb.tot_sol_mass_start = 0.;
  mb.fb = NULL;
  
  mb.print_iter = 10;
  mb.boundaries = true;

  mb.gas_mass_start = NULL;
  mb.gas_mass_bd = NULL;
  mb.gas_mass_bdnow = NULL;
  mb.gas_mass = NULL;

  mb.sol_mass_start = NULL;
  mb.sol_mass = NULL;
  
  sprintf(mb.name, "balances-%d", maxlevel);
  
  foreach()
    U[] = 1.;
}

event init (i = 0) {
#ifdef MULTICOMPONENT
  mb.gas_mass_start = (double*) malloc(NGS*sizeof(double));
  mb.gas_mass_bd =    (double*) malloc(NGS*sizeof(double));
  mb.gas_mass_bdnow = (double*) malloc(NGS*sizeof(double));
  mb.gas_mass =       (double*) malloc(NGS*sizeof(double));

  for (int jj=0; jj<NGS; jj++) {
    mb.gas_mass_start[jj] = 0.;
    mb.gas_mass_bd[jj] = 0.;
    mb.gas_mass_bdnow[jj] = 0.;
    mb.gas_mass[jj] = 0.;
  }

  mb.sol_mass_start = (double*) malloc(NSS*sizeof(double));
  mb.sol_mass =       (double*) malloc(NSS*sizeof(double));
  
  for (int jj=0; jj<NSS; jj++) {
    mb.sol_mass_start[jj] = 0.;
    mb.sol_mass[jj] = 0.;
  }
#endif

  mb.fb = fopen(mb.name, "w");
  //total mass header
  fprintf(mb.fb, "%-18s %-18s %-18s %-18s ", "t(1)", "tot_sol_mass(2)", "tot_gas_mass(3)", "tot_mass(4)");
  
#ifdef MULTICOMPONENT
  int counter = 5;
  //individual gas species mass header
  for (int jj=0; jj<NGS; jj++) {
    char buffer[80];
    sprintf(buffer, "mG_%s(%d)", OpenSMOKE_NamesOfSpecies(jj), counter++);
    fprintf(mb.fb, "%-18s", buffer);
  }

  //solid species mass header
  for (int jj=0; jj<NSS; jj++) {
    char buffer[80];
    sprintf(buffer, "mS_%s(%d)",  OpenSMOKE_NamesOfSolidSpecies(jj), counter++);
    fprintf(mb.fb, "%-18s", buffer);
  }
#endif

  fprintf(mb.fb, "\n");
  fflush(mb.fb);
  
  compute_initial_state();
}

event cleanup(t = end) {
  if (mb.fb != NULL)
    fclose(mb.fb);

  free(mb.gas_mass_start);
  free(mb.gas_mass_bd);
  free(mb.gas_mass_bdnow);
  free(mb.gas_mass);
  free(mb.sol_mass_start);
  free(mb.sol_mass);
}

event chemistry (i++) {

  compute_balances();

  if (i % mb.print_iter == 0)
    write_balances();
}