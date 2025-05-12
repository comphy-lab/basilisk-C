/**
*   ##Continuity waves computation
*
* Here we compute the quantity $\alpha = \int_S \chi_d dS$
*
*/

void compute_alpha_layer(){
  int nL = round(Ls/(D/LDn)) + 1;
  double alpha[nL];
  double surf[nL];
  double surff[nL];
  for (int i = 0; i < (nL); i++)
    alpha[i]=surf[i]=surff[i]=0;

  foreach(reduction(+:alpha[:nL]) reduction(+:surf[:nL])
          reduction(+:surff[:nL])){
    int idx = round(  (y+Ls/2.)  /   (D/LDn) );
    alpha[idx] += f[]* dv();
    surf[idx] += dv();
    surff[idx] += dv()*(1 - f[]);
  }

  if(main_process()){
    falpha = fopen("falpha.csv","a");
    fprintf(falpha,"%g",t);
    for (int j = 0; j < (nL-1); j++)
      if(surf[j])
        fprintf(falpha,",%g",alpha[j]/surf[j]);
      else
        fprintf(falpha,",%g",0.);
    fprintf(falpha,"\n");
    fclose(falpha);
  }
}

void Eulerian_nearest_stat(Drop n[],int step){
  coord dist_nst, vel_nst,vel_rel_nst, a_nst;
  double p_nst, f_nst;
  coord  g_vel_rel_nst,g_vel_nst, g_a_nst;
  double g_p_nst, g_f_nst;
  
  // We first find the nearest particle of the stencil 
  for (int i = 0; i < (N+1); i++)
  for (int j = 0; j < (N+1); j++)
  for (int k = 0; k < (N+1); k++){
    double x = Ls/(N+2) + i * Ls/(N+2)  - Ls/2;
    double y = Ls/(N+2) + j * Ls/(N+2)  - Ls/2;
    double z = Ls/(N+2) + k * Ls/(N+2)  - Ls/2;
    coord Pos = {x,y,z};
    double d = Ls;
    for (int l = 0; l < n[0].tagmax; l++){
      double r = dist_perio(Pos,n[l].Y);
      if(r < d){
        d = r;
        dist_nst = dist_perio_coord(Pos,n[l].Y);
        p_nst = interpolate(p,x,y,z);
        f_nst = interpolate(f,x,y,z);
        foreach_dimension(){
          vel_nst.x = interpolate(u.x,x,y,z);
          vel_rel_nst.x = interpolate(u.x,x,y,z) - n[l].U.x;
          a_nst.x = interpolate(a.x,x,y,z);
        }
      }
    }
    #if _MPI
    MPI_Allreduce(&p_nst, &g_p_nst, 1, MPI_DOUBLE, MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&f_nst, &g_f_nst, 1, MPI_DOUBLE, MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&vel_nst.x, &g_vel_nst.x, dimension, MPI_DOUBLE, MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&vel_rel_nst.x, &g_vel_rel_nst.x, dimension, MPI_DOUBLE, MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&a_nst.x, &g_a_nst.x, dimension, MPI_DOUBLE, MPI_MIN,MPI_COMM_WORLD);
    #else
    g_vel_nst = vel_nst;
    g_vel_rel_nst = vel_rel_nst;
    g_a_nst = a_nst;
    g_p_nst =  p_nst;
    g_f_nst = f_nst;
    #endif
  

    if(main_process()){
      fCAnst = fopen("fCAnst.csv","a");
      fprintf(fCAnst,"%d,%g,",step,t);
      einstein_sum(i){
        fprintf(fCAnst,"%g,",dist_nst.i);
        fprintf(fCAnst,"%g,",g_vel_nst.i);
        fprintf(fCAnst,"%g,",g_vel_rel_nst.i);
        fprintf(fCAnst,"%g,",g_a_nst.i);
      }
      fprintf(fCAnst,"%g,",g_f_nst);
      fprintf(fCAnst,"%g",g_p_nst);
      fprintf(fCAnst,"\n");
      fclose(fCAnst);
    }
  }

  //those must be array of my bins 
}
