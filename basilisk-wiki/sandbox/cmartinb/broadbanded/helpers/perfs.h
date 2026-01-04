/**
# Performance monitoring (for the Navier--Stokes solvers)

This logs simple statistics available for the various [Navier--Stokes
solvers](/src/README#navierstokes). */

event perfs (i += 1) {
  if( pid() == 0) {
    fflush(stderr);
    char file[80];
    sprintf (file,"perfs.out");
    FILE * fp = fopen(file,"a");
    if (i == 0)
      fprintf (fp,
               "t iter dt mgp.i mgp.nrelax mgpf.i mgpf.nrelax mgu.i mgu.nrelax "
               "grid->tn perf.t perf.speed npe\n");
    fprintf (fp, "%.10e %d %.10e %d %d %d %d %d %d %ld %.10e %.10e %d\n", 
             t, i, dt, mgp.i, mgp.nrelax, mgpf.i, mgpf.nrelax, mgu.i, mgu.nrelax,
             grid->tn, perf.t, perf.speed, npe());
    fclose(fp);
  }
  if (from_pr == 1 && is_first) {
    i_sp = i;
    t_sp = t;
    is_first = false;
    if (pid() == 0) { 
      fflush(stderr);
      char name_1[80];
      sprintf (name_1,"info_sp.out");
      FILE * log_sim = fopen(name_1,"a");
      fprintf (log_sim, "%d %.10e %d %.10e\n", i, t, i_sp, t_sp);
      fclose(log_sim);
    }
  }
}
