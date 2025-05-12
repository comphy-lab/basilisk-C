event movie(t += 0.005){
  char fname[99];
  for (scalar s in{f}){
    view(camera = "front");
    draw_vof("f", lc = {1, 0, 0}, lw = 3.);
    squares(s.name, linear = false, map = cool_warm, spread = -1);
    sprintf(fname, "front_%s_U%g_f%gHz_l%d.mp4", s.name, force.V0, force.freq0, MAXLEVEL);
    save(fname);
  }
}


event time_series(t += 0.005){
  double ekin = 0., epot = 0.;
  foreach (reduction(+ : ekin), reduction(+ : epot)){
    epot += dv() * rhov[] * y;
    foreach_dimension(){
      ekin += dv() * 0.5 * rhov[] * sq(u.x[]);
    }
  }

  if (pid() == 0){
    FILE *fp = fopen("globals.asc", "a");
    fprintf(fp, "%f %.9g %.9g \n", t, ekin, epot);
    fclose(fp);
  }
}

event save_profiles(t += 0.1){

  scalar ek[], ff[], cc[], epsv[];
  foreach(){
    ek[] = rhov[]*sq(u.x[]) + sq(u.y[]);
    ff[] = sq(f[]);
    cc[] = sq(rhov[]);
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);    
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double Sxx = dudx;
    double Sxy = 0.5*(dudy + dvdx);
    double Syx = Sxy;
    double Syy = dvdy;
    epsv[] = 2.*mu1/rhov[]*(sq(Sxx) + sq(Sxy) + sq(Syx) + (Syy))*dv();
  }

  int n = (1 << MAXLEVEL);
  scalar * list = {f,rhov,ek,ff,cc,epsv};
  profile_foreach_region(list, xmin=-L0/2., xmax=L0/2., hmin=-L0/2.+_mindel, hmax=L0/2.-_mindel, n=n, m1=n);
}


event save_facets(t += 0.1){
  scalar kappa[];
  curvature(f, kappa);

  char fname[99];
  sprintf(fname, "interface_t%.5f", t);
  output_facets_vtu(f, kappa, fname);
}

event save_snapshots(t += 0.1)
{
  char filename[99];
  sprintf(filename, "snapshot_t%.5f", t);
  output_xmf({f, p}, {u}, filename);
}

extern scalar *interfaces;

event properties(i = 0)
{
  for (scalar c in interfaces)
    if (c.height.x.i)
      heights(c, c.height);
}

event tracer_advection(i++)
{
  for (scalar c in interfaces)
    if (c.height.x.i)
      heights(c, c.height);
}