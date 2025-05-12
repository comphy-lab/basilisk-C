#if TREE
static void momentum_refine (Point point, scalar u) {
  refine_bilinear (point, u);
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rho(f2[],f3[])*u[];
  double du = u[] - rhou/((1 << dimension)*cm[]*rho(f2[],f3[]));
  foreach_child()
    u[] += du;
}

static void momentum_restriction (Point point, scalar u)
{
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rho(f2[],f3[])*u[];
  u[] = rhou/((1 << dimension)*cm[]*rho(f2[],f3[]));
}
#endif // TREE

event defaults (i = 0)
{
  stokes = true;

#if TREE

  int i = 0;
  while (all[i].i != f1.i) i++;
  while (all[i].i)
    all[i] = all[i-1], i--;
  all[i] = f1;

  i = 0;
  while (all[i].i != f2.i) i++;
  while (all[i].i)
    all[i] = all[i-1], i--;
  all[i] = f2;

  i = 0;
  while (all[i].i != f3.i) i++;
  while (all[i].i)
    all[i] = all[i-1], i--;
  all[i] = f3;
  foreach_dimension() {
    u.x.refine = u.x.prolongation = momentum_refine;
    u.x.restriction = momentum_restriction;
  }
#endif
}


event stability (i++)
  dtmax = timestep (uf, dtmax);

foreach_dimension()
static double boundary_q2_x (Point neighbor, Point point, scalar q2, void * data)
{
  return clamp(f2[],0.,1.)*rho2*u.x[];
}

foreach_dimension()
static double boundary_q3_x (Point neighbor, Point point, scalar q3, void * data)
{
  return clamp(f3[],0.,1.)*rho3*u.x[];
}

foreach_dimension()
static double boundary_q1_x (Point neighbor, Point point, scalar q1, void * data)
{
  return clamp(1.-f2[]-f3[],0,1)*rho1*u.x[];
}

#if TREE
foreach_dimension()
static void prolongation_q2_x (Point point, scalar q2) {
  foreach_child()
    q2[] = clamp(f2[],0.,1.)*rho2*u.x[];
}

foreach_dimension()
static void prolongation_q3_x (Point point, scalar q3) {
  foreach_child()
    q3[] = clamp(f3[],0.,1.)*rho3*u.x[];
}

foreach_dimension()
static void prolongation_q1_x (Point point, scalar q1) {
  foreach_child()
    q1[] = clamp(1.-f2[]-f3[],0,1)*rho1*u.x[];
}
#endif

static scalar * interfaces1 = NULL;

event vof (i++) {
  vector q1[], q2[], q3[];
  for (scalar s in {q1,q2,q3})
    foreach_dimension()
      s.v.x.i = -1; // not a vector
  for (int i = 0; i < nboundary; i++)
    foreach_dimension() {
      q1.x.boundary[i] = boundary_q1_x;
      q2.x.boundary[i] = boundary_q2_x;
      q3.x.boundary[i] = boundary_q3_x;
    }
#if TREE
  foreach_dimension() {
    q1.x.prolongation = prolongation_q1_x;
    q2.x.prolongation = prolongation_q2_x;
    q3.x.prolongation = prolongation_q3_x;
  }
#endif

  foreach()
    foreach_dimension() {
      double fc2 = clamp(f2[],0,1);
      double fc3 = clamp(f3[],0,1);
      double fc1 = clamp(1.-f2[]-f3[],0,1);            
      q2.x[] = fc2*rho2*u.x[];
      q3.x[] = fc3*rho3*u.x[];      
      q1.x[] = fc1*rho1*u.x[];
    }
  boundary ((scalar *){q1,q2,q3});

  foreach_dimension() {
//    q2.x.inverse = true;
    q1.x.inverse = false;
    q2.x.inverse = false;
    q3.x.inverse = false;            
    q1.x.gradient = q2.x.gradient = q3.x.gradient =  u.x.gradient;
  }

//  f.tracers = (scalar *){q1,q2};
  f1.tracers = (scalar *){q1};
  f2.tracers = (scalar *){q2};
  f3.tracers = (scalar *){q3};
      
  vof_advection ({f2,f3,f1}, i);

  foreach()
    foreach_dimension()
      u.x[] = (q1.x[] + q2.x[] + q3.x[])/rho(f2[],f3[]);
  boundary ((scalar *){u});

  interfaces1 = interfaces, interfaces = NULL;
}

event tracer_advection (i++) {
  interfaces = interfaces1;
}