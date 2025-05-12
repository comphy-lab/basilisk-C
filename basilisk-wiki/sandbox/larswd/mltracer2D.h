/**
2D algorithm for transport of tracer 

The general structure is: 

1. Initializing the linear system
2. Solve linear system
3. Check for convergence 
4. Cleanup

The initialization requires us declaring a lot of arrays, which we initialize once at the start. Some of these are likely superflous, and will be cleaned up in a future patch. 
*/

Ml_Tracer transport_tracer(Ml_Tracer _i_tracer){
  // Initialization
  double tol = tml().tol;
  scalar c   = tml().c;
  double rho_c = tml().rho_c;
  scalar c_0 = new scalar[nl];
  scalar c_prev = new scalar[nl];
  double a_c[N], a_w[N], a_e[N],b[N], c_[N];
  double total_c0 = 0;
  TDMA_fluxes cf;
  foreach(serial){
    foreach_layer() {
      c_0[] = c[];
      c_prev[] = c[];
      total_c0 += Delta*Delta*h[]*sq(c[]);
    }
  }
  printf("C total before = %g\n", total_c0);
  double residual = HUGE;
  int iterations = 0;

  for (int i = 0; i < N; i++){
    a_c[i]  = 0.0;
     a_w[i]  = 0.0;
    a_e[i]  = 0.0;
    b[i]    = 0.0;
    c_[i]   = 0.0;
  }
  while (residual > tol){
    residual = 0;
    foreach_point(0,0, nowarning){
      for (int l = 0; l < nl; l++){
        // Initializing Linear system
        for (point.i = GHOSTS; point.i < point.n + GHOSTS; point.i++){
          cf = compute_fluxes(point, l, _i_tracer, c_0);
          int idx = point.i - GHOSTS;
          a_w[idx]  = - cf.aW;
          a_c[idx]  =  rho_c/dt + cf.aW + cf.aE + cf.aT + cf.aB + cf.div - cf.Sp;
          a_e[idx]  = - cf.aE;
          b[idx]    = rho_c/dt * c_0[0,0,l] + cf.Su;
          if (l < nl-1)
            b[idx] += cf.aT*c_prev[0,0,l+1];
          if (l > 0)
            b[idx] += cf.aB*c_prev[0,0,l-1];
          double Sdc = TVD_deferred_correction(point, c_prev, l, _i_tracer);
          b[idx] += Sdc/Delta;
        }
        //sleep(3);
        TDMA(a_w, a_c, a_e, b, c_, point.n);
        for (point.i = GHOSTS; point.i < point.n + GHOSTS; point.i++){
          c[0,0,l] = c_[point.i - GHOSTS];
        }
      }
    }

    double total_c = 0;
    foreach(serial){
      foreach_layer(){
        residual += sq(c[] - c_prev[])*Delta*Delta*h[];
        total_c += Delta*Delta*h[]*sq(c[]);
      }
    }
    foreach(){
      foreach_layer()
        c_prev[] = c[];
    }
    printf("total_c = %g residual = %g.\n",total_c, residual);
    if (nl > 1)
      residual = sqrt(residual);
    else
      residual = 0;
    iterations++;
  }
  tml().c = c;
  printf("Solution converged after %d iterations\n", iterations);
  delete((scalar *){c_0, c_prev});
  return _i_tracer;
}