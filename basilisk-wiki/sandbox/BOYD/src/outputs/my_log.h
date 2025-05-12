event logfile(i++) {
  if (i == 0){
    fprintf(stderr, "t dt grid->tn/1e6 perf.t/3600 perf.speed\n");
  }
  fprintf(stderr, "%g %g %.2fe+06 cells %g hrs %g cells/s\n", t, dt,
         ((float)(grid->tn)/1.0e6), (float)(perf.t)/3600.0, perf.speed);
}

// this sets the simulation duration
event finish_simulation(i++, t <= T_END) {
}
