/**
# Coupling interface for master/slave simulations

This is used in combination with [slave.h]().

The user interface is just `slave_start` which defines the starting
time of the slave simulation and the interface of coupling functions:
here only `slave_interpolate()`. */

double slave_start = 0.;

extern double slave_interpolate (const char * name, double xp = 0, double yp = 0, double zp = 0,
				 bool linear = false);

/**
## Synchronization events */

event slave (i++) {
  extern void slave_step (double t);  
  slave_step (t + slave_start);
}

event slave_end (t = end) {
  extern void slave_stop();
  slave_stop();
}
