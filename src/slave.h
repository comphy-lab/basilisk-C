/**
# Coupled slave 

This can be used in combination with [master.h]() to couple two
Basilisk solvers.

## Example of coupling function

This is a simple function which can be called by the master to
retrieve (interpolated) field values from the slave. 

Using this model, it should be easy to write more sophisticated
coupling functions, if necessary. 

Note that new coupling functions must be added to the list of symbols
which are exported. See the recipe for `slave.o` in the examples
[Makefile](/src/examples/Makefile) for details. */

double slave_interpolate (const char * name, double xp, double yp, double zp, bool linear)
{
  if (!grid) {
    fprintf (stderr, "slave_interpolate: error: no grid! this may be a "
	     "master/slave synchronization issue\n");
    exit (1);
  }
  scalar s = lookup_field (name);
  if (s.i < 0) {
    fprintf (stderr, "slave_interpolate: error: unknown field '%s'\n", name);
    exit (1);
  }
  return interpolate (s, xp, yp, zp, linear);
}

/**
## Synchronization functions

These functions are called by the master to synchronize the
master/slave times. */

#if _OBJECT
event slave_event (t <= HUGE);

void slave_step (double t0)
{
  static Event * slave = NULL;
  if (!slave) {
    _init_solver();
    int slave_init();
    slave_init();
    init_grid (N);
  
    for (Event * ev = Events; !ev->last; ev++)
      if (!strcmp (ev->name, "slave_event")) {
	slave = ev;
	break;
      }
    assert (slave);

    iter = 0, t = 0., dt = 1.;
    events (true);
    iter = inext, t = tnext;
  }

  slave->t = t0;
  while (t0 - t > TEPS*t0 && events (true))
    iter = inext, t = tnext;
}

void slave_stop()
{
  free_grid();
}

/**
When this file is included in a simulation (i.e. compiled as an
object), the main() and run() functions are overloaded, since the
master main() function will drive the simulation. */

# define main() slave_init()
# define run() return 0
#endif // _OBJECT

