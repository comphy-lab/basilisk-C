/**
Checks "undefined" scalar reductions. */

int main()
{
  init_grid (1);
  scalar f[];
  {
    interpreter_verbosity (2);
    double min = HUGE, max = - HUGE, sum = 0., sum2 = 0., volume = 0.;
    foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
	    reduction(max:max) reduction(min:min))
      if (dv() > 0.) {
	volume += dv();
	if (f[] > max) max = f[];
      }
    if (sqrt(volume));
    if (max);
  }
}
