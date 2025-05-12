/** Introduce a field for concentration*/

event defaults (i = 0)
{
  assert (cl == NULL);
  assert (nl > 0);
  for (int l = 0; l < nl; l++) {
    scalar c = new scalar;
    cl = list_append (cl, c);
  }
  reset (cl, 0.);
  int l = 0;
  for (scalar c in cl) {
    tracers[l] = list_append (tracers[l], c);
    l++;
  }
}
