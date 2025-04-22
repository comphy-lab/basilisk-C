event event1 (t = 4; t <= 12.23; t += 0.1)
{
  fprintf (stderr, "t: %g i: %d\n", t, i);
}

event event2 (i = 5) fprintf (stderr, "i: %d\n", i);

event event3 (i = 2; i *= 2) fprintf (stderr, "i^2: %d\n", i);

event event4 (i++) fprintf (stderr, "i++: %d\n", i);

event event5 (t += 1; t <= 20) fprintf (stderr, "t++: %g\n", t);

event event6 (t = end) fprintf (stderr, "end: %g\n", t);

event event7 (t = 10; t <= 12; i++)
  fprintf (stderr, "event7: %d %g\n", i, t);

event event8 (i = end)
  fprintf (stderr, "event8: %d %g\n", i, t);

int main (int argc, char * argv[])
{
  iter = 0, t = 0.;
  while (events (true))
    iter = inext, t = tnext;
}
