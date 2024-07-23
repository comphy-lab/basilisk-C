int main()
{
  init_grid (64);
  scalar a[];
  fprintf (stderr, "depth: %d mem: %ld\n", depth(), pmtrace.total);
  refine (level < 8);
  fprintf (stderr, "depth: %d mem: %ld\n", depth(), pmtrace.total);
  unrefine (level > 5);
  tree_check();
  long tnc = 0;
  foreach()
    tnc++;
  fprintf (stderr, "depth: %d leaves: %ld mem: %ld\n",
	   depth(), tnc, pmtrace.total);
}
