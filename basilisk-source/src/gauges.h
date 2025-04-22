/**
# Tide gauges

An array of *Gauge* structures passed to *output_gauges()* will create
a file (called *name*) for each gauge. Each time *output_gauges()* is
called a line will be appended to the file. The line contains the time
and the value of each scalar in *list* in the cell containing
*(x,y)*. The *desc* field can be filled with a longer description of
the gauge. */

typedef struct {
  char * name;
  double x, y;
  char * desc;
  FILE * fp;
} Gauge;

trace
void output_gauges (Gauge * gauges, scalar * list)
{
  int n = 0;
  for (Gauge * g = gauges; g->name; g++, n++);
  coord a[n];
  n = 0;
  for (Gauge * g = gauges; g->name; g++, n++) {
    double xp = g->x, yp = g->y;
    unmap (&xp, &yp);
    a[n].x = xp, a[n].y = yp;
  }
  int len = list_len(list);
  double v[n*len];
  interpolate_array (list, a, n, v, false);

  if (pid() == 0) {
    n = 0;
    for (Gauge * g = gauges; g->name; g++) {
      if (!g->fp) {
	g->fp = fopen (g->name, "w");
	if (g->desc)
	  fprintf (g->fp, "%s\n", g->desc);
      }
      if (v[n] != nodata) {
	fprintf (g->fp, "%g", t);
	for (scalar s in list)
	  fprintf (g->fp, " %g", v[n++]);
	fputc ('\n', g->fp);
	fflush (g->fp);
      }
      else
	n += len;
    }
  }
}
