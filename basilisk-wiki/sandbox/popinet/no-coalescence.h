/**
# Avoiding coalescence of VOF interfaces

When two interfaces defined by the same VOF tracer are close enough
(of the order of the cell size), they automatically merge. This is one
of the strength and weakness of the VOF method.

In some cases, it may be desirable to avoid coalescence entirely, for
example in the case of foams, emulsions, bubble clouds etc...

A simple way to do this is to use a different VOF tracer for each
bubble/droplet. When one wants to simulate more than a few bubbles,
this can of course become very expensive (both in CPU and memory).

This simple idea can be improved by noting that it would be sufficient
to use different VOF tracers only for bubbles which are "too close" to
one another. Determining the minimum number of VOF tracers required,
for a given arrangement of bubbles, is clearly a variant of the [graph
coloring](https://en.wikipedia.org/wiki/Graph_coloring) problem. In
two dimensions, the famous [four color
theorem](https://en.wikipedia.org/wiki/Four_color_theorem) states that
a maximum of four VOF tracers are required. Note however that finding
this optimal coloring can be very difficult
([NP-complete](https://en.wikipedia.org/wiki/Graph_coloring#Computational_complexity)). The
important point here is that one can expect that even a non-optimal
number of VOF tracers will be much smaller than the number of bubbles.

## User interface

This file is typically combined with the [two-phase
solver](/src/two-phase.h). From the user point-of-view, the only thing
to be aware of is that the default *f* volume fraction field is not
transported using VOF anymore, but is the sum of all VOF
tracers. Individual VOF tracers are named *f0*, *f1*, *f2*, ... and
are stored in the *interfaces* list (defined by
[vof.h](/src/vof.h)). They should be used in particular to display the
actual interfaces, as displaying interfaces using *f* will result in
coalescence artefacts.

## Utility functions

We will need to "tag" individual bubbles. EPS is the threshold used
for tagging. */

#include "tag.h"
#define EPS 1e-6

/**
This function returns a renamed clone of the default volume fraction
field *f*, with *i* its index. */

static scalar fclone (int i)
{
  scalar c = new scalar;
  scalar_clone (c, f);
  free (c.name);
  char s[80];
  sprintf (s, "%s%d", f.name, i);
  c.name = strdup (s);
  return c;
}

/**
This function does a fast detection of cases which may correspond to
two interfaces being close to one another. */

static bool tracer_is_close (Point point, scalar c)
{
  if (c[] > EPS)
    return false;
  foreach_dimension()
    if (c[-1] > EPS && c[1] > EPS)
      return true;
  if (c[-1,-1] > EPS && c[1,1] > EPS)
    return true;
  if (c[-1,1] > EPS && c[1,-1] > EPS)
    return true;
  return false;
}

/**
This is similar to the function above, but now takes into account
whether the two interfaces belong to different bubbles (identified by
the tag field *b*). If they do then the indices of the two bubbles are
returned in *b1* and *b2*. */

static bool bubbles_are_close (Point point, scalar c, scalar b,
			       int * b1, int * b2)
{
  if (c[] > EPS)
    return false;
  foreach_dimension()
    if (c[-1] > EPS && c[1] > EPS && b[-1] != b[1]) {
      *b1 = b[-1] - 1; *b2 = b[1] - 1;
      return true;
    }
  if (c[-1,-1] > EPS && c[1,1] > EPS && b[-1,-1] != b[1,1]) {
    *b1 = b[-1,-1] - 1; *b2 = b[1,1] - 1;
    return true;
  }
  if (c[-1,1] > EPS && c[1,-1] > EPS && b[-1,1] != b[1,-1]) {
    *b1 = b[-1,1] - 1; *b2 = b[1,-1] - 1;
    return true;
  }
  return false;
}

#if _MPI
static void reduce_bubbles_op (void * pin, void * pout, int * len,
			       MPI_Datatype * dptr)
{
  int * in = pin, leni, * out = pout, leno;
  for (leni = 0; leni < *len && in[leni] >= 0; leni++);
  for (leno = 0; leno < *len && out[leno] >= 0; leno++);
  int add = leno;
  for (int i = 0; i < leni; i += 2) {
    bool found = false;
    for (int j = 0; j < leno && !found; j++)
      if (in[i] == out[j] || in[i + 1] == out[j])
	found = true;
    if (!found) {
      assert (add < *len);
      out[add++] = in[i];
      assert (add < *len);
      out[add++] = in[i + 1];
    }
  }
}

trace
static void reduce_bubbles (Array * tc)
{
  if (npe() > 1) {
    int len1 = tc->len/sizeof(int), len = len1;
    mpi_all_reduce (len, MPI_INT, MPI_SUM);
    if (len > 0) {
      tc->max = len*sizeof(int);
      tc->p = realloc (tc->p, tc->max);
      for (int i = len1; i < len; i++)
	((int *)tc->p)[i] = -1;
      MPI_Op op;
      MPI_Op_create (reduce_bubbles_op, false, &op);
      MPI_Allreduce (MPI_IN_PLACE, tc->p, len, MPI_INT, op, MPI_COMM_WORLD);
      MPI_Op_free (&op);
      for (len1 = 0; len1 < len && ((int *)tc->p)[len1] >= 0; len1++);
      tc->len = len1*sizeof(int);
    }
  }
}
#endif // _MPI

/**
## Algorithm */

trace
void no_coalescence()
{

  /**
  We first make a quick test to check which VOF tracers may correspond
  to bubbles which are too close to one another. This is essentially
  an optimisation which avoids calling the relatively expensive
  *tag()* function if it is obviously not necessary. */

  int nvar = datasize/sizeof(double), too_close[nvar];
  for (int i = 0; i < nvar; i++)
    too_close[i] = false;
  foreach_leaf() // no openMP
    for (scalar c in interfaces)
      if (tracer_is_close (point, c))
	too_close[c.i] = true;
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, too_close, nvar, MPI_INT, MPI_MAX,
		 MPI_COMM_WORLD);
#endif
  scalar * maybe_close = NULL;
  for (scalar c in interfaces)
    if (too_close[c.i])
      maybe_close = list_append (maybe_close, c);
  
  for (scalar c in maybe_close) {

    /**
    For each VOF tracer which may be too close, we first tag the
    corresponding bubbles. */
    
    scalar b[];
    foreach()
      b[] = c[] > EPS;
    tag (b);

    /**
    The next step is to build the array *tc* of the bubbles which are
    indeed too close to one another. */
    
    Array * tc = array_new();
    foreach_leaf() { // no openMP
      int b1, b2;
      if (bubbles_are_close (point, c, b, &b1, &b2)) {
	for (int l = 0, * p = tc->p; l < tc->len/sizeof(int); l++, p++)
	  if (*p == b1 || *p == b2) {
	    // one of the bubbles is already in the list
	    b1 = -1; break;
	  }
	// Add these bubbles to the list if they are not already there
	if (b1 != -1) {
	  array_append (tc, &b1, sizeof (int));
	  array_append (tc, &b2, sizeof (int));
	}
      }
    }
#if _MPI
    reduce_bubbles (tc);
#endif
       
    int len = tc->len/sizeof(int);
    if (len > 0) {

#if 1
      fprintf (ferr, "rep %g", t);
      for (int i = 0; i < len; i++)
	fprintf (ferr, " %d", ((int *)tc->p)[i]);
      fputc ('\n', ferr);
      fflush (ferr);
#endif
      
      /**
      ### Neighboring tracers
	 
      We need to know which tracer fields are neigboring each
      bubble. If the tracer of index *j* is neighboring the bubble
      of index *i* (in *tc*), then *adj[i*nvar + j]* is set to *true*. */
      
      int nvar = datasize/sizeof(double), adj[len*nvar];
      for (int i = 0; i < len*nvar; i++)
	adj[i] = false;
      
      /**
      Since we are updating *adj*, we cannot use openMP. */
      
      foreach_leaf() // no openMP
	if (b[])
	  for (int i = 0, * p = tc->p; i < len; i++, p++)
	    if (b[] == *p + 1)
	      
	      /**
	      We check whether bubble b[] neighbors cells containing
	      another tracer. */
		
	      foreach_neighbor()
		for (scalar s in interfaces)
		  if (s.i != c.i && s[] > EPS)
		    adj[i*nvar + s.i] = true;
      
#if _MPI
      MPI_Allreduce (MPI_IN_PLACE, adj, len*nvar, MPI_INT, MPI_MAX,
		     MPI_COMM_WORLD);
#endif
      
      /**
      ### Finding a replacement tracer

      If this is the first bubble we need to replace, then the only
      existing VOF interface is *f*. We create a new interface, add it
      to the list and remove *f* from the list of interfaces (since it
      is not advected by VOF anymore). */
      
      if (c.i == f.i) {
	scalar f1 = fclone (0);
	foreach()
	  f1[] = f[];
	boundary ({f1});
	interfaces = list_copy ({f});
	swap (char *, f.name, f1.name);
	f.i = f1.i;
      }
      
      /**
      Array *rep* will contain the index of the replacement VOF tracer
      for each bubble. */
      
      int rep[len/2];
      for (int i = 0, * p = tc->p; i < len; i += 2) {
	
	/**
	The indices of the pair of neighboring bubbles are stored in
	`p[i]` and `p[i+1]`. We need to replace only one of the two
	bubbles. We choose to replace the bubble with the smallest
	number of neighboring tracers. */

	int n1 = 0, n2 = 0, j = i;
	for (scalar s in interfaces)
	  if (s.i < nvar) {
	    if (adj[i*nvar + s.i]) n1++;
	    if (adj[(i + 1)*nvar + s.i]) n2++;
	  }
	if (n2 < n1) {
	  p[i] = p[i + 1];
	  j++;
	}
	
	/**
	We look for a replacement VOF tracer which is not already
	neighboring the bubble. */

	rep[i/2] = -1;
	for (scalar s in interfaces)
	  if (s.i != c.i && (s.i >= nvar || !adj[j*nvar + s.i])) {
	    rep[i/2] = s.i; break;
	  }

	/**
	If we didn't find any, we create a new one. */
	
	if (rep[i/2] < 0) {
	  scalar t = fclone (list_len (interfaces));
	  reset ({t}, 0.);
	  interfaces = list_append (interfaces, t);
	  rep[i/2] = t.i;
	}
      }

      /**
      ### Replacing tracers
      
      We perform the replacement for each bubble (which is too
      close) and apply boundary conditions. */
      
      foreach()
	for (int i = 0, * p = tc->p; i < len; i += 2, p += 2)
	  if (b[] == *p + 1) {
	    scalar t = {rep[i/2]};
	    t[] = c[]; c[] = 0.;
	  }
      scalar * list = list_copy ({c});
      for (int i = 0; i < len/2; i++)
	list = list_add (list, (scalar){rep[i]});
      boundary (list);
      free (list);
    }
    
    /**
    Finally, we free the arrays and lists. */
    
    array_free (tc);
  }
  free (maybe_close);
}

/**
## Coupling with the solver

We apply the no coalescence function just before VOF advection. */

event vof (i++)
  no_coalescence();

/**
After VOF advection, the default volume fraction field *f* is updated
as the sum of all VOF tracer fields. */

event tracer_advection (i++)
{
  foreach() {
    double fsum = 0.;
    for (scalar c in interfaces)
      fsum += c[];
    f[] = clamp(fsum,0,1);
  }  
  boundary ({f});
}

/**
We need to free the list of interfaces if it has been dynamically
allocated. */

event cleanup (i = end)
{
  if (interfaces[0].i != f.i) {
    free (interfaces);
    interfaces = NULL;
  }
}
