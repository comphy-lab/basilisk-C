/**
#Controlling the coalescence using VOF

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
one another. Which is acheived by no-coalescence.h. 

In many cases there is coalescence, which can't be achieved by
no-coalescence.h. Eventhough using single VOF tarcer for all the
bubbles allow coalescence between them, to match with physical
coalescence time scales we might need to use fine mesh which will make
it computationally very expensive. Here by taking advantage of
no-coalescence.h, we can say when the bubbles are in contact for the
desired amount of time, then we will make the VOF tracers of these
bubbles in contact to be same VOF tracer, allowing them to coalesce
after the desired amount of time. Here on contrary to DNS we can use
coarse mesh to delay the coalescence.

## User interface

This file is typically combined with the [two-phase
solver](/src/two-phase.h) and [no-coalescence
file](/sandbox/mani/no-coalescence.h).

## Utility functions

We will need to "tag" individual bubbles or drops. EPS1 is the
threshold used for tagging. */

#include "tag.h"
#define EPS1 1e-6

/** We create a data structure that stores the details of the pairs
of bubbles that we are interetsed to track over time. */

typedef struct
{
  coord cv[2];
  int b;
  double t, ft;
  int vof1, vof2;
}clsncdat;

/** This function update the information of pairs of bubbles */

static void dat_change(clsncdat * a, clsncdat * datc)
{
  a->b = datc->b;
  a->cv[0].x = datc->cv[0].x;
  a->cv[0].y = datc->cv[0].y;
  a->cv[1].x = datc->cv[1].x;
  a->cv[1].y = datc->cv[1].y;
 
#if dimension > 2
  a->cv[0].z = datc->cv[0].z;
  a->cv[1].z = datc->cv[1].z;
#endif
  a->t = datc->t;
  a->ft += datc->ft;
  a->vof1 = datc->vof1;
  a->vof2 = datc->vof2;
}

/** This function lits the tracers corresponding to the bubbles that
    are one to two cell distance away
*/
static bool tracer_is_too_close (Point point, scalar a, scalar c)
{
  if(a.i != c.i && a[]>EPS1)
    foreach_neighbor(1)
      if(c[]>EPS1)
	return true;
   return false;
}
/** This function detects the pairs of bubbles that are one to two cell 
distance away from each other (and have different VOF tracers for each)
*/
static bool interfaces_in_cell(Point point, scalar a, scalar c, scalar T)
{
  int tag_val = T[];
  if(a.i != c.i && a[]>EPS1 && T[])
    foreach_neighbor(1)
      if(c[]>EPS1 && T[] == tag_val)
	return true;
  return false;
}
/** This function helps to create new VOF tracers */
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

/** Here cntrl_nclsnc function list the pairs of bubbles that came
close to each other (and are already processed by no-coalescence.h,
there by having different VOF tracer for each), track them over time
and if the time of contact (time for which they have stayed close,
less than one to two cell distance away from each other) exceeds the
prescribed time (film drainage time) then list those pairs in a new
list and remaining pairs in another list. The list that contains the
pairs to be coalesced is used to iterate through each pair of bubble
and replicate VOF tracer of one of its to other. At the end delete
this list and return the list that contains bubble pairs who stil
didn't reach the drainage time. 

##control coalescence function
*/

trace
Array * cntrl_nclsnc(Array * dat, double drainage_time)
{

  /** Detect the VOF tracers corresponding to pairs of bubbles that
      are one to two cell distance away*/
  
  int nvar = datasize/sizeof(double), too_close[nvar];
  for (int i = 0; i < nvar; i++)
    too_close[i] = false;
  foreach_leaf() // no openMP
    for (scalar a in interfaces)
      for (scalar c in interfaces)
	if (tracer_is_too_close (point, a, c)){
	  too_close[a.i] = true;
	  too_close[c.i] = true;
	}
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, too_close, nvar, MPI_INT, MPI_MAX,
		 MPI_COMM_WORLD);
#endif
  scalar * maybe_close = NULL;
  for (scalar c in interfaces)
    if (too_close[c.i])
      maybe_close = list_append (maybe_close, c);

  /** Here we tag the bubbles that are one to two cell distance away. Each
      pair of bubbles will have the same tag value. */

  scalar T = new scalar;  //fix me, when scalar T[] is used instead, it gives an error
  foreach(){
    T[]=0;
    for(scalar a in maybe_close)
      if(a[]>EPS1)
	T[] = 1;
  }
  tag(T);
  
  /** Here we detect the pairs of bubbles that have different VOF
      tracers and add them to the list bub */
  
  Array * bub = array_new();
  
  foreach_leaf() { // no openMP
    clsncdat b1;
    for(scalar a in maybe_close)
      for(scalar c in maybe_close)
	if(interfaces_in_cell(point, a, c, T))
	  {
	    b1.b = T[]-1;
	    b1.vof1 = a.i;
	    b1.vof2 = c.i;
	    
	    clsncdat * p = bub->p;
	    
	    for (int l = 0; l < (bub->len/sizeof(clsncdat)); l++)
	      {
		if ( p->b == b1.b) {
		  // one of the bubbles is already in the list
		  b1.b = -1;
		  break;
		}
		p++;
	      }
	    // Add these bubbles to the list if they are not already there
	    if (b1.b != -1) {
	      assert (b1.b >= 0);
	      array_append (bub, &b1, sizeof(clsncdat));		  	      
	    }
	  }
  }  
  int length = bub->len/sizeof(clsncdat);
#if 1
  if(length)
    fprintf(stderr,"length = %d, time = %g\n", length, t);
#endif
  
  scalar T1 = new scalar;
  foreach()
    T1[] = 0;

  /** Here we store the centroid and VOF tracer information for each
      pair of the bubbles */
  clsncdat v1[length];
  clsncdat * p1 = bub->p;
  for(int j = 0; j < length; j++){
    foreach_leaf()
      if(T[] == (p1->b) + 1)
	T1[] = T[]+j;
    v1[j].cv[0].x = v1[j].cv[0].y =  0.;
    v1[j].cv[1].x = v1[j].cv[1].y = 0.;
#if dimension > 2
    v1[j].cv[0].z = 0.;
    v1[j].cv[1].z = 0.;
#endif
    v1[j].b = v1[j].t = v1[j].ft = 0.;
    v1[j].vof1 = p1->vof1, v1[j].vof2 = p1->vof2 ;
    p1++;
  }
  
  foreach_leaf()
    if (T1[] > 0) {
      int j = T1[] - T[];
      coord p = {x,y,z};
      for(scalar a in maybe_close)
	{
	  if(a.i == v1[j].vof1){
	    v1[j].t += dv()*a[];
	    foreach_dimension()
	      v1[j].cv[0].x += dv()*a[]*p.x;
	  }
	  else if(a.i == v1[j].vof2){
	    v1[j].ft += dv()*a[];
	    foreach_dimension()
	      v1[j].cv[1].x += dv()*a[]*p.x;
	  }
	}
    }
 
/** Here we iterate through each pair of bubble in the list bub and
check if pair is already in the list clsncdat (list that have the
pairs of bubbles who are being tracked over time). If they are not in
the list, then add them to list named dat. If they are in the list
then update the information of pair of bubbles like time of contact,
center of volume, their VOF tracer index etc. */
  clsncdat * p = bub->p;
  for(int i = 0; i < length; i++)
    {
      clsncdat datc;
      
      if(v1[i].t > 1e-6 && v1[i].ft > 1e-6) {
	
	  datc.cv[0].x = v1[i].cv[0].x/v1[i].t;
	  datc.cv[0].y = v1[i].cv[0].y/v1[i].t;
#if dimension > 2
	  datc.cv[0].z = v1[i].cv[0].z/v1[i].t;
#endif
	  datc.cv[1].x = v1[i].cv[1].x/v1[i].ft;
	  datc.cv[1].y = v1[i].cv[1].y/v1[i].ft;
#if dimension > 2
	  datc.cv[1].z = v1[i].cv[1].z/v1[i].ft;
#endif
      }
      datc.b = (p->b) + 1;
      datc.vof1 = v1[i].vof1;
      datc.vof2 = v1[i].vof2;
      datc.t = t;
      datc.ft = dt;

      int ctt = 0;
      int ii = (v1[i].t < v1[i].ft) ? 0 : 1 ;

      /** Here we check if the bubble already exists at previous time
	  step, if yes then we will update the information
	  corresponding to same pair in the list otherwise we will add
	  to it. To do this we check the difference between centroid
	  of the bubbles at pevious time step and the bubbles in the
	  current list, if this value is less than velocity of that
	  co-ordinate*time interval, then this bubble pair is updated,
	  otherwise added as the new element */
      double cell_dist = 3*L0/(1<<LEVEL);
      clsncdat * a = dat->p;
      for(int j=0; j<(dat->len/sizeof(clsncdat)); j++){
	if((a->t != t) && (((fabs(a->cv[0].x - datc.cv[ii].x) < cell_dist) && \
			    (fabs(a->cv[0].y - datc.cv[ii].y) < cell_dist)) || \
			   ((fabs(a->cv[1].x - datc.cv[ii].x) < cell_dist) && \
			    (fabs(a->cv[1].y - datc.cv[ii].y) < cell_dist)))) {
	    dat_change(a, &datc);
	    ctt = -1;
	    break;
	  }
	a++;
      }
      if(ctt != -1)
	array_append(dat, &datc, sizeof(clsncdat));
      p++;
    }
  
  int lendat = dat->len/sizeof(clsncdat);;
  // fprintf(stderr,"lendat = %d, time = %g\n", lendat, t);

  array_free(bub);

 /** Here we loop through the list named dat and find the elements
     whose film drainage time is greater than prescribed time and add
     them to the new list named datlist. Add all other elements to
     list named dat1. */

  Array * dat1 = array_new();
  Array * datlist = array_new();
  
  clsncdat * a1 = dat->p;
  for(int i = 0; i < lendat; i++)
    {
      if((a1->ft) > drainage_time && (a1->t) == t) 
	array_append(datlist, a1, sizeof(clsncdat));
      else if((a1->t) == t)
	array_append(dat1, a1, sizeof(clsncdat));
      a1++;
    }
  int len = datlist->len/sizeof(clsncdat);
  fprintf(stderr,"len = %d, time = %g\n", len, t);
  int len1= dat1->len/sizeof(clsncdat);
  fprintf(stderr,"len1 = %d, time = %g\n", len1, t);
  
  /** Replace the VOF tracer for each pair of bubbles in the datlist*/
  if(len)
    {
      int nvar = datasize/sizeof(double), adj[len*nvar];
      int rep[len];
      for (int i = 0; i < len*nvar; i++)
	adj[i] = false;
      
      /** Since we are updating *adj*, we cannot use openMP. */
      clsncdat * a2 = datlist->p;
      for (int i = 0; i < len; i++){
	foreach_leaf() // no openMP
	  if (T[] == a2->b + 1)
	        /** We check which tracer neighbors each bubbles pair T[]
		inlcuding their own tracers */
	    foreach_neighbor()
	      for (scalar s in interfaces)
		if ((s.i != a2->vof1 || s.i != a2->vof2) && s[] > EPS1)
		  adj[i*nvar + s.i] = true;
	
	/** We look for a replacement VOF tracer which is not already
	    neighboring the pair of bubbles and its own tracer. */
	
	rep[i] = -1;
	for (scalar s in interfaces)
	  if ((s.i != a2->vof1 && s.i != a2->vof2) && (s.i >= nvar || !adj[i*nvar + s.i])){
	    rep[i] = s.i; break;
	  }

	/** If we didn't find any, we create a new one. */
	
	if (rep[i] == -1) {
	  scalar t = fclone (list_len (interfaces));
	  reset ({t}, 0.);
	  interfaces = list_append (interfaces, t);
	  rep[i] = t.i;
	}

      /** ### Replacing tracers
	  
	  We perform the replacement for each bubble pair */
	  foreach()
	    if(T[] == a2->b + 1)
	      for(scalar c in maybe_close){
		if((c.i == a2->vof1 || c.i == a2->vof2) && c[]>EPS1)
		  {
		    scalar t = {rep[i]};
		    t[] = c[];
		    c[] = 0.;
		  }
		scalar * list = list_copy ({c});
		for (int i = 0; i < len; i++)
		  list = list_add (list, (scalar){rep[i]});
		boundary (list);
		free (list);
	      }
	
	  a2++;
      }
      
      fprintf(stderr,"coalesced at t=%g len = %d\n", t, len);
    }

/** We will free all the arrays and lists
 */
  free (maybe_close);
  array_free(dat);
  array_free(datlist);
  delete({T,T1});
  
/** Return the list named dat1 where the information of pairs of
bubbles that didn't reach the drainage time is stored. This will be
passed through the control coalescence function in the next time
step */

  return dat1;
}

/** Calling the cntrl_nclsnc function after updating VOF. Here at the
first time step NULL pointer is passed.

##function calling
*/

Array * dat = NULL;
event tracer_advection(i++)
{
  if(!i)
    dat = array_new();
  else
  dat = cntrl_nclsnc(dat, drainage_time);
}

event data_free(i=end)
{ 
  array_free(dat);
}

/*Restore from dump file*/

/*
event defaults (i=0)
{
  if(restore)
    {


    }
}
*/
