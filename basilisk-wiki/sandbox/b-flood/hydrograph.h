/**
## Hydrograph

This function return the inflow to impose by linearly interpolating the data. */

double hydrograph (const char * name, double mult)
{
  FILE * fp;
  if ((fp = fopen (name , "r")) == NULL) {
    fprintf (stderr,"cannot open hydrograph data file.\n name = %s ",name);
    assert (false);
  }
  double time = 0, timea, q = 0, qa, alpha = 0;
  //Read the data at each call => Can definitively be optimised !
  do {
    qa = q;
    timea = time;
    if (fscanf (fp, "%lf \t %lf \n", &time, &q) == EOF) break;
    time *= mult;
    if (time - timea != 0) alpha = (q - qa)/(time - timea);
    else alpha = 0;
  } while (time < t);
  fclose (fp);
  /*
    If the solver time is sup to the larger time of the data, return the last value   */
  if (timea  >= t)  return q;
  // Else, return the interpolation
  else return alpha*(t - timea) + qa;
}


/**
## Line structure
*/

struct linedis{
  double xi,yi,xf,yf; // Starting and end point of the line
};

/**
## Computing the inflow passing through a line
*/
double compdischarge(struct linedis l){
#if IMPLICIT
  vector u[];
  foreach(){
    if( h[] > dry ){
      foreach_dimension()
	u.x[] = q.x[]/h[];
    }
  }
#endif
  double ds = L0/N;
  double dis = 0;
  // Trivial in one dimension
#if dimension == 1
    return interpolate(u.x,l.xi,0)*interpolate(h,l.xi,0);
#endif
    
  if(!l.yf){ // Horizontal line
    for(double x = l.xi; x <=l.xf; x+= ds) 
      if (interpolate(h,x,l.yi) > dry ) dis += interpolate(u.y,x,l.yi)*interpolate(h,x,l.yi)*ds;    
    return fabs(dis);
  }

  if(!l.xf) {// Vertical line
    for(double y = l.yi; y <=l.yf; y+= ds) 
      if (interpolate(h,l.xi,y) > dry )  dis += interpolate(u.x,l.xi,y)*interpolate(h,l.xi,y)*ds;
    return fabs(dis);
  }

  return 0; // error
}

/**
## Link to the homepage
* [Homepage](http://basilisk.fr/sandbox/B-flood/Readme)
*/
