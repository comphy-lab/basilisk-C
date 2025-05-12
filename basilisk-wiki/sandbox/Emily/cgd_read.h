#include <stdio.h>
#include <string.h>

struct cgd_read {
  scalar * var;
  FILE *fid; 
  unsigned int it;
  double *dt;
};
// Reads in an opened .cgd file that is specified by p.fid (renamed fid in the function) and interpolates it onto p.var
void basilisk_cgd_read ( struct cgd_read p) {
  scalar var_out = p.var[0];
  unsigned int dimnum,tdimlen,xdimlen,ydimlen;
  char cc1[2], cc2[2], cc3[2];
  dimnum=0;
  fscanf(p.fid,"%u",&dimnum);
  fprintf(stdout,"%% dimnum = %d\n",dimnum);
  // If 3d reads in timeslice 'it'
  if (dimnum == 3){
    fscanf(p.fid,"%s %s %s",cc1,cc2,cc3);
    fprintf(stdout,"%% Dimensions are: %s %s %s\n",cc1,cc2,cc3);
    if (!strcmp(cc1,"t") && !strcmp(cc2,"y") && !strcmp(cc3,"x")) {
      fscanf(p.fid,"%u %u %u",&tdimlen,&ydimlen,&xdimlen);
    }
    else if  (!strcmp(cc1,"t") && !strcmp(cc2,"x") && !strcmp(cc3,"y")) {
      fscanf(p.fid,"%u %u %u",&tdimlen,&xdimlen,&ydimlen);
    }
    double t0[tdimlen];
    for (unsigned int it=0 ; it<tdimlen ; it++)
      fscanf(p.fid, "%lf", &t0[it]);
    if (p.it >= tdimlen-1)
      *p.dt=1.e37;
    else
      *p.dt=t0[p.it+1]-t0[p.it];
  }
  else if(dimnum == 2) {
    fscanf(p.fid,"%s %s",cc1,cc2);
    fprintf(stdout,"%% Dimensions are %s %s\n",cc1,cc2);
    if (!strcmp(cc1,"y") && !strcmp(cc2,"x")) {
      fscanf(p.fid,"%u %u",&ydimlen,&xdimlen);
      fprintf(stdout,"%% %s dim length %d, %s dim length %d\n",cc1,ydimlen,cc2,xdimlen);
    }
    else if  (!strcmp(cc1,"x") && !strcmp(cc2,"y")) {
      fscanf(p.fid,"%u %u",&xdimlen,&ydimlen);
      fprintf(stdout,"%% %s dim length %d, %s dim length %d\n",cc1,xdimlen,cc2,ydimlen);
    }
  }
  double y0[ydimlen],x0[xdimlen];
  double var[ydimlen][xdimlen];
  if (!strcmp(cc1,"t")) {
    if (!strcmp(cc2,"y") && !strcmp(cc3,"x")) {
      for (unsigned int iy=0 ; iy<ydimlen ; iy++)
	fscanf(p.fid, "%lf", &y0[iy]);
      for (unsigned int ix=0 ; ix<xdimlen ; ix++)
	fscanf(p.fid, "%lf", &x0[ix]);
      for (unsigned int i=0 ; i<= p.it ; i++) {
	for (unsigned int iy=0 ; iy<ydimlen && !feof(p.fid); iy++){
	  for (unsigned int ix=0 ; ix<xdimlen && !feof(p.fid); ix++){
	    fscanf(p.fid, "%lf", &var[iy][ix] );
	  }
	}
      }
    }
    else if  (!strcmp(cc2,"x") && !strcmp(cc3,"y")) {
      for (unsigned int ix=0 ; ix<xdimlen ; ix++)
	fscanf(p.fid, "%lf", &x0[ix]);
      for (unsigned int iy=0 ; iy<ydimlen ; iy++)
	fscanf(p.fid, "%lf", &y0[iy]);
      for (unsigned int i=0 ; i<= p.it ; i++) {
	for (unsigned int ix=0 ; ix<xdimlen && !feof(p.fid); ix++){
	  for (unsigned int iy=0 ; iy<ydimlen && !feof(p.fid); iy++){
	    fscanf(p.fid, "%lf", &var[iy][ix] );
	  }
	}
      }
    }
  }
  else if (!strcmp(cc1,"y") && !strcmp(cc2,"x")) {
    for (unsigned int iy=0 ; iy<ydimlen ; iy++)
      fscanf(p.fid, "%lf", &y0[iy]);
    for (unsigned int ix=0 ; ix<xdimlen ; ix++)
      fscanf(p.fid, "%lf", &x0[ix]);
    for (unsigned int iy=0 ; iy<ydimlen && !feof(p.fid); iy++){
      for (unsigned int ix=0 ; ix<xdimlen && !feof(p.fid); ix++){
	fscanf(p.fid, "%lf", &var[iy][ix] );
      }
    }
  }
  else if  (!strcmp(cc1,"x") && !strcmp(cc2,"y")) {
    for (unsigned int ix=0 ; ix<xdimlen ; ix++)
      fscanf(p.fid, "%lf", &x0[ix]);
    for (unsigned int iy=0 ; iy<ydimlen ; iy++)
      fscanf(p.fid, "%lf", &y0[iy]);
    for (unsigned int ix=0 ; ix<xdimlen && !feof(p.fid); ix++){
      for (unsigned int iy=0 ; iy<ydimlen && !feof(p.fid); iy++){
	fscanf(p.fid, "%lf", &var[iy][ix] );
      }
    }	
  }
  // Below is the code to interpolate from the input file onto the variable
  double dx0 = x0[1] - x0[0];
  double dy0 = y0[1] - y0[0];
  foreach()  {
      double out;
      int ix=(int) floor( ( x - x0[0] ) / dx0 );
      int iy=(int) floor( ( y - y0[0] ) / dy0 );
      double dx=( x - x0[0] ) / dx0 - floor( ( x - x0[0] ) / dx0 );
      double dy=( y - y0[0] ) / dy0 - floor( ( y - y0[0] ) / dy0 );
      if ( ( ix >= 0 && ix < xdimlen-1 && iy >= 0 && iy < ydimlen-1 ) )	{
	  out=var[iy][ix] * ( 1. - dx ) * ( 1. - dy ) +
	    var[iy][ix+1] * dx * ( 1. - dy ) + 
	    var[iy+1][ix] * ( 1. - dx ) * dy + 
	    var[iy+1][ix+1] * dx * dy;
	}
      else 
	{ 
	  out = 0.;
	}
      var_out[] = out;
  };
}

struct Deformation_cgd_read {
  scalar d;
  double x, y;
  unsigned int xdimlen,ydimlen;
  double *x0,*y0,*var;
  FILE *fid;
  int (* iterate) (void);
};

void cgd_interpolate (struct Deformation_cgd_read p) {
  double dx0 = *(p.x0+1) - *(p.x0);
  double dy0 = *(p.y0+1) - *(p.y0);
  double outval;
  /*for (unsigned int ix=0 ; ix<p.xdimlen; ix++){
    for (unsigned int iy=0 ; iy<p.ydimlen; iy++){
    }
    }*/
  //fprintf(stdout,"%% Step values are dx0 = %10.5f dy0 = %10.5f\n",dx0,dy0);
  //fprintf(stdout,"%% Lowest values are x0 = %10.5f y0 = %10.5f\n",*(p.x0),*(p.y0));
  //fprintf(stdout,"%% Centering x = %10.5f y = %10.5f\n",p.x,p.y);
  foreach()  {
    int ix=(int) floor( ( x - p.x - *(p.x0) ) / dx0 );
    int iy=(int) floor( ( y - p.y - *(p.y0) ) / dy0 );
    double dx=( x - p.x - *(p.x0) ) / dx0 - (double)ix;
    double dy=( y - p.y - *(p.y0) ) / dy0 - (double)iy;
    if ( ( ix >= 0 && ix < p.xdimlen-1 && iy >= 0 && iy < p.ydimlen-1 ) ) {
      outval=*(p.var+ix*p.ydimlen+iy) * ( 1. - dx ) * ( 1. - dy ) +
	*(p.var+(ix+1)*p.ydimlen+iy) * dx * ( 1. - dy ) +
	*(p.var+ix*p.ydimlen+iy+1) * ( 1. - dx ) * dy +
	*(p.var+(ix+1)*p.ydimlen+iy+1) * dx * dy;
    }
    else {
      outval = 0.;
    }
    val(p.d,0,0)=outval;
  }
}
// Reads in an opened .cgd file that is specified by p.fid (renamed fid in the function) and interpolates it onto p.var
void deformation_cgd_read (struct Deformation_cgd_read p) {
  unsigned int dimnum;
  char cc1[2], cc2[2];
  dimnum=0;
  fscanf(p.fid,"%u",&dimnum);
  fprintf(stdout,"%% dimnum = %d\n",dimnum);
  fscanf(p.fid,"%s %s",cc1,cc2);
  fprintf(stdout,"%% Dimensions are %s %s\n",cc1,cc2);
  if (!strcmp(cc1,"y") && !strcmp(cc2,"x")) {
    fscanf(p.fid,"%u %u",&p.ydimlen,&p.xdimlen);
    fprintf(stdout,"%% %s dim length %d, %s dim length %d\n",cc1,p.ydimlen,cc2,p.xdimlen);
  }
  else if  (!strcmp(cc1,"x") && !strcmp(cc2,"y")) {
    fscanf(p.fid,"%u %u",&p.xdimlen,&p.ydimlen);
    fprintf(stdout,"%% %s dim length %d, %s dim length %d\n",cc1,p.xdimlen,cc2,p.ydimlen);
  }
  p.x0=(double*)malloc(p.xdimlen*sizeof(double));
  p.y0=(double*)malloc(p.ydimlen*sizeof(double));
  p.var=(double*)malloc((p.ydimlen*p.xdimlen)*sizeof(double));;
  if (!strcmp(cc1,"y") && !strcmp(cc2,"x")) {
    for (unsigned int iy=0 ; iy<p.ydimlen ; iy++)
      fscanf(p.fid, "%lf", p.y0+iy);
    for (unsigned int ix=0 ; ix<p.xdimlen ; ix++)
      fscanf(p.fid, "%lf", p.x0+ix);
    for (unsigned int iy=0 ; iy<p.ydimlen && !feof(p.fid); iy++){
      for (unsigned int ix=0 ; ix<p.xdimlen && !feof(p.fid); ix++){
	fscanf(p.fid, "%lf", p.var+ix*p.ydimlen+iy );
      }
    }
  }
  else if  (!strcmp(cc1,"x") && !strcmp(cc2,"y")) {
    for (unsigned int ix=0 ; ix<p.xdimlen ; ix++){
      fscanf(p.fid, "%lf", p.x0+ix);
    }
    for (unsigned int iy=0 ; iy<p.ydimlen ; iy++) {
      fscanf(p.fid, "%lf", p.y0+iy);
    }
    for (unsigned int ix=0 ; ix<p.xdimlen; ix++){
      for (unsigned int iy=0 ; iy<p.ydimlen; iy++){
	fscanf(p.fid, "%lf", p.var+ix*p.ydimlen+iy);
      }
    }
  }
  fprintf(stdout,"%% Done reading in\n");

  /**
     This is taken from the okada fault routine.  Hopefully it iterates this step until the deformation is well resolved */
  scalar hold[],def[];
  // save the initial water depth
  scalar_clone (hold, h); // clone h into hold
  foreach(){
    hold[] = h[];
    def[]=0.;
  }
  boundary ({hold});
  //p.d = def;
  p.d = h;
  int nitermax = 20;
  do {
    scalar l[];
    // Below is the code to interpolate from the input file onto the variable
    cgd_interpolate (p);
    foreach(){
      l[]=level;
    }
      printf ("%% file: A-%d\n", 20-nitermax);
      output_field ({eta, h, zb, def, l}, stdout, n = 1 << 8, linear = false);
    foreach(){
      //h[]=def[];
      h[] = hold[] > dry ? max (0., hold[] + h[]) : hold[];
      eta[] = zb[] + h[];
    }
    foreach(){
      l[]=level;
    }
      printf ("%% file: B-%d\n", 20-nitermax);
      output_field ({eta, h, zb, def, l}, stdout, n = 1 << 8, linear = false);
  } while (p.iterate && p.iterate() && nitermax--);
}
