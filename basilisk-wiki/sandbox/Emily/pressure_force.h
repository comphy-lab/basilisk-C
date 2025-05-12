#include <stdio.h>
#include <string.h>
/** Two dimenional pressure forcing for atmospheric forcing. Inputs a 2D .cgd file of the pressure anomaly with time and radius as the two dimensions. The function panom_calc calculates the pressure anomaly at time t and point x,y assuming that the centre of the radial anomaly is at x0,y0
*/

double tp0,dtp,rp0,drp;
double* panom0;
size_t tdimlen,rdimlen;

double panom_calc( double t, double x, double y , double x0, double y0 );

void read_panom (FILE *fid){
  size_t dimnum;
  char cc1[2], cc2[2];
if (fid == NULL){
    fprintf(stderr,"Error: CGD file not found");
    exit(1);
  }
  fprintf(stderr, "*** Pressure Anomaly CGD Read ***\n");
  dimnum=0;
  fscanf(fid,"%zu",&dimnum);
  fprintf(stdout,"%% dimnum = %zu\n",dimnum);
  if(dimnum == 2) {
    fscanf(fid,"%s %s",cc1,cc2);
    fprintf(stdout,"%% Dimensions are %s %s\n",cc1,cc2);
    if (!strcmp(cc1,"t") && !strcmp(cc2,"r")) {
      fscanf(fid,"%zu %zu",&tdimlen,&rdimlen);
      fprintf(stdout,"%% %s dim length %zu, %s dim length %zu\n",cc1,tdimlen,cc2,rdimlen);
    }
    else {
      fscanf(fid,"%zu %zu",&rdimlen,&tdimlen);
      fprintf(stdout,"%% %s dim length %zu, %s dim length %zu\n",cc1,tdimlen,cc2,rdimlen);
    }
  }
  fprintf(stdout,"%s %s\n",cc1,cc2);
  fflush(stdout);
  double* t0cgd=(double*)malloc(tdimlen*sizeof(double));
  double* r0cgd=(double*)malloc(rdimlen*sizeof(double));
  panom0=(double*)malloc((rdimlen*tdimlen)*sizeof(double));
  fprintf(stdout,"%s %s\n",cc1,cc2);

 if (!strcmp(cc1,"t") && !strcmp(cc2,"r")) {
   fprintf(stdout,"%s %s\n",cc1,cc2);
    for (size_t it=0 ; it<tdimlen ; it++)
     fscanf(fid, "%lf", t0cgd+it);
    for (size_t ir=0 ; ir<rdimlen ; ir++)
      fscanf(fid, "%lf", r0cgd+ir);
    for (size_t it=0 ; it<tdimlen && !feof(fid); it++){
      for (size_t ir=0 ; ir<rdimlen && !feof(fid); ir++){
        fscanf(fid, "%lf", panom0+it*rdimlen+ir );
      }
    }
  }
 dtp = *(t0cgd+1) - *t0cgd;
 drp = *(r0cgd+1) - *r0cgd;
 tp0 = *t0cgd;
 rp0 = *r0cgd;
 fclose(fid);
}

double panom_calc( double t, double x, double y , double x0, double y0 ) {
  double r,dr,dt,panom_val;
  int ir,it;

  r = sqrt( ( x - x0 ) * ( x - x0 ) + ( y - y0 ) * ( y - y0 ) );
  //fprintf(stdout,"panom_calc %7.0f %7.0f %6.4f  %7.0f\n",x,y,t,r);
  it = (int) floor( ( t - tp0 ) / dtp );
  dt = ( t - tp0 ) / dtp - floor( ( t - tp0 ) / dtp );
  ir = (int) floor( ( r - rp0 ) / drp );
  dr = ( r - rp0 ) / drp - floor( ( r - rp0 ) / drp );
  if ( ( ir >= 0 && ir < rdimlen-1 && it >= 0 && it < tdimlen-1 ) )     {
        panom_val = *(panom0 + it * rdimlen + ir) * ( 1. - dr ) * ( 1. - dt ) +
          *(panom0 + it * rdimlen + ir + 1) * dr * ( 1. - dt ) +
          *(panom0 + (it + 1) * rdimlen + ir) * ( 1. - dr ) * dt +
          *(panom0 + (it + 1) * rdimlen + ir + 1) * dr * dt;
        }
      else
        {
          panom_val = 0.;
        }

  return panom_val;
}
