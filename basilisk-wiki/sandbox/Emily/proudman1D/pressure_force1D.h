#include <stdio.h>
#include <string.h>

double tp0,dtp,rp0,drp;
double* panom0;
size_t tdimlen,rdimlen;

double panom_calc( double t, double x, double x0 );

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
      for (size_t ir=0 ; ir= 0 && ir < rdimlen-1 && it >= 0 && it < tdimlen-1 ) )     {
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
