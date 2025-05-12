#include <complex.h>
#include "view.h"
#define PI 3.141592653589793

struct Polygon {
  int nv;
  coord * Vertices;
};

int ReadPolygon ( char * polyfile, struct Polygon * _Poly)
{
  char buffer[200]; 
  int nv;
  FILE * fp;
  double xv,yv;

  if(!(fp = fopen (polyfile, "r")))
  {
    printf("Failed to open polygon file!\n");
    return -1;
  }

  _Poly->Vertices   = (coord *) malloc(0*sizeof(coord));
  nv=0;

  while ( fgets(buffer, 200, fp)!=NULL )
  {	  	
    sscanf(buffer,"%lf %lf",&xv,&yv);
    _Poly->Vertices = (coord *) realloc( _Poly->Vertices , (nv+1)*sizeof(coord));
    (_Poly->Vertices)[nv].x = xv;
    (_Poly->Vertices)[nv].y = yv;
    nv++;    
  }

  _Poly->nv = nv;

  return(0);
}

bool isInPolygon(coord P, struct Polygon * _Poly )
{
   int i;
   coord Pcurr;
   double complex zprev, zcurr;
   double sum,dsum;

   Pcurr = (_Poly->Vertices)[0];
   zcurr = (Pcurr.x-P.x) + I*(Pcurr.y-P.y);

   if(cabs(zcurr)<1e-6) // P almost coincides with current vertex, leave
     return(true);

   sum = 0.;

   for(i=1;i<(_Poly->nv);i++)
   {
        zprev = zcurr;
        Pcurr = (_Poly->Vertices)[i];
        zcurr = (Pcurr.x-P.x) + I*(Pcurr.y-P.y);

        if(cabs(zcurr)<1e-6) // P almost coincides with current vertex, leave
        return(true);

        dsum = cimag(clog(zcurr)-clog(zprev));
        dsum += fabs(dsum)>PI ? -2.0*PI*dsum/fabs(dsum) : 0.; 
        sum += dsum; 
   }

   // Ensure polygon is closed

   i = (_Poly->nv)-1;
   zprev = zcurr;
   Pcurr = (_Poly->Vertices)[i];
   zcurr = (Pcurr.x-P.x) + I*(Pcurr.y-P.y);

   dsum = cimag(clog(zcurr)-clog(zprev));
   dsum += fabs(dsum)>PI ? -2.0*PI*dsum/fabs(dsum) : 0.; 
   sum += dsum; 

   return( fabs(sum)>PI ); 
}

int draw_polygon ( struct Polygon * _Poly , float lc[3] = {0}, float lw = 1.)
{
  bview * view = draw();
  draw_lines (view, lc, lw) {
  glBegin (GL_LINE_LOOP);
  for(int k=1;k<(_Poly->nv);k++)
    glvertex2d (view, (_Poly->Vertices)[k].x, (_Poly->Vertices)[k].y);
  glEnd();
  view->ni++;
  }
  return 0;
}

int draw_polygon_on_surface ( struct Polygon * _Poly , scalar s, float lc[3] = {0}, float lw = 1., float z_offset = 0.)
{
  bview * view = draw();
  draw_lines (view, lc, lw) {
  glBegin (GL_LINE_LOOP);
  for(int k=1;k<(_Poly->nv);k++)
  {
    Point point = locate ( (_Poly->Vertices)[k].x, (_Poly->Vertices)[k].y );
    glvertex3d (view, x, y , s[] + z_offset);
  }
  glEnd();
  view->ni++;
  }
  return 0;
}