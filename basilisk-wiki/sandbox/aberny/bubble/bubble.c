/**  
This small code is just an example for the library bubbleShape.h and
findBond.h. We are solving the Young-Laplace equation in a static case.  

Please note that if you are using the library "findBond.h", then you don't need
to include "bubbleShape.h".*/


/**
In this example, we want to output all the information about the spherical cap. For that, we define outCap at 1*/
#define outCap 1

#define WRONG 0

double velocityNoise(double amplitude) {
  return( (2*(rand()/(double)RAND_MAX)-1)*amplitude );
}

#define InterfaceNoise 0

#include "bubbleShape.h"
#include "findBond.h"

#define dimension 2

/**
Computation of the surface of the interface

We give an axi-symmetric segment, and it returns the surface of the conical
shape if it's 3D.*/

double surfaceCyl(coord p1, coord p2) {
  double s = 0;
  if(fabs(p2.y-p1.y)>0.){
    double tanPhi = (p2.x-p1.x)/(p2.y-p1.y);
    s = M_PI*fabs(sq(p2.y)-sq(p1.y))/cos(atan(tanPhi));
  }
  else
    s =2*M_PI*fabs(p2.x-p1.x)*p1.y;
  return s;
}

double surfaceC (coord * bubble) {
  double s = 0;
  int i = 0;

  /**
  We are in cylindrical coordinate.

  x stands for the heihgt z
  y stands for the radius r

  We want to compute the surface of the bubble shape, for a radius below 10.
  */
  while (bubble[i].x != nodata && bubble[i].y<10.){
    s += surfaceCyl(bubble[i], bubble[i+1]);
    // s += M_PI*fabs(bubble[i].y+bubble[i-1].y)*fabs(bubble[i].x-bubble[i-1].x);
    // s += M_PI*fabs(sq(bubble[i].y)-sq(bubble[i-1].y));
    i+=2;
  }
  return fabs(s);
}


int main(int argc, char const *argv[])
{
  double Bond = 0.2;
  int caseBond = 1;
  if (argc >= 2){
      Bond = atof (argv[1]);
  }
  if (argc >= 3) {
    caseBond = atoi (argv[2]);
  }

  /**
  The function shape will return the data in a coord structure type.*/

  coord* data;

  /**
  We also want to compute 2 other Bond number. One is an elliptical
  approximation. The second is based on the total volume of the bubble. The 
  reason is that, when working on experimental photography of buble, we can 
  very easily measure the elliptical Bond number. We call this Bond number 
  BondGhabache.*/

  double* BondGhabache = NULL;
  double BoGha = 0;

  double* BoEff = NULL;
  double Bo = 0;

  BondGhabache = &BoGha;
  BoEff = &Bo;

  /** 
  If we want to use the Effective Bond as input data instead of the
  normal one, we will call the function findBond, with the argument caseBond
  equal to one. If we want to find the input giving the BondGhabache, we will
  call the same function, findBond, with case Bond equal to 2. */
  double BondFunction = Bond;

  if(caseBond<3){  
    Bond = findBond(BondFunction,caseBond);
  }

  /**    
  The function shape take 3 input parameter. The Bond number based on the
  initial curvature, the Bond number used by Ghabache et al (unknown at this
  point), and the Bond number, based on the volume of the drop (unknown at this
  point). The last 2 parameters will be compute by the shape function.*/

  int* size = NULL;
  int taille = 0;
  size = &taille;

  Circle* hollow = NULL;
  Circle fillet;
  hollow = &fillet;

  Circle* cap = NULL;
  Circle topCap;
  cap = &topCap;

  double* phi = NULL;
  double phiCirle = 0;
  phi = &phiCirle;

  // double* radius = NULL;
  // double rMax = 0;
  // radius = &rMax;

  int type = 0;

  data = shape(Bond, BondGhabache, BoEff, phi, size, 00, hollow, cap, type, 0);
  // data = bubbleTot (Bond, 0.05, size);
  // coord * doubleCap;
  // doubleCap = partialCircle(topCap, 0.01, phiCirle);

  coord * capTop;
  capTop = topCircle(topCap, phiCirle);

  // coord * bubbleB;
  // int * size;
  // bubbleB = bubbleTot(BondFunction, 0.01, size);



  #if WRONG
  fprintf(stdout, "Carefull: Bond number based on surface,\n not on volume as it should be\n"); 
  #endif

  fprintf(stdout, "Bond en Input = %f\n", Bond);
  fprintf(stdout, "BondGhabache = %f\n", *BondGhabache);
  fprintf(stdout, "BondEffective = %f\n", *BoEff);
  // fprintf(stdout, "Circle: x_c= %f, y_c= %f, r= %f\n",
  //          fillet.x, fillet.y, fillet.r);
  fprintf(stdout, "Cap: x_cap = %f, y_cap = %f, r-cap = %f\n", topCap.x, topCap.y, topCap.r);
  // fprintf(stdout, "radius of bubble: %f\n", rMax);
  // fflush(stdout);
  int i = 0;

  while(data[i].x != nodata){
  // for (int i = 0; i<taille-1; i++){
    fprintf(stderr, "%f %f %i\n",data[i].y, data[i].x, i);
    i++;
    if (i%2==0)
      fprintf(stderr, "\n");
  }
  

  // FILE * fpCap = fopen("doubleCap", "w");
  // i = 0;
  // while(doubleCap[i].x != nodata) {
  //   fprintf(fpCap, "%f %f\n",doubleCap[i].x, doubleCap[i].y);
  //   if (i % 2 ==1)
  //     fprintf(fpCap, "\n");
  //   i++;
  // }

  FILE * fpCap = fopen("cap", "w");
  i = 0;
  while(capTop[i].x != nodata) {
    fprintf(fpCap, "%f %f %d\n",capTop[i].x, capTop[i].y,i);
    if (i % 2 ==1)
      fprintf(fpCap, "\n");
    i++;
  }

if (type == -1){
  double sCap = surfaceC(capTop);
  double sTot = surfaceC(data);
  fprintf(stdout, "surface bubble %f\n", sTot - sCap);
  fprintf(stdout, "tot %f, cap %f\n", sTot, sCap);
}

if (type == 0){
  double sTot = surfaceC(data);
  double sOrigine = M_PI*sq(10);
  fprintf(stdout, "surface perturb %f\n", sTot);
  fprintf(stdout, "surface origine %f\n", sOrigine);
  fprintf(stdout, "surface diff %f\n",sTot - sOrigine);
}
/** Here, "size" correspond to the number of points we have compute. In the 2D
array, we have twice more points, minus 2. Indeed, the first point is not
repeated, so is the last one.*/

#if dimension > 2

  static FILE * fp = fopen("shape3D", "w");

  coord* shapeVolume = shape3D(data, *size);

  i = 0;

  // for (int i = 0; i< ((* size)*360); i++) {
  while (shapeVolume[i].x != nodata) {
    if (i % 3 == 0) 
      fprintf(fp, "\n");
    // fprintf(stdout, "%d %g\n",i, theta);
    // if (fabs(shapeVolume[i].y) <7 && fabs(shapeVolume[i].z) <7 && fabs(shapeVolume[i].x)< 7)
    fprintf(fp, "%g %g %g\n", shapeVolume[i].x, shapeVolume[i].y, shapeVolume[i].z);
    i++;
  }
#endif

  return 0;
}


/**
Once the code is over, we can plot the shape of the bubble, without the top
spherical cap

~~~gnuplot Shape of the bubble
set term svg size 1920,1080
set size ratio -1
set xrange [-5:5]
set yrange [:1]
plot 'log' u (-$1):2 w l t "" lt rgb "red",\
     'log' u 1:2 w l t "shape of the bubble" lt rgb "red",\
     'cap' u (-$1):2 w l t "" lt rgb "red",\
     'cap' u 1:2 w l t "" lt rgb "red"
~~~ 

*/
