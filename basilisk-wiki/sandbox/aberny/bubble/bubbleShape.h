/**
#Way of solving the equation

The goal of this library is to solve the static interface of an air bubble 
close to a liquid-air interface. This is solved in 3 parts. First, we solve 
the water-bubble interface from x=0 up to dzdx=1. Then we switch the axis to 
continue to solve the water-bubble interface up to the end. We change the axis 
because for a given x value, we can have 2 z values. We are not solving this 
problem with x as a function of z because it would mean to start with 
infinite derivative. Finally, we solve the tail of this function. To solve the 
tail, we will iterate until the starting point of the tail, between 2 steps 
goes belong a certain value.


This library will solve the Young-Laplace equations. The equation link to the
bottom part of the bubble is the following one:

$$
\frac{1}{\kappa_1}+\frac{1}{\kappa_2} = Bo\times y + 2R_0 
$$

Where $\kappa_i$ is a principal radius of curvature of the shape. $Bo$ is the 
Bond number. 

$Bo = \frac{\Delta\rho L^2 g}{\gamma}$

Here, $L$ is a characteristic length scale, $g$ is the gravity acceleration,
$\gamma$ is the surface tension between the liquid and the air. We choose, as
characteristic length scale, the curvature of the bubble at the bottom.

The Young Laplace equation is solved on an axisymmetric bubble shape. The 
differential equation we have to solve is then:

$$
\frac{\frac{d^2y}{dx^2}}{\left(1+\left(\frac{dy}{dx}\right)^2\right)^{3/2}} + \frac{\frac{dy}{dx}}{x\left(1+\left(\frac{dy}{dx}\right)^2\right)^{1/2}} = Bo\times y+2
$$

We won't solve the upper part of the bubble. We still need to solve the tail 
part of the bubble shape. The tail is also obtained by solving a dimensionless 
Young-Laplace's equation. This time, we have to solve the following equation: 

$$
\frac{1}{\kappa_1}+\frac{1}{\kappa_2} = Bo\times (z-L)
$$ 

Here, $L$ is the dimensionless difference between the water level at $x=\infty$ and the lowest part of the bubble.

To find the starting point of the tail, we will select one starting point, 
solve the differential equation, and apply a test on the solution to know if 
our point is below or above the solution.

The method implemented here to solve the equation of the static droplet was 
partially described in 1962 bu H.M. Princen.

For the differential equation solver, we will use a Runge Kutta algorithm, at 
order 4.

#Geneal definition
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
fstep and gstep are the step used in the resolution of the bubble shape. f is 
the first part (until dz/dx reach 1). For g, we reduce the step. Indeed, for 
the last point of the drop, we will reach a high derivative in absolute value. 
In order to have a good precisio, we reduce the step.*/

#define fstep (1./1000.)
#define gstep (fstep/10.)
#define tstep (fstep/5.)

/**
The size of the domain is 10.

However, due to an accoustician request, we can increase it up to 30.*/

#if ACOUSTIC
#define SIZE 30.
#else
#define SIZE 10.
#endif

/**
The case f, g and tail are used to select the good function during the 
computation of the k_i coefficient.*/

#define fCase 1
#define gCase 2
#define tailCase 3

/**
The stopCondition is the precision we want to have for the beginning of the 
tail.*/

#define stopCondition 0.00000000000001

/**
We define a "stupid" value for the post-process of the data*/

#define stupid (-1000)

/**
Sometimes, we want to know where is the spherical cap, and its radius. By default, we don't want this output. However, this can be change by redifining out*/

#ifndef outCap
#define outCap 0
#endif

/**
We define a struct that can store all the information for a circle which are
the coordinate of */

typedef struct Circle Circle;
struct Circle{
  double x;
  double y;
  double r;
};

/**
#Some useful function

The min function is not defined. We define it. It returns the min value 
between x and y*/

double minFunction(double x, double y){
  if (x>y)
    return y;
  else
    return x;
}

/**
Definition of trigonometric function that take angle in degrees as input*/

double cosd(double angle) {
  return cos(angle * M_PI/180.);
}

double sind(double angle) {
  return sin(angle * M_PI/180.);
}
/**
The find function take 3 arguments: the number we want to find, an array and 
the value we are looking for. It will return the position of the closest 
value, by default.

If all the value of the array are above the value we are looking for, we 
return the error code 200.*/

int find(double ToFind, double *array, int arraySize) {
  double found = -HUGE;
  int step = -1;
  double val = array[0];
  int i = 0;
  while(val != stupid){
    if (val < ToFind) {
      if (found<val) {
        found = val;
        step = i;
      }
    }
    i++;
    val = array[i];
  }
  if (step<0)
    exit(200);
  return step;
}

/**
#Definition of the functions for the RK-4 algorithm

The function f is the evolution of the curvature, considering z as a function 
of x.*/

double fz(double x, double z, double dz, 
  double Bo) {
  double term1 = (-dz*(1. + sq(dz)))/x;
  double term2 = (Bo*z+2.)*pow( (1.+sq(dz)) , (1.5) );
  return term1+term2;
}

/**
The function g is the evolution of the curvature, considering x as a function 
of z.*/

double gx(double z, double x, double dx, 
  double Bo) {
  double term1 = (1.+sq(dx))/x;
  double term2 = -(Bo*z+2.)*pow( (1.+sq(dx)) , (1.5) );
  return term1+term2;
}

/**
The function tail is the differential equation rulling the tail interface.*/

double tail(double x, double z, double dz, 
  double Bo, double L) {
  double term1 = (Bo*(z-L)*pow( (1.+sq(dz)), (1.5) ));
  double term2 = -dz/x*(1.+sq(dz));
  return term1+term2;
}

/**
Since we want to use a Runge Kutta algorithm, we have to compute the 
coefficient k_i, for the different function, f, g and the tail.

1 corresponds to the first part: shape of the bubble, z as a function of x.
2 corresponds to the second part: shape of the bubble, x as a function of z.
3 corresponds to the tail part.

The integer fgt correspond to the selector between the 3 different case. If 
we didn't implement 1, 2 or 3, the code will return an error code (1 to 4) 
corresponding to the coefficient k_i.*/

double k1(double tn, double yn, double dyn, 
  double Bo, int fgt, double L) {
  if (fgt == 1)
    return fz(tn, yn, dyn, Bo);
  else if (fgt == 2)
    return gx(tn, yn, dyn, Bo);
  else if (fgt == 3)
    return tail(tn, yn, dyn, Bo, L);
  else
    exit(1);
}

double k2(double tn, double yn, double dyn, 
  double Bo, double k1, int fgt, double L) {
  if (fgt == 1)
    return fz( (tn+fstep/2.), (yn+fstep/2.*dyn), (dyn+fstep/2.*k1), Bo);
  else if (fgt == 2)
    return gx( (tn+gstep/2.), (yn+gstep/2.*dyn), (dyn+gstep/2.*k1), Bo);
  else if (fgt == 3)
    return tail( (tn+tstep/2.), (yn+tstep/2.*dyn), (dyn+tstep/2.*k1), Bo, L);
  else
    exit(2);
}

double k3(double tn, double yn, double dyn, 
  double Bo, double k1, double k2, int fgt, double L) {
  if (fgt == 1)
    return fz( (tn+fstep/2.), 
      (yn+fstep/2.*dyn+sq(fstep)/2.*k2), 
      (dyn+fstep/2.*k2), Bo);
  else if (fgt == 2)
    return gx( (tn+gstep/2.), 
      (yn+gstep/2.*dyn+sq(gstep)/2.*k2), 
      (dyn+gstep/2.*k2), Bo);
  else if (fgt == 3)
    return tail( (tn+tstep/2.), 
      (yn+tstep/2.*dyn+sq(tstep)/2.*k2), 
      (dyn+tstep/2.*k2), Bo, L);
  else
    exit(3);
}

double k4(double tn, double yn, double dyn, 
  double Bo, double k2, double k3, int fgt, double L) {
  if (fgt == 1)
    return fz( (tn+fstep), 
      (yn+fstep*dyn+sq(fstep)/2.*k2), 
      (dyn+fstep*k3), Bo);
  else if (fgt == 2)
    return gx( (tn+gstep), 
      (yn+gstep*dyn+sq(gstep)/2.*k2), 
      (dyn+gstep*k3), Bo);
  else if (fgt == 3)
    return tail( (tn+tstep), 
      (yn+tstep*dyn+sq(tstep)/2.*k2), 
      (dyn+tstep*k3), Bo, L);
  else 
    exit(4);
}

/** 
After the computation of the k_i coefficients, we can now obtain the value
y_( n+1) and dy_(n+1).

Again, we use the integer fgt to select the case 1, 2 or 3. Again, if we didn't
select the good case, the program stop with the error code 100 when the problem
is in the function yn1 and 101 when the problem is in the function dyn1*/

double yn1(double yn, double dyn, double k_1, double k_2, double k_3,
 int fgt) {
  if (fgt == 1)
    return yn+fstep*dyn+(sq(fstep)/6.)*(k_1+k_2+k_3);
  else if (fgt == 2)
    return yn+gstep*dyn+(sq(gstep)/6.)*(k_1+k_2+k_3);
  else if (fgt == 3)
    return yn+tstep*dyn+(sq(tstep)/6.)*(k_1+k_2+k_3);
  else
    exit(100);
}

double dyn1(double dyn, double k_1, double k_2, double k_3, double k_4,
 int fgt) {
  if (fgt == 1)
    return dyn+(fstep/6.)*(k_1+2.*k_2+2.*k_3+k_4);
  else if (fgt == 2)
    return dyn+(gstep/6.)*(k_1+2.*k_2+2.*k_3+k_4);
  else if (fgt == 3)
    return dyn+(tstep/6.)*(k_1+2.*k_2+2.*k_3+k_4);
  else
    exit(101);
}

/** 
This function compute the surface define by a segment and the origin axis.
*/

double surface(double xn, double xn1, double yn, double yn1){
  return (1./2.*(xn1-xn)*(yn1+yn));
}

/** 
This function compute the volume of revolution define by a segment of 2 points.
*/

double Volume (double xn, double xn1, double yn, double yn1) {
  return ( sq((yn1+yn)/2.)*(xn1-xn)*M_PI );
}

/**
At some point, we will need to remove part of the solution, to create a small
circle at the end intersection of the tail and the bubble. We first select a
point non the tail. Then, we find the nearest point to this one on the
bubble. Lastly, we find the circle that cross those 2 points and we generate
it.*/

double distancePoint (coord A, coord B) {
  return sqrt( sq(A.x - B.x) + sq(A.y - B.y));
}

/** findB is a function that will find the array position of the nearest point
of a list of coordinate.*/

int findB (coord A, double * xB, double * yB, int sizeB) {
  int compt = sizeB;
  coord B;

  double distanceABPrev = HUGE;

  B.x = xB[compt];
  B.y = yB[compt];
  double distanceAB = distancePoint(A, B);

  while (distanceAB<distanceABPrev){
    distanceABPrev = distanceAB;
    compt--;
    B.x = xB[compt];
    B.y = yB[compt];
    distanceAB = distancePoint(A, B);
  }

  return (compt);

}

/**
The function findCircle will find the circle coordinate (center position and
radius) based on two point A and B, but also the direction between A and the
center of the circle (alpha) (assuming that the other value of the vector is
1.*/

Circle findCircle (coord A, coord B, double alpha) {
  Circle C;
  double r = sqrt(sq (alpha)+1)*
          (sq (B.x - A.x) + sq (B.y - A.y))/
          (2.*((A.x - B.x)*alpha + B.y - A.y));

  C.x = A.x - r*alpha/(sqrt (sq(alpha)+1));
  C.y = A.y + r/(sqrt (sq(alpha)+1));
  C.r = r;

  return C;
}

/**
Computation of the spherical cap above the bubble. This function will return 
the volume of the cap. It will also return the center cooridnate and the 
radius of the cap.*/

double sphericalCap(double xi, double yi, double xa, double ya, int out, Circle * cap) {
  /**
  xi and yi: beginning of the cap. xa and ya: to have the curvature at this point.*/

  double alpha = xi - xa;
  double beta = yi - ya;

  /**
  c use to compute center and radius of circle*/
  double c = -alpha*xi - beta*yi;

  double xc = 0;
  double yc = -c/beta;
  double rc = sqrt(sq(xi)+sq(yi - yc));

  cap -> x = xc;
  cap -> y = yc;
  cap -> r = rc;

  double thetaInit = atan(-alpha/beta);

  /**
  Initialisation of the circle resolution*/

  int n = 1000;

  double x[n];
  double y[n];

  x[0] = xi;
  y[0] = yi;
  double thetaStep = (M_PI/2.-thetaInit)/(n-1);
  double theta = thetaInit;

  double sCapBub = 0;

  FILE * fp = fopen("cap", "w+");

  if (out)
    fprintf(fp, "%g %g\n", x[0], y[0]);
  // fprintf(stdout, "départ cap: %g %g\n", x[0], y[0]);


  for (int i = 1; i<n; i++){
    theta+=thetaStep;
    x[i] = rc*cos(theta)+xc;
    y[i] = rc*sin(theta)+yc;
    sCapBub+= Volume(y[i-1], y[i], x[i-1], x[i]);
    if (out)
      fprintf(fp, "%g %g\n", x[i], y[i]);
  }
  fclose(fp);
  return sCapBub;
}

/**
The function sphericalCapCoord work like the sphericalCap function. However, 
instead of returning the volume of the cap, it returns all the coordinate of 
the cap.

The int direction will switch the way we output the coordinate. If direction is
negative, we will return the coordinate starting from the axis of symmetry up 
to the end of the cap. Otherwise, the coordinate will be return from the end 
of the triple point (bubble, cap and free surface intersection) to the axis.

*/

coord * sphericalCapCoord(double xi, double yi, double xa, double ya, 
  int direcion) {
  /**
  xi and yi: beginning of the cap. xa and ya: to have the curvature at this point.*/

  double alpha = xi - xa;
  double beta = yi - ya;

  /**
  c use to compute center and radius of circle*/
  double c = -alpha*xi - beta*yi;

  double xc = 0;
  double yc = -c/beta;
  double rc = sqrt(sq(xi)+sq(yi - yc));

  double thetaInit = atan(-alpha/beta);

  /**
  Initialisation of the circle resolution*/

  int n = 1000;

  double x[n];
  double y[n];

  x[0] = xi;
  y[0] = yi;
  double thetaStep = (M_PI/2.-thetaInit)/(n-1);
  double theta = thetaInit;
  
  // fprintf(stdout, "départ cap: %g %g\n", x[0], y[0]);


  for (int i = 1; i<n; i++){
    theta+=thetaStep;
    x[i] = rc*cos(theta)+xc;
    y[i] = rc*sin(theta)+yc;
  }
  Array * b = array_new();
  if (direcion <0){
    for (int i = 1; i <= n; i++){
      coord p;
      p.x = x[n-i];
      p.y = y[n-i];
      array_append(b, &p, sizeof(coord));
    }
  }
  else {
    for (int i = 1; i <= n; i++){
      coord p;
      p.x = x[i];
      p.y = y[i];
      array_append(b, &p, sizeof(coord));
    }
  }

  coord q = {nodata};
  array_append(b, &q, sizeof(coord));
  return (coord *) array_shrink(b);
}


/**
#General algorithm

Now that we have defined all the function of the RK-4 algorithm, we can solve
the bubble shape. */

coord * input(double x[], double y[]) {
  Array * a = array_new();
  coord p = {0}, last, * la = NULL;
  long int ii =0;
  double test = 0.;
  while (test != stupid) {
    p.x = x[ii];
    p.y = y[ii];
    ii++;
    if (la) {
      array_append (a, la, sizeof(coord));
      array_append (a, &p, sizeof(coord));
    }
    last = p, la = &last;
    test = x[ii];
  }
  p.x = nodata;
  array_append (a, &p, sizeof(coord));
  return (coord *) array_shrink (a);
}

/**
This function take a lot of input  argument:
 * Bo: the Bond number define by using the curvature at the bottom of the 
 bubble as characteristic length
 * BondGhabache: a pointer that will allow us to obtain the Bond number based 
 on an elliptical approximation of the bubble for the bubble radius
 * effectiveBond: a pointer that will allow us to get the Bond number based on 
 the radius of the bubble
 * size: a pointer that allow us to obtain the size of the coordinate array
 * pointLess: the number a point we want to remove from the resulting shape (
 can be usefull if you want to add a spherical allow )
 * hollow: a pointer on a circle data structure. It allow you to obtain the 
 center and the radius of a small circle that you add where you remove the 
 points (should be use with the previous argument)
 * cap: a pointer that allow you to obtain the spherical cap position and its 
 radius
 * decal: an int that you can use if you want to obtain either the bubbleShape 
 and the free surface after (default case, decal = 0), only the free surface 
 and the spherical cap (decal = 1) or just the bubble and its cap (decal = -1)
 * h: if decal is equal to -1, you can translate the bubble of -h allong the 
 axis of symmetry
*/

coord * shape(double Bo, double * BondGhabache, double * effectiveBond, 
  double * phiTriple,  int * size, int pointLess, Circle * hollow,
  Circle * cap,  int decal, double h) {

  if (decal < 0)
    decal = -1;
  else if (decal>0)
    decal = 1;

  int fsize = ceil(2./fstep);
  int gsize = ceil(2./gstep);

  /**
  Initialization of the resolution. x0 should be 0. But we have a division by 
  x in the first step. So instead of having an error because of a division by 
  0, we set the initial value to something very small, close to 0.*/

  double x0 = 0.00000001;
  double z0 = 0.;
  double dz0 = 0.;

  double xn=x0;
  double zn=z0;
  double dzn=dz0;

  /**
  Before starting solving the differential equation, we create the variables.*/

  double k_1, k_2, k_3, k_4;
  double zn1, xn1, dzn1;

  /**
  Before we solve the equation, we pre-allocate the memory to store the
  solution. This will just concern the second part of the problem.*/

  double* xValue = NULL;
  double* zValue = NULL;
  double* dzdxValue = NULL;
  double* phi = NULL;

  double* xBegin = NULL;
  double* zBegin = NULL;

  xBegin = malloc(fsize *sizeof(double));
  zBegin = malloc(fsize *sizeof(double));
  xValue = malloc(gsize *sizeof(double));
  zValue = malloc(gsize *sizeof(double));
  dzdxValue = malloc(gsize *sizeof(double));
  phi = malloc(gsize *sizeof(double));

  /**
  As we mentioned above, we will use a Runge Kutta algorithm at the rank 4.*/

  int i = 0;
  xBegin[i] = 0.;
  zBegin[i] = 0.;

  double xBubbleMax = -HUGE;
  double zBubbleMax;

  while ( dzn<1.) {
    k_1 = k1(xn, zn, dzn, Bo, fCase, 0);
    k_2 = k2(xn, zn, dzn, Bo, k_1, fCase, 0);
    k_3 = k3(xn, zn, dzn, Bo, k_1, k_2, fCase, 0);
    k_4 = k4(xn, zn, dzn, Bo, k_2, k_3, fCase, 0);
    zn1 = yn1(zn, dzn, k_1, k_2, k_3, fCase);
    dzn1 = dyn1(dzn, k_1, k_2, k_3, k_4, fCase);
    xn1 = xn+fstep;

    i++;
    xBegin[i] = xn1;
    zBegin[i] = zn1;
    
    xn = xn1;
    zn = zn1;
    dzn = dzn1;
  }

  /**
  After the first part of the solveur, we add a "stupid" parameter. This will
  be used later when we would output the coordinate of the bubble shape into a 
  coord structure. */

  i++;
  xBegin[i] = stupid;
  zBegin[i] = stupid;

  /**
  We solved the first part of the problem. Now, we switch the axis.

  We initialize the new problem. We define the new parameters, dx/dz. The 
  other values are already defined from the previous loop. Because we changed 
  the reference axis when dz/dx = 1, the value of the derivative function 
  dx/dz is also 1.*/
  
  double dxn = 1.;
  double dxn1;

  xValue[0] = xn;
  zValue[0] = zn;
  dzdxValue[0] = 1.;

  /**
  phi corresponds to the angle between the local slope of the shape and the x 
  axi.*/

  phi[0] = atan(dzdxValue[0])*180./M_PI;

  int compt = 1;

  while (dxn >= (-1000.)) {
    k_1 = k1(zn, xn, dxn, Bo, gCase, 0);
    k_2 = k2(zn, xn, dxn, Bo, k_1, gCase, 0);
    k_3 = k3(zn, xn, dxn, Bo, k_1, k_2, gCase, 0);
    k_4 = k4(zn, xn, dxn, Bo, k_2, k_3, gCase, 0);
    xn1 = yn1(xn, dxn, k_1, k_2, k_3, gCase);
    dxn1 = dyn1(dxn, k_1, k_2, k_3, k_4, gCase);
    zn1 = zn+gstep;

    /**
    We store the value in our array. We will use them later to solve the tail 
    part of the geometry. The integer compt (from the French "compter" which 
    means to count) is used to identify the position in the array.*/

    xValue[compt] = xn1;
    zValue[compt] = zn1;
    dzdxValue[compt] = 1./dxn1;
    phi[compt] = atan(dzdxValue[compt])*180./M_PI;

    /**
    If the value of phi is below 0, we add 180 to have a positive angle. The
    angle phi corresponds to the local angle between the interface and the x
    axi.*/

    if (phi[compt]<0)
      phi[compt]+=180;

    /**     
    We want to have the equivalent Bond number of Ghabache et al. For
    that, we need the maximum value of x. Then, we will use the find function 
    to obtain the associated z value. The Bond we will use is linked to the 
    numerical Bond with the following relationship:

    $$
    Bo_{Ghabache} = Bo_{num}\times \frac{(x_{BubbleMax}^2z_{BubbleMax})^{2/3}}{R_0^2}
    $$

    Where $R_0$ is the curvature at the bottom of the bubble. Since we assume, 
    for this code that $R_0=1$, we are left with:

    $$
    Bo_{Ghabache} = Bo_{num}\times (x_{BubbleMax}^2z_{BubbleMax})^{2/3}
    $$*/

    if (xValue[compt]>= xBubbleMax) {
      xBubbleMax = xValue[compt];
      zBubbleMax = zValue[compt];
      *BondGhabache = Bo*pow(sq(xBubbleMax)*zBubbleMax,(2./3.));
    }
    /**
    Initialization of the parameters for the next loop.*/

    xn = xn1;
    zn = zn1;
    dxn = dxn1;
    compt++;
  }

  xValue[compt] = stupid;
  zValue[compt] = stupid;
  phi[compt] = stupid;

  /**
  As we did for the first part of the resolution, we add a "stupid" value at 
  the end of the array for x and z.

  At this point, we have solved the bubble interface. We now have to solve the 
  tail part. The tail interface is continuous with the bubble interface. We 
  don't know where the tail is supposed to start.

  For the initialization, we take 180 as max value and 90 as min value. They 
  correspond to the local angle of the interface where the tail start, with 
  the x-axis. This angle is always above 90° and below 180°. phiC is the value 
  of the angle we will use to solve the tail part of the function.*/

  double phiMax = 180.;
  double phiMin = 90.;
  double phiC = (phiMin+phiMax)/2.;

  int position = find(phiC, phi, gsize);

  double xInit = (xValue[position]+xValue[position+1])/2.;
  double zInit = (zValue[position]+zValue[position+1])/2.;
  double dzdxInit = (dzdxValue[position]+dzdxValue[position+1])/2.;

  /**
  The loop will work on phiMax-phiMin. We will solve the tail equation with 
  various initial condition until the difference between the 2 consecutive 
  initial value of x and z is below 0.5*10^-5. Since we don't have yet a value 
  for xInitOld and zInitOld, we will set them to 0*/

  double xInitOld = 0.;
  double zInitOld = 0.;

  /**
  We will store the value of x and z in an array. We will allocate the memory 
  in the while loop. We will also free this memory in the while loop, at the 
  beginning of the loop. For the initialization, we will create the pointer.*/

  

  double tailSize = SIZE/tstep;
  int tsize;
  tsize = ceil(tailSize);

  double* xTail = NULL;
  double* zTail = NULL;

  xTail = malloc(tsize *sizeof(double));
  zTail = malloc(tsize *sizeof(double));
  double ratio;

  double R, L = 0.;
  int j;

  while ((fabs(xInit-xInitOld)>stopCondition) && (fabs(zInit-zInitOld)>stopCondition)) {

    /**
    Initialization of the parameters R and L. R corresponds to the equivalent 
    radius of the bubble. L is the height of the tail at infinity.*/

    R = xInit/sind(phiC);
    L = 2./Bo*(2./R-1);

    /**
    Initialization of the parameters before the resolution. The step is the 
    same as for the case f.*/

    zn = zInit;
    xn = xInit;
    dzn = dzdxInit;

    xTail[0] = xn;
    zTail[0] = zn;
    double minZ = zn;

    j = 0;


    while (xn < 10) {
      j++;
      k_1 = k1(xn, zn, dzn, Bo, tailCase, L);
      k_2 = k2(xn, zn, dzn, Bo, k_1, tailCase, L);
      k_3 = k3(xn, zn, dzn, Bo, k_1, k_2, tailCase, L);
      k_4 = k4(xn, zn, dzn, Bo, k_2, k_3, tailCase, L);
      zn1 = yn1(zn, dzn, k_1, k_2, k_3, tailCase);
      dzn1 = dyn1(dzn, k_1, k_2, k_3, k_4, tailCase);
      xn1 = xn+tstep;

      /**
      Just after the resolution of the step, we compute the min value of the z 
      values, in order to already have it after the resolution of the tail. 
      This will be the test value for the next step. */

      minZ = minFunction(minZ, zn1);

      /**
      As previous, we store the value of the function in an array.*/

      xTail[j] = xn1;
      zTail[j] = zn1;
      xn = xn1;
      zn = zn1;
      dzn = dzn1;

      /** To avoid the problem related to an infinite number, we implement a
      test on dz, when dz is too big or too small.*/

      if (dzn1>1000.)
        xn = SIZE+1;

      if (dzn1<-1000.)
        xn = SIZE+1;
    }

    /**
    We test the min value of the tail. If it's below L, then the tail is 
    dropping under the axi z = L. We need to increase the angle. If not, it 
    means that the tail is going up to +infty. We need to reduce the angle*/

    if (minZ<L) {
      phiMin = phiC;
      phiC = (phiC+phiMax)/2.;
    }
    else {
      phiMax = phiC;
      phiC = (phiC+phiMin)/2.;
    } 

    /**
    We have solved the problem, we can now initialize the begining of the next 
    loop.*/

    xInitOld = xInit;
    zInitOld = zInit;

    position = find(phiC, phi, gsize);

    /**
    The variable "ratio" is used to determine how far is the new phiC from the 
    greatest lower value obtain via find. We will then use this ratio to 
    locate the x and z position of the beginning of the tail.*/

    ratio = ( phiC - phi[position] )/( phi[position+1] - phi[position] );

    xInit = xValue[position]+ratio*( xValue[position+1] - xValue[position] );
    zInit = zValue[position]+ratio*( zValue[position+1] - zValue[position] );
    dzdxInit = dzdxValue[position]+ratio*( dzdxValue[position+1] -
    dzdxValue[position]);

    /**
    At the end of the loop, we also add the "stupid" value at the end of the x 
    and z coordinate.*/

    xTail[j+1] = stupid;
    zTail[j+1] = stupid;
  }

  *phiTriple = phiC;

  /**
  The resolution is over. All we have to do now is to clear the useless value 
  in the array xValue and zValue and to return them to the main code.

  But first, we will compute some useful value for the post-treatment

  Those value will be return to the main code by using pointers.*/


  int tSize = ceil(SIZE/tstep);

  double* xFinal = NULL;
  double* yFinal = NULL;
  xFinal = malloc((tSize+fsize+gsize) *sizeof(double));
  yFinal = malloc((tSize+fsize+gsize) *sizeof(double));

  i = 1;

  /**
  The half volume of the bubble is the half volume of the spherical cap (
  computed below) and the volume of the bubble in the water.

  We first start with the spherical cap. */



  double xi,yi,xa,ya;
  xi = xValue[position];
  yi = zValue[position];
  xa = xValue[position-1];
  ya = zValue[position-1];

  #if WRONG
  double sBubble = -1/2.*sq(R)*sind(2*phiC);
  #else
  double sCap = sphericalCap(xi, yi, xa, ya, 0, cap);
  double sBubble = sCap;
  #endif

  /**
  We will reconstruct the final output. We will use the first 2 loops to 
  compute the volume of the bubble, in order to have the effective Bond 
  number of the geometry. */

  yFinal[0] = xBegin[0];
  xFinal[0] = zBegin[0]-L;

  while (xBegin[i] != stupid) {
    
    yFinal[i] = xBegin[i];
    xFinal[i] = zBegin[i]-L;

    #if WRONG
    sBubble+= surface(xFinal[i-1], xFinal[i], yFinal[i-1], yFinal[i]);
    #else
    sBubble+= Volume(xFinal[i-1], xFinal[i], yFinal[i-1], yFinal[i]);
    #endif
    i++;
  }

  int k;

  coord A;
  A.x = xTail[pointLess];
  A.y = zTail[pointLess];

  int jMoins = findB (A, xValue, zValue, position-1);

  for (k=0; k<(jMoins+1); k++) {
    yFinal[k+i] = xValue[k];
    xFinal[k+i] = zValue[k]-L;

    #if WRONG
    sBubble+= surface(xFinal[k+i-1], xFinal[k+i], yFinal[k+i-1], yFinal[k+i]);
    #else
    sBubble+= Volume(xFinal[k+i-1], xFinal[k+i], yFinal[k+i-1], yFinal[k+i]);
    #endif
  }

  int sizeBubbleTab = jMoins+i;
  int jPlus = pointLess;
  
  j = jPlus;
  #if InterfaceNoise
  int finalBegin = j+k+i-jPlus;
  // fprintf(stdout, "toto\n");
  #endif

  while (xTail[j] != stupid) {
    
    yFinal[j+k+i-jPlus] = xTail[j];
    xFinal[j+k+i-jPlus] = zTail[j]-L;
    j++;  
  }

  xFinal[j+k+i-jPlus] = stupid;
  yFinal[j+k+i-jPlus] = stupid;

  /**
  We will find the center of the circle C (and his radius). Before that, we
  shift the input data allong the $y$ axis from $-L$*/

  A.y = A.y - L;

  coord B;
  B.x = xValue[jMoins];
  B.y = zValue[jMoins] - L;

  /**
  alpha is the sloap value on the bubble.*/
  double alpha = dzdxValue[jMoins];

  /**   
  The function findCircle take 3 argument, a point on the bubble, a point
  on the tail and the sloap on the bubble.*/

  Circle C = findCircle (B, A, alpha);

  /**
  We compute the effective Bond number based on the surface of the bubble (
  assuming the geometry is axisymmetric. The value of the effective Bond is:

  $$
  Bo_{effective} = Bo (\frac{r_{effective}}{R_0})^2
  $$
  */

  #if WRONG
  double effectiveRadius = sqrt(2*sBubble/M_PI);
  #else
  double effectiveRadius = pow( (3./(4.*M_PI)*sBubble), 1./3.);
  #endif

  *effectiveBond = Bo*sq(effectiveRadius);

  /**   
  Before outputting the data, we have to rescale the shape. Actually, the
  dimensionless scale is based on the bottom curvature. We want to have a
  dimensionless scale based on the bubble radius. The equivalent radius of the
  bubble (if the bubble were a perfect sphere) should be 1.

  For that, we will apply the following transformation to the x and y 
  coordinate:

  $$
  x_{new} = x\frac{R_0}{r_{effective}}
  $$*/

  #if WRONG
  
  #else
  sCap = sphericalCap(xi/effectiveRadius, (yi-L)/effectiveRadius, 
    xa/effectiveRadius, (ya-L)/effectiveRadius
        , outCap, cap);
  #endif
  coord * capCoord = sphericalCapCoord(xi/effectiveRadius,
     (yi-L)/effectiveRadius, xa/effectiveRadius, (ya-L)/effectiveRadius, 
     -decal);


  /**
  We must not forgot to rescale the circle*/
  C.r = C.r/(effectiveRadius);
  C.x = C.x/effectiveRadius;
  C.y = C.y/effectiveRadius;

  /**
  We can now output this circle in the corresponding pointer, hollow.*/
  hollow -> x = C.x;
  hollow -> y = C.y;
  hollow -> r = C.r;
  
  int m = 0;
  
  if (decal == 0) {

    while (xFinal[m] != stupid) {
      yFinal[m] = yFinal[m]/(effectiveRadius);
      xFinal[m] = xFinal[m]/(effectiveRadius);
      m++;
    }
      // fprintf(stdout, "normalisation, you must take the last %f\n",effectiveRadius);
    
      /**
       We transformed our shape. We may have reduced his size. In order to
have a tail that will go at least at xn = 10 (which means yFinal must reach 10,
at least), we will prolong the tail, with its last value.*/
    
      double endTest = yFinal[m-1];
      double xEnd = xFinal[m-1];
    
      while (endTest < SIZE){
        xFinal[m] = xEnd;
        yFinal[m]= yFinal[m-1] + tstep;
        endTest+=tstep;
        m++;
      }


      *size = m;
      xFinal[m]=stupid;
      yFinal[m]=stupid;
    }
  else if (decal == 1) { 
      m = sizeBubbleTab;
      int i = 0;
      double xFinal2[1000 + tSize];
      double yFinal2[1000 + tSize];
  
        
      /**
      First, we add the spherical cap*/
      while(capCoord[i].x != nodata) {
        xFinal2[i] = capCoord[i].y;
        yFinal2[i] = capCoord[i].x;
        i++;
      }
      m+=3;
      while (xFinal[m] != stupid){
        yFinal2[i] = yFinal[m]/(effectiveRadius);
        xFinal2[i] = xFinal[m]/(effectiveRadius);
        m++;
        i++;
      }
      for (int j = 0; j<i-3; j++){
        xFinal[j] = xFinal2[j];
        yFinal[j] = yFinal2[j];
      }
  
      double endTest = yFinal[i-4];
      double xEnd = xFinal[i-4];
      
      while (endTest < SIZE){
        xFinal[i] = xEnd;
        yFinal[i]= yFinal[m-1] + tstep;
        endTest+=tstep;
        i++;
      }
      *size = i-4;
      xFinal[i-3]=stupid;
      yFinal[i-3]=stupid;
    }
  else if (decal == -1) {
  /**
  For this case, we only want the bubble and its cap. 
  You can translate this bubble by a distance h*/
  
    for (int i =0; i<sizeBubbleTab; i++){
      xFinal[i] = xFinal[i]/effectiveRadius-h;
      yFinal[i] = yFinal[i]/effectiveRadius;
    }
    int j = 0;
    while (capCoord[j].x != nodata) {
      xFinal[j+sizeBubbleTab] = capCoord[j].y-h;
      yFinal[j+sizeBubbleTab] = capCoord[j].x;
      j++;
    }
    *size = sizeBubbleTab+j;
    xFinal[j+sizeBubbleTab] = stupid;
    yFinal[j+sizeBubbleTab] = stupid;
  }

  #if InterfaceNoise
    double delta = SIZE/512.;
    // fprintf(stdout, "delta = %f\n",delta );
    double amplitude = 1e-3;
    i = finalBegin;
    double yNoiseInit = yFinal[i];
    double epsilon = velocityNoise(amplitude);
    while (xFinal[i] != stupid){
      if (yFinal[i]-yNoiseInit> delta/2.){
        epsilon = velocityNoise(amplitude);
        yNoiseInit = yFinal[i];
      }
      xFinal[i] = xFinal[i] + epsilon;
      i++;
    }
  #endif
  

  coord* sortie = NULL;
  
  sortie = input(xFinal, yFinal);  

  /**
  Before stopping the code, we free the memory :)*/

  free(xBegin);
  free(zBegin);
  free(xValue);
  free(zValue);
  free(dzdxValue);
  free(phi);
  free(xFinal);
  free(yFinal);
  free(xTail);
  free(zTail);

  return sortie; 
}

/**
#Use of this library

We use this library in a test case, to compute the shape of a bursting bubble:

 * [a bubble shape](bubble.c)

A more advance test case also use this library:

 * [volume of a bubble](bubbleRef.c)

This library is also called by:

 * [the bursting bubble](burstingBubble.c) 
 * [Input for a given Bond](findBond.h)

#Known issues

## Divergence of the tail

While solving the differential equations, if the input Bond number is too high,
to code will not reach the end of the solving range (the tail will diverge
before x=10). This is completely normal since the tail shape is highly 
dependent of the initial condition.

A possible solution is to decrease the value of the stop condition, but this
will increase the computation time. It's not really a problem: it's taking less
than 2 seconds on a laptop. */