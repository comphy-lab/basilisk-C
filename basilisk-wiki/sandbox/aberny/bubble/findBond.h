/**  
When we solve the equation of Young-Laplace for a bubble shape, we need to
define a Bond number. This number is based on the curvature of the bottom of 
the bubble, for calculus simplicity. However, we work with another Bond 
number, based on the volume of the bubble.

with this code, we introduce two function. The first one, findBond give us the 
Bond number based on the curvature at the bottom of a bubble, that correspond 
to the Bond number based on the volume of the bubble (case 1) or to the Bond 
number define by E. Ghabache (case 2)

The second function is computing the interface of a bubble, for the Bond 
number based on the bubble volume (the one we are working with, in the common 
case). */

#include "bubbleShape.h"

/** We will proceed by dichotomy. To avoid that the program is running for
nothing, we define the number of maximum iteration*/

#define MaxCompt 100

/**
The tol value corresponds to a relative tolerance in error for the value of 
the Bond number based on the bubble volume.*/

double tolBond = 0.02/100.;

double findBond(double Bond, int caseBond) {
  
  /**
  We will use the shape function (see bubble.c for the use of this function).*/

  double* BondGhabache = NULL;
  double BoGha = 0;

  double* BoEff = NULL;
  double Bo = 0;

  double* phi = NULL;
  double phiCirle = 0;
  phi = &phiCirle;

  int* size = NULL;
  int taille = 0;
  size = &taille;

  Circle* hollow = NULL;
  Circle fillet;
  hollow = &fillet;

  Circle* cap = NULL;
  Circle capTop;
  cap = &capTop;

  /**  
  We will work, by dichotomy. So, we define a maximum for the input Bond
  number. There is also a minimum (at 0, obviously). The first input value 
  will be in the mean value of the min and the max.

  It's important to note that the we cannot find value above 5 */

  double BoMax = 1000;
  double BoMin = 0;
  double BoInput = (BoMax + BoMin)/2.;

  /**
  We stop the function if the wanted Bond number is above 5
  */

  if (Bond > 30.){
    fprintf(stderr, "Wanted Bond number is to high (above 30). Error 50\n");
    exit(50);
  }

  BondGhabache = &BoGha;
  BoEff = &Bo;


  /**   
  We define the test value of the dichotomy algorithm. At the beginning,
  we fix it to 1 (to be sure that we will go at least once in the while loop)*/

  double test = 1;
  int compt = 0;

  if (caseBond == 3)
    return Bond;


  while (  test > tolBond) {
    coord* data;
    data = shape(BoInput, BondGhabache, BoEff, phi, size, 0, hollow, cap, 0, 0);
    compt++;
    free(data);
    // fprintf(stderr, "%f\n",BoInput );
    /**     
    Redefinition of the Input and the min/max value. If the Effective
    Bond is too high, we have to decrease the input value. Else, we have to 
    increase the input value.

    It's important to notice the bifurcation here. There are indeed two cases. 
    Or we want to find the Effective Bond number (which is most of our case 
    when doing numerical simulations). We could also want to find the Bond 
    number define by E. Ghabache (which is easier to measure when doing 
    experiments). For that, we just change the test, so we will work with the 
    Bond Ghabache instead of the effective Bond number.*/

    if (caseBond == 1) {
      test = fabs(Bo-Bond)/Bond;
      if (Bo > Bond) {
        BoMax = BoInput;
        BoInput = (BoInput+BoMin)/2.;
      }
      else {
        BoMin = BoInput;
        BoInput = (BoInput + BoMax)/2.;
      }
    }
    else if (caseBond == 2) {
      test = fabs(BoGha-Bond)/Bond;
      if (BoGha > Bond) {
        BoMax = BoInput;
        BoInput = (BoInput+BoMin)/2.;
      }
      else {
        BoMin = BoInput;
        BoInput = (BoInput + BoMax)/2.;
      } 
    }
    else {
      fprintf(stderr, "Bad case selection. End of the programm, error 52.\n");
      exit(52);
    }
    /**     
    If we didn't reach the Bond we are looking for after 50 iterations,
    we stop the calcul with an error message. The exit code is 51*/
    if (compt>MaxCompt){
      fprintf(stderr, "Bond not reach after %d iterations. End of the programm, error 51\n", compt);
      exit(51);
    }
  }
  
  /**   
  Once the calculation is over, we output the value we have to input in
  the solver to have a bubble with the wanted Bond number.

  The case 3 is when we want to use directly the Bond number define with the
  radius of curcature at the bottom of the bubble*/

  return BoInput;
}

/** The function shapeBond will return, in a coord structure, the shape of a
bubble, for a given Bond number. The Bond number is based on the volume of the
bubble.*/

coord * shapeBond(double Bond, Circle * hollow, Circle * cap) {
  double BoInput = findBond(Bond, 1);
  
  double* BondGhabache = NULL;
  double BoGha = 0;

  double* BoEff = NULL;
  double Bo = 0; 

  BondGhabache = &BoGha;
  BoEff = &Bo;

  double* phi = NULL;
  double phiCirle = 0;
  phi = &phiCirle;

  int* size = NULL;
  int taille = 0;
  size = &taille;

  coord* data;
  data = shape(BoInput, BondGhabache, BoEff, phi, size, 200, hollow, cap, 0, 0);

  return data;
}

coord * partialCircle(Circle cap, double h, double phi) {
  double thetaInit = (phi-90)*M_PI/180.;
  Array * a = array_new();
  int steps = 1000;
  double x[steps+1];
  double y[steps+1];
  for (int i = 0; i<=steps; i++) {
    double step = steps;
    double theta = thetaInit+i/step*(M_PI/2.-thetaInit);
    x[i] = cap.r*cos(theta)+cap.x;
    y[i] = cap.r*sin(theta)+cap.y;
  }

  for (int i = steps; i >0; i--) {
    coord p;
    p.x = x[i];
    p.y = y[i]-h;
    coord q;
    q.x = x[i-1];
    q.y = y[i-1]-h;
    array_append(a, &p, sizeof(coord));
    array_append(a, &q, sizeof(coord));
  }

  for (int i = 1; i < steps; i++) {
    coord p;
    p.x = x[i-1];
    p.y = y[i-1];
    coord q;
    q.x = x[i];
    q.y = y[i];
    array_append(a, &p, sizeof(coord));
    array_append(a, &q, sizeof(coord));
  }
  coord p = {nodata};
  array_append(a, &p, sizeof(coord));
  return (coord *) array_shrink(a);
}

coord * topCircle(Circle cap, double phi) {
  double thetaInit = (phi-90)*M_PI/180.;
  Array * a = array_new();
  int steps = 1000;
  double x[steps+1];
  double y[steps+1];
  for (int i = 0; i<=steps; i++) {
    double step = steps;
    double theta = thetaInit+i/step*(M_PI/2.-thetaInit);
    x[i] = cap.r*cos(theta)+cap.x;
    y[i] = cap.r*sin(theta)+cap.y;
  }

  for (int i = 1; i < steps; i++) {
    coord p;
    p.x = x[i-1];
    p.y = y[i-1];
    coord q;
    q.x = x[i];
    q.y = y[i];
    array_append(a, &p, sizeof(coord));
    array_append(a, &q, sizeof(coord));
  }
  coord p = {nodata};
  array_append(a, &p, sizeof(coord));
  return (coord *) array_shrink(a);
}

coord * bubbleOnly(double Bond, double h, int * size) {
  double BoInput = findBond(Bond, 1);
  
  double* BondGhabache = NULL;
  double BoGha = 0;

  double* BoEff = NULL;
  double Bo = 0; 

  double* phi = NULL;
  double phiCirle = 0;
  phi = &phiCirle;

  BondGhabache = &BoGha;
  BoEff = &Bo;

  Circle* hollow = NULL;
  Circle fillet;
  hollow = &fillet;

  Circle* cap = NULL;
  Circle topCap;
  cap = &topCap;

  coord* data;
  data = shape(BoInput, BondGhabache, BoEff, phi, size, 0, hollow, cap, -1, h);

  return data;
}

coord * capOnly(double Bond, double h) {
  double BoInput = findBond(Bond, 1);
  
  double* BondGhabache = NULL;
  double BoGha = 0;

  double* BoEff = NULL;
  double Bo = 0; 

  BondGhabache = &BoGha;
  BoEff = &Bo;

  double* phi = NULL;
  double phiCirle = 0;
  phi = &phiCirle;

  Circle* hollow = NULL;
  Circle fillet;
  hollow = &fillet;

  Circle* cap = NULL;
  Circle topCap;
  cap = &topCap;

  int* size = NULL;
  int taille = 0;
  size = &taille;

  coord* data;
  data = shape(BoInput, BondGhabache, BoEff, phi, size, 0, hollow, cap, -1, h);
  coord* sortie;
  sortie = partialCircle(topCap, h, phiCirle);

  return sortie;
}

coord * freeSurface(double Bond, int * size) {
  double BoInput = findBond(Bond, 1);
  
  double* BondGhabache = NULL;
  double BoGha = 0;

  double* BoEff = NULL;
  double Bo = 0; 

  BondGhabache = &BoGha;
  BoEff = &Bo;

  double* phi = NULL;
  double phiCirle = 0;
  phi = &phiCirle;

  Circle* hollow = NULL;
  Circle fillet;
  hollow = &fillet;

  Circle* cap = NULL;
  Circle topCap;
  cap = &topCap;

  coord * data;
  data = shape(BoInput, BondGhabache, BoEff, phi, size, 0, hollow, cap, 1, 0);

  return data;
}

coord * bubbleFree(double Bond, int * size) {
  double BoInput = findBond(Bond, 1);
coord * fusion(coord * list1, int size1, coord * list2, int size2) {
  Array * final = array_new();
  int i = 0;
  while(list1[i].x != nodata){
    array_append(final, &list1[i], sizeof(coord));
    i++;
  }

  if (i%2 == 1)
    array_append(final, &list1[2*size1-3], sizeof(coord));
  // array_append(final, &list2[0], sizeof(coord));

  // int debut = 0;
  // if ((list2[0].x == final[i+1].x) && list2[0].x == final[i+1].y)

  for (i = 0; i < 2*size2; i++) {
    if (list2[i].x != nodata && list2[i].y >= 0)
      array_append(final, &list2[i], sizeof(coord));
  }

  coord p = {nodata};
  array_append (final, &p, sizeof(coord));
  return (coord *) array_shrink (final);

}

coord * bubbleTot(double Bond, double thick, int * size) {
  int * sizeBubble;
  coord * bubble;
  int tailleBubble = 0;
  sizeBubble = &tailleBubble;

  bubble = bubbleOnly(Bond, thick, sizeBubble);

  coord * topSurf;
  int* sizeFreeSurface = NULL;
  int  tailleSurface = 0;
  sizeFreeSurface = &tailleSurface; 

  topSurf = freeSurface(Bond, sizeFreeSurface);

  int sizeTot = tailleBubble + tailleSurface+1;
  *size = sizeTot;

  return fusion(bubble, tailleBubble, topSurf, tailleSurface);
}


  double* BondGhabache = NULL;
  double BoGha = 0;

  double* BoEff = NULL;
  double Bo = 0;

  BondGhabache = &BoGha;
  BoEff = &Bo;

  double* phi = NULL;
  double phiCirle = 0;
  phi = &phiCirle;

  Circle* hollow = NULL;
  Circle fillet;
  hollow = &fillet;

  Circle* cap = NULL;
  Circle topCap;
  cap = &topCap;

  coord * data;
  data = shape(BoInput, BondGhabache, BoEff, phi, size, 0, hollow, cap, 0, 0);

  return data;
}

coord * fusion(coord * list1, int size1, coord * list2, int size2) {
  Array * final = array_new();
  int i = 0;
  while(list1[i].x != nodata){
    array_append(final, &list1[i], sizeof(coord));
    i++;
  }

  if (i%2 == 1)
    array_append(final, &list1[2*size1-3], sizeof(coord));
  // array_append(final, &list2[0], sizeof(coord));

  // int debut = 0;
  // if ((list2[0].x == final[i+1].x) && list2[0].x == final[i+1].y)

  for (i = 0; i < 2*size2; i++) {
    if (list2[i].x != nodata && list2[i].y >= 0)
      array_append(final, &list2[i], sizeof(coord));
  }

  coord p = {nodata};
  array_append (final, &p, sizeof(coord));
  return (coord *) array_shrink (final);

}

coord * bubbleTot(double Bond, double thick, int * size) {
  int * sizeBubble;
  coord * bubble;
  int tailleBubble = 0;
  sizeBubble = &tailleBubble;

  bubble = bubbleOnly(Bond, thick, sizeBubble);

  coord * topSurf;
  int* sizeFreeSurface = NULL;
  int  tailleSurface = 0;
  sizeFreeSurface = &tailleSurface; 

  topSurf = freeSurface(Bond, sizeFreeSurface);

  int sizeTot = tailleBubble + tailleSurface+1;
  *size = sizeTot;

  return fusion(bubble, tailleBubble, topSurf, tailleSurface);
}

/** 
When the simulation is in 3D, we need to generate a 3D bubble shape. We
first need to rotate the 2D shape. We then need to generate triangles. This 
will be done simultaneously. */

coord * shape3D(coord* dataShape) {
  int size = 0;

  while (dataShape[size].x != nodata){
    size++;
  }



  // coord rotateOrigin[size];
  Array * rotate = array_new();
  Array * sortie = array_new();

  for (int i = 0; i < size; i++)
    array_append(rotate, &dataShape[i], sizeof(coord));
    // rotateOrigin[i] = dataShape[i];
  coord p = {nodata};
  array_append (rotate, &p, sizeof(coord));
  coord * rotateOrigin = array_shrink(rotate);

  fprintf(stderr, "size %d\n",size);

  // int m = 0;
  // int j = 0;

  // while (dataShape[j].x != nodata ) {
  //   j++;
  // }

  // m = j/2;
  // coord rotateOrigin[m];


  // for (int i =0; i < m; i++)
  //   rotateOrigin[i] = dataShape[2*i];
  int j = 0;

  coord rotateCalc[size];
  coord rotatePrevious[size];

  for (int i = 0; i<=180; i++) {
    
    double theta = 2*i*M_PI/180.;
    double thetaPrev = 2*(i-1)*M_PI/180.;
    // fprintf(stdout, "%g %g\n", theta, cos(theta));

    for (j = 0; j<size; j++) {
      /**
      Generation of the rotating shape.*/
      rotateCalc[j].x = rotateOrigin[j].x;
      rotateCalc[j].y = rotateOrigin[j].y*cos(theta) 
                       -rotateOrigin[j].z*sin(theta);
      rotateCalc[j].z = rotateOrigin[j].y*sin(theta)
                       +rotateOrigin[j].z*cos(theta);

      /**
      Generation of the previous line*/
      rotatePrevious[j].x = rotateOrigin[j].x;
      rotatePrevious[j].y = rotateOrigin[j].y*cos(thetaPrev)
                       -rotateOrigin[j].z*sin(thetaPrev);
      rotatePrevious[j].z = rotateOrigin[j].y*sin(thetaPrev)
                       +rotateOrigin[j].z*cos(thetaPrev);
    }

    for (int j = 0; j<size; j+=2) {
      /**
      We just want to keep the point in box of 10 by 10 by 10. We know that 
      the point on the x axis are alwais in that box, but it's not the case 
      for the point on the axis $y$ and $z$*/
      #if Quarter
      bool testY = (fabs(rotateCalc[j].y)<6. && fabs(rotatePrevious[j].y)<6.);
      bool testZ = (fabs(rotateCalc[j].z)<6. && fabs(rotatePrevious[j].z)<6.);
      #else      
      bool testY = (fabs(rotateCalc[j].y)<3. && fabs(rotatePrevious[j].y)<3.);
      bool testZ = (fabs(rotateCalc[j].z)<3. && fabs(rotatePrevious[j].z)<3.);
      #endif

      if (testY && testZ){

        /**
        First triangle*/
        
        array_append (sortie, &rotatePrevious[j], sizeof(coord));
        array_append (sortie, &rotatePrevious[j+1], sizeof(coord));
        array_append (sortie, &rotateCalc[j], sizeof(coord));
        
        /**
        Second triangle*/
        
        array_append (sortie, &rotateCalc[j], sizeof(coord));
        array_append (sortie, &rotatePrevious[j+1], sizeof(coord));
        array_append (sortie, &rotateCalc[j+1], sizeof(coord));
      }
    }
  }
  // coord p = {nodata};
  array_append (sortie, &p, sizeof(coord));
  return (coord *) array_shrink (sortie);
}

// coord * shape3D(double Bond) {
//   int size = 0;



//   coord * bubble;
//   int* sizeBubble = NULL;
//   int tailleBubble = 0;
//   sizeBubble = &tailleBubble;

//   #if FILM 
//     bubble = bubbleTot (Bond, thick, sizeBubble);
//   #else
//     bubble = bubbleFree(Bond, sizeBubble);
//   #endif

//   while (bubble[size].x != nodata){
//     size++;
//   }

//   coord rotateOrigin[size];
//   Array * sortie = array_new();

//   for (int i = 0; i < size; i++)
//     rotateOrigin[i] = bubble[i];

//   fprintf(stderr, "size %d\n",size);

//   // int m = 0;
//   // int j = 0;

//   // while (dataShape[j].x != nodata ) {
//   //   j++;
//   // }

//   // m = j/2;
//   // coord rotateOrigin[m];


//   // for (int i =0; i < m; i++)
//   //   rotateOrigin[i] = dataShape[2*i];
//   int j = 0;

//   coord rotateCalc[size];
//   coord rotatePrevious[size];

//   for (int i = 0; i<=180; i++) {
    
//     double theta = 2*i*M_PI/180.;
//     double thetaPrev = 2*(i-1)*M_PI/180.;
//     // fprintf(stdout, "%g %g\n", theta, cos(theta));

//     for (j = 0; j<size; j++) {
//       /**
//       Generation of the rotating shape.*/
//       rotateCalc[j].x = rotateOrigin[j].x;
//       rotateCalc[j].y = rotateOrigin[j].y*cos(theta) 
//                        -rotateOrigin[j].z*sin(theta);
//       rotateCalc[j].z = rotateOrigin[j].y*sin(theta)
//                        +rotateOrigin[j].z*cos(theta);

//       /**
//       Generation of the previous line*/
//       rotatePrevious[j].x = rotateOrigin[j].x;
//       rotatePrevious[j].y = rotateOrigin[j].y*cos(thetaPrev)
//                        -rotateOrigin[j].z*sin(thetaPrev);
//       rotatePrevious[j].z = rotateOrigin[j].y*sin(thetaPrev)
//                        +rotateOrigin[j].z*cos(thetaPrev);
//     }

//     for (int j = 0; j<size; j+=2) {
//       /**
//       We just want to keep the point in box of 10 by 10 by 10. We know that 
//       the point on the x axis are alwais in that box, but it's not the case 
//       for the point on the axis $y$ and $z$*/

//       bool testY = (fabs(rotateCalc[j].y)<3. && fabs(rotatePrevious[j].y)<3.);
//       bool testZ = (fabs(rotateCalc[j].z)<3. && fabs(rotatePrevious[j].z)<3.);

//       if (testY && testZ){

//         /**
//         First triangle*/
        
//         array_append (sortie, &rotatePrevious[j], sizeof(coord));
//         array_append (sortie, &rotatePrevious[j+1], sizeof(coord));
//         array_append (sortie, &rotateCalc[j], sizeof(coord));
        
//         /**
//         Second triangle*/
        
//         array_append (sortie, &rotateCalc[j], sizeof(coord));
//         array_append (sortie, &rotatePrevious[j+1], sizeof(coord));
//         array_append (sortie, &rotateCalc[j+1], sizeof(coord));
//       }
//     }
//   }
//   coord p = {nodata};
//   array_append (sortie, &p, sizeof(coord));
//   return (coord *) array_shrink (sortie);
// }