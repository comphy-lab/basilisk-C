struct BasiliskDataStructure {
  /* Number of particles */
  long int particlessize;

  /* Geometric aspect */
  double rayon;
  double centreX, centreY, centreZ;
  int ncorners;

  /* Coordinates of the corners, number of faces and their numerotation */
  long int allPoints;
  double ** PointsCoord;
  long int allFaces; 
  long int * numPointsOnFaces;
  long int ** FacesIndex;
  
  /* Velocities */
  double vitesseTX, vitesseTY, vitesseTZ;
  double vitesseRX, vitesseRY, vitesseRZ; 

  /* Physical properties */
  double masseVol, masse, inertie[6];
};
