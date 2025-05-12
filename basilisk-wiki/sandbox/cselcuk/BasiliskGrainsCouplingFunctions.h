/** 
# Helper functions for interfacing Basilisk/Grains3D 
*/
void printBasiliskDataStructure (struct BasiliskDataStructure * b) {
  printf("Grains nparticles = %lu\n", b->particlessize);
  
  /* Geometric aspect */
  printf("Grains radius = %f\n", b->rayon);
  printf("Grains center.x = %f\n", b->centreX);
  printf("Grains center.y = %f\n", b->centreY);
  printf("Grains center.z = %f\n", b->centreZ);
  printf("Grains ncorners = %d\n", b->ncorners);


  /* Coordinates of the corners, number of faces and their numerotation */

  printf("Grains allPoints = %lu\n", b->allPoints);
  printf("Grains allFaces = %lu\n", b->allFaces);

  for (int i = 0; i < (b->allPoints); i ++) {

# if dimension == 2
    printf ("Grains Coordinates of point %d is (%f,%f)\n", i, b->PointsCoord[i][0],b->PointsCoord[i][1]);
    
# elif dimension ==3

    printf ("Grains Coordinates of point %d is (%f,%f,%f)\n", i, b->PointsCoord[i][0], b->PointsCoord[i][1], b->PointsCoord[i][2]);
# endif
  } 
   
  /* Velocities */
  printf("Grains U.x = %f\n", b->vitesseTX);
  printf("Grains U.y = %f\n", b->vitesseTY);
  printf("Grains U.z = %f\n", b->vitesseTZ);
  
  printf("Grains w.x = %f\n", b->vitesseRX);
  printf("Grains w.y = %f\n", b->vitesseRY);
  printf("Grains w.z = %f\n", b->vitesseRZ);

  /* Physical properties */
  printf ("Grains mass = %f\n", b->masse);
  printf ("Grains rho_s = %f\n", b->masseVol);
  for (int k = 0; k < 6; k++) {
      printf ("Grains inertia = %g\n", b->inertie[k]);
    }
}

void UpdateParticlesBasilisk (struct BasiliskDataStructure * b, particle * p, const int m) {

  int r;
  for (int k = 0; k < m; k++) {
    
#if _MPI

    MPI_Bcast (&(b[k].particlessize), 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Geometric aspect */
    MPI_Bcast (&(b[k].centreX), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].centreY), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].centreZ), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].rayon), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].ncorners), 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Coordinates of the corners, number of faces and their numerotation */
    MPI_Bcast (&(b[k].allPoints), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].allFaces), 1, MPI_INT, 0, MPI_COMM_WORLD);

    r = b[k].allPoints;
  
    /* Allocate the structure in the other threads (Grains is serial) */
    if (pid() > 0) {
      b[k].PointsCoord = (double **) malloc (r*sizeof(double *));
#if dimension == 2 
      for (int i = 0; i < r; i++) 
	b[k].PointsCoord[i] = (double *) malloc (2 * sizeof(double));
#elif dimension == 3
      for (int i = 0; i < r; i++) 
	b[k].PointsCoord[i] = (double *) malloc (3 * sizeof(double));
#endif
    }

    for (int i = 0; i < r; i++) {
      MPI_Bcast (&(b[k].PointsCoord[i][0]), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast (&(b[k].PointsCoord[i][1]), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if dimension ==3
      MPI_Bcast (&(b[k].PointsCoord[i][2]), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    }

    /* Allocate the structures in the other threads (Grains is serial) */
    r = b[k].allFaces;
    if (pid() > 0) {
      b[k].FacesIndex = (long int **) malloc (r * sizeof(long int *));
      b[k].numPointsOnFaces = (long int *) malloc (r * sizeof(long int));
    }

    /* Broadcast-allocate-broadcast */
    for (int i = 0; i < r; i++) {
      MPI_Bcast (&(b[k].numPointsOnFaces[i]), 1, MPI_LONG, 0, MPI_COMM_WORLD);
      int rr = b[k].numPointsOnFaces[i];
      if (pid() > 0)
	b[k].FacesIndex[i] = (long int *)malloc(rr * sizeof(long int));
      MPI_Barrier (MPI_COMM_WORLD);

      for (int j = 0; j < rr; j++)
	MPI_Bcast (&(b[k].FacesIndex[i][j]), 1, MPI_LONG, 0, MPI_COMM_WORLD);
    }

    /* Velocities */
    MPI_Bcast (&(b[k].vitesseTX), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].vitesseTY), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].vitesseTZ), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].vitesseRX), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].vitesseRY), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].vitesseRZ), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Physical properties */
    MPI_Bcast (&(b[k].masseVol), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].masse), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i = 0; i < 6; i++)
      MPI_Bcast (&(b[k].inertie[i]), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

 
    /* Geometric aspect */ 
    GeomParameter gg; 
    coord cc = {0., 0., 0.};
    cc.x = b[k].centreX;
    cc.y = b[k].centreY;
    cc.z = b[k].centreZ;
    gg.center = cc;
    gg.radius = b[k].rayon;
    gg.ncorners = b[k].ncorners;
    gg.allPoints = b[k].allPoints;
    gg.allFaces = b[k].allFaces;

    /* We need the coordinates of the corners for a cube */

    /* Allocate on all threads */
    r = gg.ncorners;
    gg.cornersCoord = (double **) malloc (r*sizeof(double *));
    
#if dimension == 2 
    for (int i = 0; i < r; i++) 
      gg.cornersCoord[i] = (double *) malloc (2 * sizeof(double));
#elif dimension == 3
    for (int i = 0; i < r; i++) 
      gg.cornersCoord[i] = (double *) malloc (3 * sizeof(double));
#endif

    int ndim = 0;
    for (int i = 0; i < r; i++) {
#if dimension == 2 
      ndim  = 2; 
#elif dimension == 3
      ndim = 3;
#endif
      for (int j = 0; j < ndim; j++) {
	  gg.cornersCoord[i][j] = b[k].PointsCoord[i][j];
      }
    }

    /* We need the indices of the corners for a cube */
    /* Allocate */
    r = gg.allFaces;
    
    gg.cornersIndex = (long int **) malloc (r * sizeof(long int *));
    gg.numPointsOnFaces = (long int *) malloc (r * sizeof(long int));
    
    for (int i = 0; i < r; i++) {
      int rr = b[k].numPointsOnFaces[i];
      gg.cornersIndex[i] = (long int *) malloc (rr * sizeof(long int));
      gg.numPointsOnFaces[i] = rr;
      for (int j = 0; j < rr; j++) {
	gg.cornersIndex[i][j] = b[k].FacesIndex[i][j];
      }
    }
    
    p[k].g = gg;
    
#if DLM_Moving_particle    
    /* Velocities */
    cc.x = b[k].vitesseTX;
    cc.y = b[k].vitesseTY;
    cc.z = b[k].vitesseTZ;

#if TRANSLATION
    /* Save the previous translational velocity before updating */
    p[k].Unm1 = p[k].U;
    p[k].U = cc;
#endif
    cc.x = b[k].vitesseRX;
    cc.y = b[k].vitesseRY;
    cc.z = b[k].vitesseRZ;

#if ROTATION
    /* Save the previous angular velocity before updating */
    p[k].wnm1 = p[k].w;
    p[k].w = cc;
#endif
#endif
    
    /* Physical properties */
    p[k].M = b[k].masse;
    p[k].rho_s = b[k].masseVol;
    p[k].Vp = (p[k].M)/(p[k].rho_s);

    /* Inertia tensor: Ggrains store them as */
    /* inertie[0] = Ixx; */
    /* inertie[1] = Ixy; */
    /* inertie[2] = Ixz; */
    /* inertie[3] = Iyy; */
    /* inertie[4] = Iyz; */
    /* inertie[5] = Izz; */

    /* Basilisk stores these as */
    /* Ip[0] = Ixx */
    /* Ip[1] = Iyy */
    /* Ip[2] = Izz */
    /* Ip[3] = Ixy */
    /* Ip[4] = Ixz */
    /* Ip[5] = Iyz */
      
    /* Ixx */
    p[k].Ip[0] = b[k].inertie[0];
    /* Iyy */
    p[k].Ip[1] = b[k].inertie[3];
    /* Izz */
    p[k].Ip[2] = b[k].inertie[5];
    /* Ixy */
    p[k].Ip[3] = b[k].inertie[1];
    /* Ixz */
    p[k].Ip[4] = b[k].inertie[2];
    /* Iyz */
    p[k].Ip[5] = b[k].inertie[4];
  }
}


void UpdateBasiliskStructure (struct BasiliskDataStructure * b, particle * p, const int m) {

  coord U = {0., 0., 0.};
  coord w = {0., 0., 0.};
  
  for (int k = 0; k < m; k++) {
#if DLM_Moving_particle
#if TRANSLATION
    U = p[k].U;
#endif
#if ROTATION
    w = p[k].w;
#endif
#endif
    b[k].vitesseTX = U.x;
    b[k].vitesseTY = U.y;
    b[k].vitesseTZ = U.z;
    b[k].vitesseRX = w.x;
    b[k].vitesseRY = w.y;
    b[k].vitesseRZ = w.z;   
  }
}


void unallocateBasiliskDataStructure(struct BasiliskDataStructure * b, const int m) {

  double * c;
  
  for (int k = 0; k < m; k++) {
    /* Unallocate fields of BasiliskDataStructure that are no longer needed (we transferred everything on the particles structure */
    /* printf("addresse of b[k].PointsCoord = %p\n", (void *) b[k].PointsCoord); */
    
    for (signed long int i = 0; i < b[k].allPoints; i++) {
      c = b[k].PointsCoord[i];
      /* printf("addresse of b[k].PointsCoord[i] = %p\n",(void *) b[k].PointsCoord[i]); */
      /* printf("c = %p\n",(void *) c); */
      if(c)
	free(c);
    }

    /* printf("ok free(b[k].PointsCoord[i]) on thread %d\n",pid()); */
    free(b[k].PointsCoord);

    /* printf("ok free(b[k].PointsCoord) on thread %d\n",pid()); */
    
    for (signed long int i = 0; i < b[k].allFaces; i++) {
      free(b[k].FacesIndex[i]);
    }
    /* printf("ok free(b[k].FaceIndex[i]) on thread %d\n",pid()); */
    free(b[k].FacesIndex);
    free(b[k].numPointsOnFaces);
  }

}
