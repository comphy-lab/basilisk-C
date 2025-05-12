/**
# Helper functions for the fictitious-domain method with distributed
  Lagrangian multipliers.
*/


/** Structure for the coordinates of the fictitious domain's
    boundary. */
typedef struct {
  double *x;
  double *y;
  double *z;
  int m;
  int nm;
} SolidBodyBoundary;

/** Geometric parameters of the particle. */
typedef struct {
  coord  center;
  double radius;
  int ncorners, allPoints, allFaces;
  double ** cornersCoord;
  long int ** cornersIndex;
  long int * numPointsOnFaces;
  
  /*3 principals vectors of the cube */
  coord u1, v1, w1;
  coord mins, maxs;
  
} GeomParameter;

typedef struct {
  SolidBodyBoundary s;
  GeomParameter g;
  GeomParameter gnm1;
#if DLM_Moving_particle
  double kn, en, vzero, wished_ratio;
  coord normalvector;
  coord gravity;
   double M, Ip[6], rho_s, Vp, DLMFD_couplingfactor;
#if TRANSLATION
  coord U, Unm1, qU, tU;
#endif
#if ROTATION
  coord w, wnm1, qw, tw;
#endif
#endif
  size_t pnum, iswall, iscube;
  coord wallmax, wallmin, wallpos;
  coord imposedU, imposedw;
  Cache Interior;
  Cache reduced_domain;
  coord adforce;
  long tcells, tmultipliers;
} particle;

/** Function that allocates in memory the SolidBodyBoundary
    structure. */
void allocate_SolidBodyBoundary (SolidBodyBoundary *sbm, const int m) {
  sbm->x = (double*) calloc( m, sizeof(double)); 
  sbm->y = (double*) calloc( m, sizeof(double));
  
  if ( dimension == 3 )
    sbm->z = (double*) calloc( m, sizeof(double));
  else
    sbm->z = NULL;

  sbm->m = m;
}

void reallocate_SolidBodyBoundary (SolidBodyBoundary *sbm, const int m) {
  sbm->x = (double*) realloc (sbm->x, m*sizeof(double)); 
  sbm->y = (double*) realloc (sbm->y, m*sizeof(double));
  
  if (dimension == 3)
    sbm->z = (double*) realloc (sbm->z, m*sizeof(double));
  
  sbm->m = m;
}

/** Function that de-allocates in memory the SolidBodyBoundary
    structure. */
void free_SolidBodyBoundary (SolidBodyBoundary *sbm) {
  free(sbm->x); sbm->x = NULL;
  free(sbm->y); sbm->y = NULL;
  if ( dimension == 3 ){free(sbm->z); sbm->z = NULL;}
}

/** Function that allocates an initial Cache structure */
void allocate_Cache (Cache * p) {
  p->n = 0;
  p->nm = 0;
  if (p->n >= p->nm) {
    p->p = (Index *) calloc (5, sizeof(int));
  }
}

/** Set of functions for the sphere/circle as fictitious domain */
void compute_nboundary (const GeomParameter gcp, int * nb) {
  coord pos[6];
  Cache poscache = {0};
  Point lpoint;
   
  pos[0].x = gcp.center.x + gcp.radius + X0;
  pos[0].y = gcp.center.y + Y0;

  pos[1].x = gcp.center.x - gcp.radius + X0;
  pos[1].y = gcp.center.y + Y0;

  pos[2].x = gcp.center.x + X0;
  pos[2].y = gcp.center.y + gcp.radius + Y0;

  pos[3].x = gcp.center.x + X0;
  pos[3].y = gcp.center.y - gcp.radius + Y0;


  
#if dimension == 3
  pos[0].z = gcp.center.z + Z0;
  pos[1].z = gcp.center.z + Z0;
  pos[2].z = gcp.center.z + Z0;
  pos[3].z = gcp.center.z + Z0;
  pos[4].x = gcp.center.x + X0;
  pos[4].y = gcp.center.y + Y0;
  
  pos[4].z = gcp.center.z + gcp.radius + Z0;
  pos[5].x = gcp.center.x + X0;
  pos[5].y = gcp.center.y + Y0;
  pos[5].z = gcp.center.z - gcp.radius + Z0;
#endif
  
  *nb = 0;
  lpoint.level = -1;
  int ip = 0;
  
  /* We assume that around the sphere at least one point has to be at
     maximum level of refinement. */
  for (ip = 0; ip < dimension*2; ip++) {
#if dimension == 3
    lpoint = locate (pos[ip].x, pos[ip].y, pos[ip].z);
#elif dimension ==2 
    lpoint = locate (pos[ip].x, pos[ip].y);
#endif
    if (lpoint.level == depth())
      break;
  }
  
  /** Multiple threads can have a different point at the same level of
      refinement. This is fine, they will do the same computation for
      *nb. This problem does not exist in serial but the code below
      still works. */
  if (lpoint.level == depth()) {
           
    /** Only this thread creates the Cache ... */
    cache_append(&poscache, lpoint, 0);

    /** and only this thread computes the number of boundary points
	... */
    foreach_cache(poscache) {
#if dimension == 2
      *nb += (int){floor(gcp.radius*2*pi/(sqrt(2)*Delta))};
      /* Debug periodic boundaries, we use only one multiplier */
      /* *nb =1.; */
#elif dimension == 3
      *nb += (int){floor(pow( 3.809 * gcp.radius / (sqrt(3)*Delta), 2.))};
#endif
    }
    /** and finally, this thread destroys the cache. */
    
    free(poscache.p);
  }

#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
  mpi_all_reduce (*nb, MPI_INT, MPI_MAX);
#endif
    
  
  if(*nb == 0)
    printf("nboundary = 0: No boundary points !!!\n");

}

void create_FD_Interior_Particle (particle * p, vector Index_lambda, vector shift, scalar flag) {
  GeomParameter gci = p->g;
  Cache * fd = &(p->Interior);
  
  /** Create the cache for the interior points */
  foreach() {
#if dimension == 2
    if ((sq(x - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0))) < sq(gci.radius)) {
      cache_append (fd, point, 0);
      
      /* tagg cell with the number of the particle */
      if ((int)Index_lambda.y[] == -1)
	Index_lambda.y[] = p->pnum; 
    }

    if (Period.x) {
      if ((sq(x - L0 - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0))) < sq(gci.radius)) {
    	cache_append (fd, point, 0);
	
    	/* tagg cell with the number of the particle */
    	if ((int)Index_lambda.y[] == -1) 
    	  Index_lambda.y[] = p->pnum;
	if ((int)Index_lambda.y[] == p->pnum)
	  shift.x[] = L0;
      }
      if ((sq(x + L0 - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0))) < sq(gci.radius)) {
    	cache_append (fd, point, 0);

	/* tagg cell with the number of the particle */
    	if ((int)Index_lambda.y[] == -1) 
    	  Index_lambda.y[] = p->pnum;
	if ((int)Index_lambda.y[] == p->pnum)
	  shift.x[] = -L0;
      }
    }
   
    if (Period.y) {
      if ((sq(x - (gci.center.x + X0)) + sq(y - L0 - (gci.center.y + Y0))) < sq(gci.radius)) {
    	cache_append (fd, point, 0);
	
	/* tagg cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.y[] = L0;
	}
      }
      if ((sq(x - (gci.center.x + X0)) + sq(y + L0 - (gci.center.y + Y0))) < sq(gci.radius)) {
    	cache_append (fd, point, 0);
	
	/* tagg cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.y[] = -L0;
	}
      }
    }


    
#endif




    
#if dimension == 3
    if ((sq(x - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0)) + sq(z - (gci.center.z + Z0))) < sq(gci.radius)) {
      cache_append (fd, point, 0);
      /* tagg cell with the number of the particle */
      if ((int)Index_lambda.y[] == -1)
	Index_lambda.y[] = p->pnum;
    }
    if (Period.x) {
      if ((sq(x - L0 - (gci.center.x + X0) ) + sq(y - (gci.center.y + Y0)) + sq(z - (gci.center.z + Z0))) < sq(gci.radius)) {
    	cache_append (fd, point, 0);
    	/* tagg cell with the number of the particle */
    	if ((int)Index_lambda.y[] == -1) 
    	  Index_lambda.y[] = p->pnum;
	if ((int)Index_lambda.y[] == p->pnum)
	  shift.x[] = L0;
      }
	if ((sq(x + L0 - (gci.center.x + X0) ) + sq(y - (gci.center.y + Y0)) + sq(z - (gci.center.z + Z0))) < sq(gci.radius)) {
    	cache_append (fd, point, 0);
    	/* tagg cell with the number of the particle */
    	if ((int)Index_lambda.y[] == -1) 
    	  Index_lambda.y[] = p->pnum;
	if ((int)Index_lambda.y[] == p->pnum)
	  shift.x[] = -L0;
      }
    }
   
    if (Period.y) {
      if ((sq(x - (gci.center.x + X0)) + sq(y - L0 - (gci.center.y + Y0)) + sq(z - (gci.center.z + Z0))) < sq(gci.radius)) {
    	cache_append (fd, point, 0);
	/* tagg cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.y[] = L0;
	}
      }
      if ((sq(x - (gci.center.x + X0)) + sq(y + L0 - (gci.center.y + Y0)) + sq(z - (gci.center.z + Z0))) < sq(gci.radius)) {
    	cache_append (fd, point, 0);
	/* tagg cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.y[] = -L0;
	}
      }
    }

    if (Period.z) {
      if ((sq(x - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0)) + sq(z - L0 - (gci.center.z + Z0))) < sq(gci.radius)) {
    	cache_append (fd, point, 0);
	/* tagg cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.z[] = L0;
	}
      }
      if ((sq(x - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0)) + sq(z + L0 - (gci.center.z + Z0))) < sq(gci.radius)) {
    	cache_append (fd, point, 0);
	/* tagg cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.z[] = -L0;
	}
      }
    }
#endif
    
  }
    cache_shrink (fd);
}

void create_FD_Boundary_Particle (GeomParameter gcp, SolidBodyBoundary * dlm_bd, const int nsphere, vector pshift) {

  int m = nsphere;
  coord pos;
  Point lpoint;
  
  coord shift = {0., 0., 0.};
  coord ori = {X0, Y0, Z0};
 
#if dimension == 2
  int i;
  double theta[m];
  double radius  = gcp.radius;
  coord fact_theta[m];

  for (i = 0; i < m; i++) {
    theta[i] = i*2*pi/m; 
    fact_theta[i].x = cos (theta[i]);
    fact_theta[i].y = sin (theta[i]);
    
    /* Assign positions x,y on a circle's boundary */ 
    pos.x = radius*fact_theta[i].x + gcp.center.x + X0;
    pos.y = radius*fact_theta[i].y + gcp.center.y + Y0;

    /* Check if the point falls outside of the domain */
    foreach_dimension() {
      if (Period.x) {
	shift.x = 0.;
	if (pos.x > (L0 + ori.x)) {
	  pos.x -= L0;
	  shift.x = -L0;
	} 
	if (pos.x < (0. + ori.x)) {
	  pos.x += L0;
	  shift.x = L0;
	}
      }
      dlm_bd->x[i] = pos.x; 
    }

   
    Cache poscache = {0};
    lpoint = locate (pos.x, pos.y);

    if (lpoint.level > -1) {
	
      cache_append (&poscache, lpoint, 0);

      foreach_cache (poscache) {
	foreach_dimension() {
	  pshift.x[] = shift.x;
	}
      }
      
      free (poscache.p);
    }
  }
  
#elif dimension == 3

  /* Lay out points on a sphere with a Z oriented spiral scheme */
  /* This portion of code is recovered from Peligriff and adapted */
  /* More information of this method can be found at: */
  /* Saff, E. & Kuijlaars, A. Distributing many points on a sphere The
     Mathematical Intelligencer, 1997 */
  
  double spiral_spacing_correction = 0;
  foreach_level(depth()) 
     spiral_spacing_correction = 2.*Delta/sqrt(3.); 
  
  double hydro_radius = gcp.radius;
  size_t k;
  
  /* Number of points */
  size_t NSpiral = m;      
    
  /* Spiral points construction */
  double hk, thetak, phik, phikm1 = 0., TwoPi = 2. * pi, 
    Cspiral = 3.6 / sqrt( NSpiral ), dphi = 0.;

  for (k = 0; k < NSpiral; ++k) {
    hk = - 1. + 2. * (double)(k) / ( NSpiral - 1. );
    thetak = acos( hk ) ;
    if ( k == 0 ) {
      phik = 0.;
      thetak = pi - 0.5 * spiral_spacing_correction / hydro_radius ;
    }      
    else if ( k == NSpiral - 1 ) {
      phik = phikm1 + 1. * dphi; 
      if ( phik > TwoPi ) phik -= TwoPi ;
      thetak = 0.5 * spiral_spacing_correction / hydro_radius ;        
    }          
    else {
      dphi = Cspiral / sqrt( 1. - hk * hk ) ;
      phik = phikm1 + dphi ;
      if ( phik > TwoPi ) phik -= TwoPi ;  
    }   

    phikm1 = phik ;
         
    if ( k == 1 ) thetak -= 0.4 * spiral_spacing_correction / hydro_radius ;
    
    if ( k == NSpiral - 2 ) {
      phik -= 0.1 * dphi ;
      thetak += 0.25 * spiral_spacing_correction / hydro_radius ; 
    }
      
    pos.x = hydro_radius * cos( phik ) * sin( thetak ) + gcp.center.x + X0;
    pos.y = hydro_radius * sin( phik ) * sin( thetak ) + gcp.center.y + Y0;
    pos.z = hydro_radius * cos( thetak ) + gcp.center.z + Z0;

    
    
    foreach_dimension() {
      if (Period.x) {
	shift.x = 0.;
	if (pos.x > (L0 + ori.x)) {
	  pos.x -= L0;
	  shift.x = -L0;
	} 
	if (pos.x < (0. + ori.x)) {
	  pos.x += L0;
	  shift.x = L0;
	}
      }
    }
    
    
    Cache poscache = {0};
    /* find the cell of this multiplier */
    lpoint = locate (pos.x, pos.y, pos.z);

    if (lpoint.level > -1) {
	
      cache_append(&poscache, lpoint, 0);

      foreach_cache(poscache) {
	foreach_dimension() {
	  pshift.x[] = shift.x;
	}
      }
      
      free(poscache.p);
    }
    

      
    dlm_bd->x[k] = pos.x;
    dlm_bd->y[k] = pos.y;
    dlm_bd->z[k] = pos.z;

  }
#endif
  boundary((scalar*){pshift});
}

/** Set of functions for the walls as fictitious domain */
void create_FD_Interior_Wall (particle *p, vector index, scalar flag) {
  coord wallmin, wallmax; 
  Cache * c;
  
  /** Create the cache of the interior points for a wall*/
  
  c = &(p->Interior);
  wallmin = p->wallmin;
  wallmax = p->wallmax;
  foreach() {	   
#if dimension == 2
    if ((x >= wallmin.x) && (y >= wallmin.y) && (x <= wallmax.x) && (y <= wallmax.y)) {
      index.y[] = p->pnum;
      cache_append (c, point, 0);
    }
#elif dimension == 3
    if ((x > wallmin.x) && (y > wallmin.y) && (z > wallmin.z) && (x < wallmax.x) && (y < wallmax.y) && (z < wallmax.z)) {
      if ((int)index.y[] == -1 && flag[] < 1)
	index.y[] = p->pnum;
      cache_append (c, point, 0);
    }
#endif
  }
  cache_shrink (c);    
}

void create_FD_Boundary_Wall (particle *p) {
  
  SolidBodyBoundary * sbm;
  int m = 0;
  coord  wallpos;

  double h = 0;
  foreach_level(depth())
    h = Delta;
  
  sbm = &(p->s);
  m = sbm->m;
  wallpos = p->wallpos;
  for (int i = 0; i < m; i++) {

    sbm->x[i] = (3.*h/4.) + i*h;
    sbm->y[i] = wallpos.y;
  }

  for (int i = m/2; i < m; i++) {
    sbm->x[i] -= h/2.;
  }
}

void compute_nboundary_Wall (coord wallpos, int * nb) {

  double mindelta = 0.;
  foreach_level(depth())
    mindelta = Delta;

  *nb = floor(wallpos.x/mindelta);
    
}

/** Set of functions for the cube as fictitious domain */
void compute_nboundary_Cube_v2 (GeomParameter * gcp, int * nb, int * lN) {
  Cache poscache = {0};
  Point lpoint;
  *nb = 0;
  coord pos = {0., 0., 0.};
  int ip = 0;

  while (*nb == 0) {
    pos.x = gcp->cornersCoord[ip][0];
    pos.y = gcp->cornersCoord[ip][1];
    
#if dimension == 3
    pos.z = gcp->cornersCoord[ip][2];
    lpoint = locate(pos.x, pos.y, pos.z);
#elif dimension ==2
    lpoint = locate(pos.x, pos.y);
#endif
   

    /** Only one thread has the point in its domain (works in serial
	too). */
    if (lpoint.level > -1) {
    
      /** Only this thread creates the Cache ... */
      cache_append(&poscache, lpoint, 0);

      /** and only this thread computes the number of boundary points
	  ... */
    
      /* Grains sends an equivalent radius that is  multiplied by sqrt(3)/2 */
      double lengthedge = 2.*(gcp->radius)/sqrt(3.);
    
      /* printf ("lengthedge = %f on thread %d\n", lengthedge, pid() ); */
      foreach_cache (poscache) {
	*lN = floor (sqrt(3.)*lengthedge/(2.*Delta));
      }
#if dimension == 2
      /* number of points required for the 4 edges of the square */
      *nb += (*lN-2)*4;
      /* number of points required for the face */
      *nb += (*lN-2)*(*lN-2);
      /* number of points required for the 4 corners */
      *nb += 4;
#elif dimension == 3
      /* number of points required for the 12 edges of the cube */
      *nb += (*lN-2)*12;
      /* number of points required for the 6 faces of the cube */
      *nb += 6*(*lN-2)*(*lN-2);
      /* number of points required for the 8 corners */
      *nb += 8;
#endif
      
      /** and finally, this thread destroys the cache. */
      free(poscache.p);
    }
  
#if _MPI
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_all_reduce(*nb, MPI_INT, MPI_MAX);
    mpi_all_reduce(*lN, MPI_INT, MPI_MAX);
#endif
    if (ip < gcp->ncorners)
      ip++;
    else
      break;
  }
  
  /* printf("lN number of points for a complete edge:%d\n",*lN); */
  if (*nb == 0)
    fprintf(stderr,"nboundary = 0: No boundary points for the cube/square !!!\n");
}

void distribute_points_edge_Cube_v2 (coord const corner1, coord const corner2, SolidBodyBoundary * dlm_bd, int const lN, int const istart) {
# if dimension == 3 
  if (lN > 0) {
    coord dinc; 

    foreach_dimension() {
      dinc.x = (corner2.x - corner1.x)/(lN-1);
    }
    /* printf( "length edge corner to corner = %f\n", sqrt( sq(corner2.x - corner1.x) + sq(corner2.y - corner1.y)+ sq(corner2.z - corner1.z))); */
    for (int i = 1; i <= lN-2; i++) {
      dlm_bd->x[istart + i -1] = corner1.x + (double)i * dinc.x; 
      dlm_bd->y[istart + i -1] = corner1.y + (double)i * dinc.y;
      dlm_bd->z[istart + i -1] = corner1.z + (double)i * dinc.z;
    }
  }
#endif
}

bool is_it_in_cube_v2 (coord * u, coord * v, coord * w, coord * mins, coord * maxs, coord * checkpt) {
  
  /* a point x with coord (x,y,z) lies in the rectangle if the 3
     following conditions are satisfied */
  /* 1- u.p0 <= u.x <= u.p3 */
  /* 2- v.p0 <= v.x <= v.p1 */
  /* 3- w.p0 <= w.w <= w.p4 */
  double x = checkpt->x;
  double y = checkpt->y;
  double z = checkpt->z;
  double checkval = 0.;
  coord u1 = *u;
  coord v1 = *v;
  coord w1 = *w;
  coord cpt = {0., 0., 0.};
  bool isin = false;
  
  checkval = u1.x*x + u1.y*y + u1.z*z;
  cpt.x = checkval;

  checkval = v1.x*x + v1.y*y + v1.z*z;
  cpt.y = checkval; 

  checkval = w1.x*x + w1.y*y + w1.z*z;
  cpt.z = checkval;

  if ((cpt.x >= maxs->x) && (cpt.x <= mins->x) && (cpt.y >= maxs->y) && (cpt.y <= mins->y) && (cpt.z >= maxs->z) && (cpt.z <= mins->z) )
    {
      isin = true;
    }
  return isin;
}

void compute_principal_vectors_Cubes (particle * p) {

  GeomParameter * gcp = &(p->g);
  int nfaces = gcp->allFaces;
  int npoints;

  /* get the 3 directions u1, u2, u3 of the cube with such that */
  /* u1 = corner0 - corner3; */
  /* u2 = corner0 - corner1; */
  /* u3 = corner0 - corner4; */
  coord corner0 = {0, 0, 0}, corner1 = {0, 0, 0}, corner3 = {0, 0, 0}, corner4 = {0, 0, 0};
  coord u1, v1, w1;
  coord mins, maxs;
  
  /* find the coordinates of these 4 corners */
  for (int i = 0; i < nfaces; i++) {
    npoints = gcp->numPointsOnFaces[i];
    
    for (int j = 0; j < npoints; j++) {
      /* printf("index = %lu\n", gcp->cornersIndex[i][j]); */
      if (gcp->cornersIndex[i][j] == 0) {
	long int ii = gcp->cornersIndex[i][j];
	corner0.x = gcp->cornersCoord[ii][0];
	corner0.y = gcp->cornersCoord[ii][1];
	corner0.z = gcp->cornersCoord[ii][2];
	/* printf("index = %lu with coordinate (%f,%f,%f)\n", gcp->cornersIndex[i][j], corner0.x, corner0.y, corner0.z); */
      }
      if (gcp->cornersIndex[i][j] == 1) {
	long int ii = gcp->cornersIndex[i][j];
	corner1.x = gcp->cornersCoord[ii][0];
	corner1.y = gcp->cornersCoord[ii][1];
	corner1.z = gcp->cornersCoord[ii][2];
	/* printf("index = %lu with coordinate (%f,%f,%f)\n", gcp->cornersIndex[i][j], corner1.x, corner1.y, corner1.z); */
      }
      if (gcp->cornersIndex[i][j] == 3) {
	long int ii = gcp->cornersIndex[i][j];
	corner3.x = gcp->cornersCoord[ii][0];
	corner3.y = gcp->cornersCoord[ii][1];
	corner3.z = gcp->cornersCoord[ii][2];
	/* printf("index = %lu with coordinate (%f,%f,%f)\n", gcp->cornersIndex[i][j], corner3.x, corner3.y, corner3.z); */
      }
      if (gcp->cornersIndex[i][j] == 4) {
	long int ii = gcp->cornersIndex[i][j];
	corner4.x = gcp->cornersCoord[ii][0];
	corner4.y = gcp->cornersCoord[ii][1];
	corner4.z = gcp->cornersCoord[ii][2];
	/* printf("index = %lu with coordinate (%f,%f,%f)\n", gcp->cornersIndex[i][j], corner4.x, corner4.y, corner4.z); */
      }
    }  
  }

  foreach_dimension() {
    u1.x = corner0.x - corner3.x;
    v1.x = corner0.x - corner1.x;
    w1.x = corner0.x - corner4.x;
    mins.x = 0.;
    maxs.x = 0.;
  }
  
  gcp->u1 = u1;
  gcp->v1 = v1;
  gcp->w1 = w1;
  
  double minval = 0., maxval = 0.;

  foreach_dimension() {
    minval += u1.x*corner0.x;
    maxval += u1.x*corner3.x;
  }

  mins.x = (minval); maxs.x = (maxval);

  minval = 0; maxval = 0;
  foreach_dimension() {
    minval += v1.x*corner0.x;
    maxval += v1.x*corner1.x;
  }
  
  mins.y = (minval); maxs.y = (maxval);

  minval = 0; maxval = 0;
  foreach_dimension() {
    minval += w1.x*corner0.x;
    maxval += w1.x*corner4.x;
  }
  
  mins.z = (minval); maxs.z = (maxval);


  gcp->mins = mins;
  gcp->maxs = maxs;
  
}

/* void compute_principal_mins_maxs_Cubes (particle * p) { */
/*   GeomParameter * gcp = &(p->g); */
/*   coord u1 = gcp->u1; */
/*   coord v1 = gcp->v1; */
/*   coord w1 = gcp->w1; */
  
/* } */

void create_FD_Boundary_Cube_v2 (GeomParameter * gcp, SolidBodyBoundary * dlm_bd, const int m, const int lN, vector pshift) {

  int nfaces = gcp->allFaces;
  int iref, i1, i2, ichoice;

  ichoice = 0;
  int isb = 0;
  int npoints;

  /* Add first interrior points on surfaces */
  for (int i = 0; i < nfaces; i++) {
    npoints = gcp->numPointsOnFaces[i];
    
    iref = gcp->cornersIndex[i][ichoice];
    i1 = gcp->cornersIndex[i][ichoice + 1];
    i2 = gcp->cornersIndex[i][npoints-1];

    /* printf("on face %d, iref = %d, i1 = %d, i2 = %d\n", i, iref, i1, i2); */
    coord refcorner = {gcp->cornersCoord[iref][0], gcp->cornersCoord[iref][1],gcp->cornersCoord[iref][2]} ; 

    coord dir1 = {gcp->cornersCoord[i1][0], gcp->cornersCoord[i1][1],gcp->cornersCoord[i1][2]};

    coord dir2 = {gcp->cornersCoord[i2][0], gcp->cornersCoord[i2][1],gcp->cornersCoord[i2][2]};
    
    /* printf ("corresponding coord for ref = (%f,%f,%f)\n",refcorner.x, refcorner.y, refcorner.z); */
    /* printf ("corresponding coord for i1 = (%f,%f,%f)\n", dir1.x, dir1.y, dir1.z); */
    /* printf ("corresponding coord for i2 = (%f,%f,%f)\n", dir2.x, dir2.y, dir2.z); */
    
    foreach_dimension() {
      dir1.x -= refcorner.x;
      dir2.x -= refcorner.x;
      dir1.x /= (lN-1);
      dir2.x /= (lN-1);
    }
       
    for (int ii = 1; ii <= lN-2; ii++) {
      for (int jj = 1; jj <= lN-2; jj++) { 

	dlm_bd->x[isb] = refcorner.x + (double) ii * dir1.x + (double) jj * dir2.x;

	dlm_bd->y[isb] = refcorner.y + (double) ii * dir1.y + (double) jj * dir2.y;

	dlm_bd->z[isb] = refcorner.z + (double) ii * dir1.z + (double) jj * dir2.z;
	isb++;
      }
    }
  }

  int allindextable[8][8] = {{0}};
  int j1,jm1;

  
  /* Add points on the edges without the corners*/
  for (int i = 0; i < nfaces; i++) {
    npoints = gcp->numPointsOnFaces[i];
    i1 = gcp->cornersIndex[i][1];

    for (int j = 1; j < npoints; j++) {
      jm1 = gcp->cornersIndex[i][j-1];
      j1 = gcp->cornersIndex[i][j];
      if (jm1 > j1) {
	if (allindextable[jm1][j1] == 0) {
	  coord c1 = {gcp->cornersCoord[jm1][0], gcp->cornersCoord[jm1][1], gcp->cornersCoord[jm1][2]};
	  coord c2 = {gcp->cornersCoord[j1][0], gcp->cornersCoord[j1][1], gcp->cornersCoord[j1][2]};
	  distribute_points_edge_Cube_v2 (c1, c2, dlm_bd, lN, isb);
	  allindextable[jm1][j1] = 1;
	  isb +=lN-2;
	}
      }
      
      else {
	if (allindextable[j1][jm1] == 0) {
	  coord c1 = {gcp->cornersCoord[j1][0], gcp->cornersCoord[j1][1], gcp->cornersCoord[j1][2]};
	  coord c2 = {gcp->cornersCoord[jm1][0], gcp->cornersCoord[jm1][1], gcp->cornersCoord[jm1][2]};
	  distribute_points_edge_Cube_v2 (c1, c2, dlm_bd, lN, isb);
	  allindextable[j1][jm1] = 1;
	  isb +=lN-2;
	}
      }
      /* printf ("isb after face %d = %d\n", i, isb); */
    }   
  }
  /* printf ("isb after all edges = %d\n", isb); */

  /* Add the final 8 corners points */
  for (int i = 0; i  < gcp->ncorners; i++) {
    dlm_bd->x[isb] = gcp->cornersCoord[i][0];
    dlm_bd->y[isb] = gcp->cornersCoord[i][1];
    dlm_bd->z[isb] = gcp->cornersCoord[i][2];
    isb++;
  }
  /* printf ("isb after 8 corners = %d\n", isb); */
}

void create_FD_Interior_Cube_v2 (particle *p, vector Index_lambda, vector pshift) {

  Cache * c;
  /** Create the cache of the interior points for a cube*/
  c = &(p->Interior);

  /* compute the 3 principal vector of the cube */
  compute_principal_vectors_Cubes (p);

  /* a point x with coord (x,y,z) lies in the cube if the 3
     following conditions are satisfied */
  /* 1- u.p0 <= u.x <= u.p3 */
  /* 2- v.p0 <= v.x <= v.p1 */
  /* 3- w.p0 <= w.w <= w.p4  */
  coord checkpt;
  coord u1 = p->g.u1;
  coord v1 = p->g.v1;
  coord w1 = p->g.w1;
  coord mins = p->g.mins;
  coord maxs = p->g.maxs;

  /* Min/Max coordinates for the AABB (Axed-Aligned-Bounding-Box) */
  coord mincoord = {HUGE, HUGE, HUGE};
  coord maxcoord = {-HUGE, -HUGE, -HUGE};

  double ** table = p->g.cornersCoord;
  for (int ii = 0; ii < p->g.ncorners; ii++) {
    if (mincoord.x > table[ii][0])
      mincoord.x = table[ii][0];

    if (mincoord.y > table[ii][1])
      mincoord.y = table[ii][1];

    if (mincoord.z > table[ii][2])
      mincoord.z = table[ii][2];

    if (maxcoord.x < table[ii][0])
      maxcoord.x = table[ii][0];

    if (maxcoord.y < table[ii][1])
      maxcoord.y = table[ii][1];

    if (maxcoord.z < table[ii][2])
      maxcoord.z = table[ii][2];
  }
  /* if (pid() == 0) { */
  /*   printf ("mincoord = (%f,%f,%f)\n",mincoord.x, mincoord.y, mincoord.z); */
  /*   printf ("maxcoord = (%f,%f,%f)\n",maxcoord.x, maxcoord.y, maxcoord.z); */
  /* } */

  foreach() {
    checkpt.x = x;
    checkpt.y = y;
    checkpt.z = z;

    /* Check only if the point is in the AABB (Axed-Aligned-Bounding-Box) */
    if ((x > mincoord.x) && (x < maxcoord.x)) 
      if ((y > mincoord.y) && (y < maxcoord.y))
	if ((z > mincoord.z) && (z < maxcoord.z))

	  /* If yes: check if it is inside the cube now */
    	  if (is_it_in_cube_v2 (&u1, &v1, &w1, &mins, &maxs, &checkpt)) {
	    cache_append (c, point, 0);
	    /* tagg cell with the number of the particle */
	    if ((int)Index_lambda.y[] == -1)
	      Index_lambda.y[] = p->pnum;
	  }

    /* /\* For periodic boundaries, we have to check for the image-cubes as well *\/ */
    /* if (Period.x)  */
    /*   if (mincoord.x < 0 + X0) */
	
  }

  /* foreach_dimension() { */
  /*   if (Period.x) { */
      
      
  /*   } */
  /* } */
 
  cache_shrink (c);  
  
}

void free_particles (particle * pp, const int n) {
  for (int k = 0; k < n; k++) {
    
    SolidBodyBoundary * sbm = &(pp[k].s);
    free_SolidBodyBoundary(sbm);
    Cache * c = &(pp[k].Interior);
    free(c->p);
    c = &(pp[k].reduced_domain);
    free(c->p);

    /* free corner coordinates tables */
    double * cc;
    GeomParameter * gg = &(pp[k].g);
    for (int i = 0; i < gg->ncorners; i++) {
      cc = &(gg->cornersCoord[i][0]);
      if(cc) {
	free(cc);
      }
    }
    free(gg->cornersCoord);

    /* free corner's indices tables */
    for (int i = 0; i < gg->allFaces; i++) {
      free(gg->cornersIndex[i]);
    }
    free(gg->cornersIndex);
    free(gg->numPointsOnFaces);
  }
}

/** Function that creates the scalar field where the cells contain the
    indices of the dlmfd points. If there is no dlmfd point it is
    tagged -1. */
void create_index_lambda_scalar (const SolidBodyBoundary dlm_bd, vector Index_lambda, const int kk) {
  
  Point lpoint;
  int i;
  Cache *fdlocal;
   
  fdlocal = (Cache *){calloc(dlm_bd.m, sizeof(Cache))};
  
  for (i = 0; i < dlm_bd.m; i++) {

#if dimension == 2 
    lpoint = locate(dlm_bd.x[i], dlm_bd.y[i]);
#elif dimension == 3
    lpoint = locate(dlm_bd.x[i], dlm_bd.y[i], dlm_bd.z[i]);
#endif    
 /*  */
    if ((lpoint.level)  == depth()) {
	/* Create a cache for each fictitious domain's boundary point */
	cache_append(&fdlocal[i], lpoint, 0);
	
      	foreach_cache(fdlocal[i]) {
	  /* Tag cell only if it was not tagged by another particle */
	  if (Index_lambda.x[] < 0){
	    Index_lambda.x[] = i;
	    Index_lambda.y[] = kk;
	    if (level != depth()) {
	      printf("On thread %d, point dlmfd %d at (%f, %f, %f) is in a cell that has not the maximum level of refinnement %d, it is on level %d \n", pid(), i, x, y, z, depth(), level);
	    }
	  }
	}
    }
  }
   
  for (i = 0; i < dlm_bd.m; i++) {
    free(fdlocal[i].p);
  }
  
  free (fdlocal);
  
  boundary ((scalar*) {Index_lambda});

}

#include "weights_functions_backup.h"

# ifndef STENCIL_EXTERIOR
# define STENCIL_EXTERIOR 0
# endif

# ifndef STENCIL_INTERIOR
# define STENCIL_INTERIOR 0
# endif

/** Function that returns the weight associated to a Lagrange
   multipliers for a cell within a foreach() loop (ie for a "real" or
   local cell) associated to a Lagrange multiplier in it's
   neighborhood. The neighboor can be in another process (i.e be a
   ghost cell for the current thread). */
double reversed_weight (particle * pp, const coord weightcellpos, const coord lambdacellpos, const coord lambdapos, const double delta, const coord pshift, const double a, const double b) {
  /* this function has to be embedded within a double foreach_cache() and foreach_neighbor() loop  */
  coord rel = {0., 0., 0.};
  coord relnl = {0., 0., 0.};
  int NCX = 0, CX = 0, weight_id = 0;
  size_t goflag = 0;
  double weight = 0.;
  GeomParameter gcb = pp->g; 
  GeomParameter gcbdum = gcb;
  coord lambdaposdum = lambdapos;
  
  /* Shift artificially the center of the particle and position of the
     Lagrange multipliers for periodic boundary cases */
  foreach_dimension() {
    gcbdum.center.x = gcb.center.x + a*pshift.x;
    lambdaposdum.x = lambdapos.x + b*pshift.x;
  }

  /* compute relative vector X_boundary - X_local = rel from the cell (containning the boundary) position to the boundary's (analytical) position .i.e vector goes from first argument to the second */
  compute_relative_vector (lambdacellpos, lambdaposdum, &rel);

  /* reset dials integers */
  NCX = 0; CX = 0; weight_id = 0; goflag = 0;

  /* assign first normal (gives the quadrant) */ 
  assign_dial (rel, &CX);
  
  /* assign fictitious-boundary's normal (use boundary's analytical position) */
  assign_dial_fd_boundary (pp, lambdaposdum, gcbdum, delta, &NCX);

  if (pp->iswall) {
#if STENCIL_EXTERIOR
    /* Stencils oriented toward the wall */
    if (lambdapos.y > 0) {
      if (lambdapos.x < L0/2)
	NCX = 1;
      if (lambdapos.x > L0/2)
	NCX = 2;
    }
	
    if (lambdapos.y < 0) {
      if (lambdapos.x < L0/2)
	NCX = 4;
      if (lambdapos.x > L0/2)
	NCX = 3;
    }
#endif
#if STENCIL_INTERIOR
    /* Second normal oriented toward the fluid */
    if (lambdapos.y > 0) {
      if (lambdapos.x < L0/2) {
	NCX = 4;
      }
      if (lambdapos.x > L0/2) {
	NCX = 3;
      }
    }
	
    if (lambdapos.y < 0) {
      if (lambdapos.x < L0/2) {
	NCX = 1;
      }
      if (lambdapos.x > L0/2) {
	NCX = 2;
      }
    }
#endif
  }
  
  /* compute relative vector from the neighbor cell to the local cell */
  compute_relative_vector (weightcellpos, lambdacellpos , &relnl);

  /* Assign weight ids. */
  assign_weight_id_quad_outward (NCX, CX, relnl, delta, &weight_id, &goflag);
  
  
  if (goflag == 1) {
    weight = compute_weight_Quad (weight_id, lambdaposdum, lambdacellpos, NCX, CX, delta);
  }
    return weight;
}

void remove_too_close_multipliers (particle * p, vector index_lambda) {

    
  for (int k = 0; k < NPARTICLES; k++) {

    SolidBodyBoundary dlm_lambda_to_desactivate;
    int allocated = 1;
    allocate_SolidBodyBoundary (&dlm_lambda_to_desactivate, allocated*BSIZE);
    int countalloc = 0;
    int other_part;
    int is_in_other_particle;
    foreach_level(depth()) {
      
      int direct_neigh = 0;  is_in_other_particle = 0; other_part = -1;
      if (((int)index_lambda.x[] > -1)  && level == depth() && is_leaf(cell) && (p[k].pnum == (int)index_lambda.y[])) {

	/* check here if this multiplier is not in another particle's
	   domain, if yes desactive it */
	
	for (int l = 0; l < NPARTICLES; l++) {
	  if (l != p[k].pnum) {
      	    particle * other_particle = &(p[l]);
      	    /* printf("thread %d, this is particle %zu checkin on particle %zu iteration %d\n", pid(), p[k].pnum, other_particle->pnum, l); */

      	    /* Check particle's type */
      	    if (other_particle->iscube) {
      	      compute_principal_vectors_Cubes (other_particle);
      	      coord u = other_particle->g.u1;
      	      coord v = other_particle->g.v1;
      	      coord w = other_particle->g.w1;
      	      coord mins = other_particle->g.mins;
      	      coord maxs = other_particle->g.maxs;

      	      /* current cell's position */
      	      coord checkpt = {x, y, z};
      	      if (is_it_in_cube_v2 (&u, &v, &w, &mins, &maxs, &checkpt)) {
      		is_in_other_particle = 1; other_part =  other_particle->pnum;
		break;
      	      }
      	    }
      	    /* if sphere */
      	    else {
      	      GeomParameter gp = other_particle->g;
      	      if (is_it_in_sphere (x, y, z, gp)) {
      		is_in_other_particle = 1; other_part =  other_particle->pnum;
		break;
      	      }
      	    }
      	  }
	  
	}

	if (is_in_other_particle) {
	  /* printf("thread %d, this is particle %zu which has a cell in particle %d\n", pid(), p[k].pnum, other_part); */
	  index_lambda.x[] = -1; index_lambda.y[] = other_part;
	}
	
	/* Check if two (or more) Lagrange multipliers (from different
	   particles) are in the same neighborhoods (i.e a 5x5 stencil
	   in 2D), if yes desactivate them. Otherwise the cells of
	   their stencils may be located in each other's domain */
	direct_neigh = 0;
	foreach_neighbor() {
	  if (((int)index_lambda.x[] > -1)  && level == depth() && is_leaf(cell) && (p[k].pnum != (int)index_lambda.y[])) {

	    direct_neigh = 1;
	    /* We may catch and identical multiplier more than once here */
	    /* Add the multiplier(s) indices to a list because MPI imposes that we can't modify the cell within a foreach_neighbor loop */
	    if (countalloc >= allocated*BSIZE) {
	      allocated ++;
	      reallocate_SolidBodyBoundary (&dlm_lambda_to_desactivate, allocated*BSIZE);
	    }
	    dlm_lambda_to_desactivate.x[countalloc] = x;
	    dlm_lambda_to_desactivate.y[countalloc] = y;
#if dimension == 3
	    dlm_lambda_to_desactivate.z[countalloc] = z;
#endif
	    countalloc ++;
	  }
	}
	/* Desactivate the local multiplier here  */
	if (direct_neigh) {
	index_lambda.x[] = -1; index_lambda.y[] = -1;
	}
      }
    }
    /* if (countalloc > 0) */
    /*   printf("there is %d points to desactivate on thread %d\n", countalloc, pid()); */

    Cache c = {0};
    Point lpoint;
    
#if _MPI
    int size = npe();
    int counts[size];
    
    /* Each thread tells the root how many multipliers it holds  */
    MPI_Gather(&countalloc, 1, MPI_INT, &counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* if (pid() == 0) */
    /*   for (int i = 0; i < size; i++) */
    /* 	printf("particle %d, this is root receiving %d elements by thread %d\n",k, counts[i], i); */

    /* Displacements in the receive buffer for MPI_GATHERV */
    int disps[size];
     /* Displacement for the first chunk of data - 0 */
    for (int i = 0; i < size; i++)
      disps[i] = (i > 0) ? (disps[i-1] + counts[i-1]) : 0;

    double * alldatax = NULL;
    double * alldatay = NULL;
#if dimension == 3 
    double * alldataz = NULL;
#endif

    int m = 0;
    /* Threads 0 allocates and compute the total number of multipliers */
    if (pid() == 0) {
      /* disps[size-1] + counts[size-1] == total number of elements */
      /* printf("thread %d : disps[size-1]+counts[size-1]  = %d\n", pid(), disps[size-1]+counts[size-1]); */
      m = disps[size-1] + counts[size-1];
      alldatax = (double*) calloc (m, sizeof(double));
      alldatay = (double*) calloc (m, sizeof(double));
#if dimension == 3
      alldataz = (double*) calloc (m, sizeof(double));
#endif
    }

    /* Gather on thread 0 the coordinates of the multipliers */ 
    MPI_Gatherv (&dlm_lambda_to_desactivate.x[0], countalloc, MPI_DOUBLE, &alldatax[0], counts, disps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv (&dlm_lambda_to_desactivate.y[0], countalloc, MPI_DOUBLE, &alldatay[0], counts, disps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if dimension == 3 
    MPI_Gatherv (&dlm_lambda_to_desactivate.z[0], countalloc, MPI_DOUBLE, &alldataz[0], counts, disps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    /* send the total number of multipliers to all threads */
    mpi_all_reduce (m, MPI_INT, MPI_MAX);

    /* allocate now alldatax, alldatay, alldataz on threads !=0 */
    if (pid() != 0) {
      alldatax = (double*) calloc (m, sizeof(double));
      alldatay = (double*) calloc (m, sizeof(double));
#if dimension == 3
      alldataz = (double*) calloc (m, sizeof(double));
#endif
    }

  
    /* Now broadcast alldatax,alldatay,alldataz from thread 0 to other threads */
    MPI_Bcast (alldatax, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (alldatay, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if dimension == 3 
    MPI_Bcast (alldataz, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
/*     for (int i = 0; i < m; i++) { */
/* #if dimension == 2  */
/*       printf("thread %d: particle %d, coord of the multiplier %d to be removed (%g,%g,%g)\n", pid(), k, i, alldatax[i], alldatay[i], 0.); */
/* #elif dimension ==3 */
/*       printf("thread %d: particle %d, coord of the multiplier %d to be removed (%g,%g,%g)\n", pid(), k, i, alldatax[i], alldatay[i], alldataz[i]); */
/* #endif */
/*     } */
    
    for (int i = 0; i < m; i++) {
#if dimension == 2 
      lpoint = locate(alldatax[i], alldatay[i]); 
#elif dimension == 3
      lpoint = locate(alldatax[i], alldatay[i], alldataz[i]); 
#endif
      if (lpoint.level > -1)  cache_append(&c, lpoint, 0);
    }

    free(alldatax);
    free(alldatay);
#if dimension == 3 
    free(alldataz);
#endif

#elif _MPI == 0

    for (int i = 0; i < countalloc; i++) {
#if dimension == 2 
      lpoint = locate(dlm_lambda_to_desactivate.x[i], dlm_lambda_to_desactivate.y[i]); 
#elif dimension == 3
      lpoint = locate(dlm_lambda_to_desactivate.x[i], dlm_lambda_to_desactivate.y[i], dlm_lambda_to_desactivate.z[i]); 
#endif
      if (lpoint.level > -1)  cache_append(&c, lpoint, 0);
    }

#endif   /* End if _MPI */

    /* Finally desactivate the multiplicators */
    int counter = 0;
    foreach_cache(c) {
      if (index_lambda.x[] > -1 && index_lambda.y[] > -1) {
	index_lambda.x[] = -1;
	index_lambda.y[] = -1;
	counter ++;
      }
    }
    
    /* printf("thread %d: removed %d multipliers \n", pid(), counter); */

    free(c.p);
    free_SolidBodyBoundary(&dlm_lambda_to_desactivate);
  } /* End loop NPARTICLES */

  boundary((scalar*) {index_lambda});
}

/** Function that tags cells associated to a Lagrange multiplier on
    the boundary of the fictitious domain. */
void reverse_fill_flagfield (particle * p, scalar f, vector index_lambda, const int cacheflag) {

  for (int k = 0; k < NPARTICLES; k++) {
  
    coord rel = {0., 0., 0.};
    coord relnl = {0., 0., 0.};
    int NCX, CX, weight_id;
    size_t goflag = 0;
    coord lambdacellpos = {0., 0., 0.};
    coord lambdapos = {0., 0., 0.};
    coord localcellpos = {0., 0., 0.};
 
    SolidBodyBoundary dlm_lambda = p[k].s;
    GeomParameter gcb = p[k].g;
    Cache * reduced_domain = &(p[k].reduced_domain);

    GeomParameter gcbdum = gcb;
    coord lambdaposdum = lambdapos;
    
    foreach_level(depth()) {
      localcellpos.x = x;
      localcellpos.y = y;
#if dimension == 3 
      localcellpos.z = z;
#endif
    
      goflag = 0;       
      // Check if there is a Lagrange multiplier in the neigborhood of this cell;
      foreach_neighbor() {
	if (((int)index_lambda.x[] > -1)  && level == depth() && is_leaf(cell) && (p[k].pnum == (int)index_lambda.y[])) {
	 
	  lambdacellpos.x = x; lambdacellpos.y = y; 
	  lambdapos.x = dlm_lambda.x[(int)index_lambda.x[]];
	  lambdapos.y = dlm_lambda.y[(int)index_lambda.x[]];
	
#if dimension == 3
	  lambdacellpos.z = z;
	  lambdapos.z = dlm_lambda.z[(int)index_lambda.x[]];
#endif

	  /* Shift artificially the center of the particle and position of the
	     Lagrange multipliers for periodic boundary cases */
	  coord ppshift = {0., 0., 0.};
	  foreach_dimension() {
	    ppshift.x = 0.;
	    if (lambdacellpos.x < 0.) ppshift.x = -L0;
	    if (lambdacellpos.x > L0) ppshift.x = +L0;
	    gcbdum.center.x = gcb.center.x + ppshift.x;
	    lambdaposdum.x = lambdapos.x + ppshift.x;
	  }

	  
	  /* compute relative vector X_boundary - X_local = rel from the cell (containning the boundary) position to the boundary's (analytical) position */
	  compute_relative_vector (lambdacellpos, lambdaposdum, &rel);

	  /* reset dials integers */
	  NCX = 0; CX = 0; weight_id = 0; 

	  /* assign first normal (gives the quadrant) */ 
	  assign_dial (rel, &CX);

	  
	  /* assign fictitious-boundary's normal (use boundary's position) */
	  assign_dial_fd_boundary (&p[k], lambdaposdum, gcbdum, Delta, &NCX);

	  if (p->iswall) {
		  
#if STENCIL_EXTERIOR
	    /* Stencils oriented toward the wall */
	    if (lambdapos.y > 0) {
	      if (lambdapos.x < L0/2)
		NCX = 1;
	      if (lambdapos.x > L0/2)
		NCX = 2;
	    }
	
	    if (lambdapos.y <0) {
	      if (lambdapos.x < L0/2)
		NCX = 4;
	      if (lambdapos.x > L0/2)
		NCX = 3;
	    }
#endif
	
#if STENCIL_INTERIOR
	    /* Stencils oriented toward the fluid */
	    if ((lambdapos.y) > 0) {
	      if (lambdapos.x < L0/2)
		NCX = 4;
	      if (lambdapos.x > L0/2)
		NCX = 3;
	    }
	
	    if (lambdapos.y <0) {
	      if (lambdapos.x < L0/2) {
		NCX = 1;
	      }
	  
	      if (lambdapos.x > L0/2)
		NCX = 2;
	    }
#endif
	    
	  }
	
	  /* compute relative vector neighbor to the local cell */
	  compute_relative_vector (localcellpos, lambdacellpos , &relnl);
	
	  /* Assign weight ids. */
	  assign_weight_id_quad_outward (NCX, CX, relnl, Delta, &weight_id, &goflag);
	} // end if (lambda.x[] > -1)
      } // end foreach_neigboor loop

      /* If the cell belongs to the stencil of a Lagrange multiplier tagg it and create a reduced domain to optimize the stencil-traversal cost */
      if (goflag == 1) {
	f[] = 1;
	if (cacheflag == 1)
	  cache_append (reduced_domain, point, 0);
      }
    }
  
    boundary ({f, index_lambda.x, index_lambda.y});
  }// loop on particles id 
}

/** Function that initialize the particles and the scalar/vector
fields needed to the method */
void allocate_and_init_particles (particle * p, const int n, vector e, scalar g, scalar h, vector pshift) {
  
  foreach() {
    e.x[] = -1;// index_lambda.x => lambda indices
    e.y[] = -1;// particle structure/fd number;
    e.z[] = -1;// constrained cells (-1:not constrained)
    g[] = 0;   // flagfield
    h[] = 0;   // flagfield_mailleur

    /* Reset the shift vector for the periodic case */
    if (Period.x || Period.y || Period.z) { 
      foreach_dimension() {
	pshift.x[] = 0.;
      }
    }
  }
  
  boundary ({g, h});
  boundary ((scalar *){e, pshift});


  for (int k = 0; k < n; k++) {
    Cache * c = NULL;

#if debugBD == 0
    GeomParameter gci = p[k].g;
    int m = 0;
  
    
    /* default case, sphere */
    if (!(p[k].iswall) && !(p[k].iscube)) {
      compute_nboundary (gci, &m);
      allocate_SolidBodyBoundary(&(p[k].s), m);
      create_FD_Boundary_Particle (gci, &(p[k].s), m, pshift);
    }
    
    if (p[k].iswall) {
      coord wallpos = p[k].wallmax;
      compute_nboundary_Wall (wallpos, &m);
      allocate_SolidBodyBoundary(&(p[k].s), m);
      create_FD_Boundary_Wall (&p[k]);
    }

    if (p[k].iscube) {
      int lN = 0;        
      compute_nboundary_Cube_v2 (&gci, &m, &lN);
      allocate_SolidBodyBoundary(&(p[k].s), m);
      create_FD_Boundary_Cube_v2 (&gci, &(p[k].s), m, lN, pshift);
    }

    create_index_lambda_scalar ((p[k].s), e, k);
    c = &(p[k].reduced_domain);
    allocate_Cache(c);
#endif     
#if debugInterior == 0
    c = &(p[k].Interior);
    allocate_Cache(c);
#endif
  }

}

#if DLM_Moving_particle
#if TRANSLATION
void compute_contact_distance (particle * p, const coord gci, double * delta) {
  GeomParameter gc = p->g;
  double radius = gc.radius;
  *delta = gci.y - radius;
}

double compute_fwo (const double en, const double v, const double wo) {

  double k = log(en)/sqrt(sq(pi) + sq(log(en)));
  double g = atan(-sqrt(1 - sq(k))/k);
  double val;
  return val = v*exp(k*g/sqrt(1- sq(k)))*sin(g)/(wo*sqrt(1 - sq(k)));

}

double derivative_fwo (const double en, const double v, const double wo) {

  double k = log(en)/sqrt(sq(pi) + sq(log(en)));
  double g = atan(-sqrt(1 - sq(k))/k);
  double val;
  return val = -v*exp(k*g/sqrt(1 - sq(k)))*sin(g)/(sqrt(1 - sq(k))*sq(wo)); 
}

void compute_wo (particle * p) {
  
  /* compute an estimate of (kn,gamman) derived from (deltamax,en) */

  /* Guess value for wo */
  double wo = 100.;
  double epsilo = 1e-2;
  int maxiter = 200;
  double aa = 0.;
  double bb = 0.;
  GeomParameter g = p->g;
  double r = g.radius;
  double deltamax = (p->wished_ratio)*r;
  double en = p->en;
  double v = p->vzero;
    
  /* Newton iteration to find the root of delta(w0) */
  if (en < 1) {
    for (int pp = 1; pp <= maxiter; pp++) {

      aa = compute_fwo (wo, v, en) - deltamax;
      bb = derivative_fwo (wo, v, en);

      if (abs(aa-deltamax)/deltamax < epsilo) {
	break;
	  }
      wo += -aa/bb;
    }
    p->kn = (p->M)*sq(wo);
  }
  else {
    /* if en = 1 the formula simplifies as such */
    p->kn = (p->M)*sq(p->vzero)/(sq(deltamax));
  }
  
}

void compute_Fontact (coord * Fc, particle * p, coord * gci, coord * U, const double gamman) {

  double delta_colision = 0.;

  foreach_dimension() {
    (*Fc).x = 0.;
  }
  
  compute_contact_distance(p, *gci, &delta_colision);
    
  if (delta_colision < 0.) {
    double kn = p->kn;
    double M = p->M;
    coord vrel = *U;
    coord Fel = {0., 0., 0.};
    coord Fdm = {0., 0., 0.};
    coord normalvec = p->normalvector;

    foreach_dimension() {
      /* compute Hookean elastic restoring force */
      Fel.x = -kn*delta_colision*normalvec.x;
    
      /* compute viscous dynamic force */
      Fdm.x = -2.*gamman*M*fabs(vrel.x);
    
      /* Fc is the sum of these two forces */
      (*Fc).x = Fel.x + Fdm.x;
    }
  }
}

void granular_subproblem (particle * p, const int gravity_flag, const double dt, const double rho_f) {
  // Mini Granular solver, which solves (1 - rho_f/rho_s)MdU/dt = (1-rho_f/rho_s)Mg + F_c
  // (1 - rho_f/rho_s) Ip dw/dt = sum_all_particles (r x F_c)  - (1 - rho_f/rho_s)w x Ip w
  // where F_c is the contact force, g the gravity acceleration and M the particle's mass and Ip the inertia tensor

  double wo = 0. ;
  double gamman = 0.;
  
  double T_c = 0., dtg = 0., M = 0., kn = 0., en = 0.;
  
  int miter, gi;
  coord Uold, Xold;
  coord k1, k2, k3, k4, Xtemp, Utemp;
  double fsf;

  /* particle's structure pointers  */
  GeomParameter * gci;
  GeomParameter * gcinm1;
  coord * U;
  coord * Unm1;
 
  coord * gravity;
  coord decal = {X0, Y0, Z0};
  
  for (int k = 0; k < NPARTICLES; k++) {
    fsf = (1. - (rho_f)/(p[k].rho_s));
   
    miter = 0;
    M = p[k].M;
    kn = p[k].kn;
    en = p[k].en;
    
    
    /* compute wo and gamman */
    wo = sqrt(kn/M); // wo = sqrt(2kn/M) for sphere-sphere contact and wo = sqrt(kn/M) for sphere-wall contact (Powder tech. Rakotonirina 2018)
    gamman = -wo*log(en)/sqrt(sq(pi) + sq(log(en)));

    /* compute the contact time Tc */
    T_c = pi/(sqrt(sq(wo) - sq(gamman)));

    /* set granular timestep */
    dtg = T_c/20.;
    miter = ceil(dt/dtg);
    dtg = dt/miter;
    
    if (k == 0) fprintf (stderr,"Tc = %g, dt = %20.18f, dtg = %20.18f, Tc/dtg = %g, miter = %d, particle = %d\n",T_c, dt, dtg, T_c/dtg, miter, k);
    
    /* get particle structure pointers */

    /* position */
    gcinm1 = &(p[k].gnm1);
    gci = &(p[k].g);

    /* translational velocity */
    U = &(p[k].U);
    Unm1 = &(p[k].Unm1);

    gravity = &(p[k].gravity);
    
    /* fprintf (stderr,"gravity = (%f,%f,%f)\n",(*gravity).x,(*gravity).y,(*gravity).z); */
    
    if (gravity_flag) {
      /* Before solving the granular problem: save previous particle's
	 position and velocity (predictor step
	 with gravity, first subproblem after N-S) */
      
      foreach_dimension() {
	(*Unm1).x = (*U).x;

	
	/* Check if the domain is periodic, if yes shift the particle's
	   position when the end of the domain is reached */
	if (Period.x) {
	  if ((*gci).center.x > (L0 + decal.x)) {
	    (*gci).center.x -= L0;
	  }
	  if ((*gci).center.x < (0. + decal.x)) {
	    (*gci).center.x += L0;
	  }
	}

	
	(*gcinm1).center.x = (*gci).center.x;
      }
    }
    
    else {
      /* corrector step without gravity: fourth subproblem after the
	 fictitious domain problem (correction step) */
      foreach_dimension() { 
	(*gci).center.x = (*gcinm1).center.x;
      }
    }
    

    /*  integrating in time with rk4 */
    for (gi = 1; gi <= miter; gi++) {
      Uold = (*U);
      Xold = (*gci).center;

      /* compute k1 */
      compute_Fontact (&k1, &p[k], &Xold, &Uold, gamman);
      foreach_dimension() {
	k1.x /= (fsf*M);
	k1.x += gravity_flag*(*gravity).x;
	Xtemp.x = Xold.x + 0.5*dtg*Uold.x;
	Utemp.x = Uold.x + 0.5*dtg*k1.x;
      }
	
      /* compute k2 */
      compute_Fontact (&k2, &p[k], &Xtemp, &Utemp, gamman);
      foreach_dimension() {
	k2.x /= (fsf*M);
	k2.x += gravity_flag*(*gravity).x;
	Xtemp.x = Xold.x + 0.5*dtg*Uold.x + 0.25*sq(dtg)*k1.x;
	Utemp.x = Uold.x + 0.5*dtg*k2.x;
      }
	
      /* compute k3 */
      compute_Fontact (&k3, &p[k], &Xtemp, &Utemp, gamman);
      foreach_dimension() {
	k3.x /= (fsf*M);
	k3.x += gravity_flag*(*gravity).x;
	Xtemp.x = Xold.x + dtg*Uold.x + 0.5*sq(dtg)*k2.x;
	Utemp.x = Uold.x + dtg*k3.x;
      }
	
      /* compute k4 */
      compute_Fontact (&k4, &p[k], &Xtemp, &Utemp, gamman);
      foreach_dimension() {
	k4.x /= (fsf*M);
	k4.x += gravity_flag*(*gravity).x;
      }

      
      foreach_dimension () {
	(*gci).center.x = Xold.x + dtg*Uold.x + sq(dtg)*(k1.x + k2.x + k3.x)/6.;
	
	/* Check if the domain is periodic, if yes shift the particle's
	   position when the end of the domain is reached */
	if (Period.x) {
	  if ((*gci).center.x > (L0 + decal.x)) {
	    (*gci).center.x -= L0;
	  }
	  if ((*gci).center.x < (0. + decal.x)) {
	    (*gci).center.x += L0;
	  }
	}



	(*U).x = Uold.x + dtg*(k1.x + 2*k2.x + 2*k3.x + k4.x)/6.;
      }
    }
    if (k == 0) {
      if (gravity_flag) {
	fprintf (stderr,"Prediction: particle-0's velocity on thread %d is (%20.18f, %20.18f, %20.18f)\n",pid(), (*U).x, (*U).y, (*U).z);
	fprintf (stderr,"Prediction: particle-0's position on thread %d is (%20.18f, %20.18f, %20.18f)\n",pid(), (*gci).center.x, (*gci).center.y, (*gci).center.z);
      
      }
      else {
	fprintf(stderr,"Correction: particle-0's velocity on thread %d is (%20.18f, %20.18f, %20.18f)\n",pid(), (*U).x, (*U).y, (*U).z);
	fprintf(stderr,"Correction: particle-0's position on thread %d is (%20.18f, %20.18f, %20.18f)\n",pid(), (*gci).center.x, (*gci).center.y, (*gci).center.z);
      }
    }
  }
}
#endif
#endif

/** Diagnostic functions. */
void writer_headers (FILE * pdata, FILE * sl) {
#if _MPI
  if(pid() == 0)
#endif
    {
      /* Write header for particles positions and velocities */
#if dimension == 2
      fprintf(pdata, "#t\t position.x\t position.y\t U.x\t U.y\t w.z\n");
      fprintf(sl, "#t\t F.x\t F.y\t T.z\n");
      
#elif dimension == 3
      fprintf (pdata, "#t\t position.x\t position.y\t position.z\t U.x\t U.y\t U.z\t w.x\t w.y\t w.z\n");
      fprintf (sl, "#t\t F.x\t F.y\t F.z\t T.x\t T.y\t T.z\n");
#endif
      fflush(pdata);
      fflush(sl);
    }
}

void particle_data (particle * p, const double t, const int i, FILE ** pdata) {
  
  for (int k = 0; k < NPARTICLES; k++) {
    GeomParameter * GCi = &(p[k].g);
#if DLM_Moving_particle
#if TRANSLATION 
    coord * U = &(p[k].U);
#endif
#if ROTATION
    coord * w =  &(p[k].w);
#endif
#endif
#if dimension == 2
    if(pid() == 0) {
      fprintf(pdata[k], "%20.18f\t %20.18f\t %20.18f %20.18f\t %20.18f\t %20.18f \n", t, (*GCi).center.x, (*GCi).center.y
#if DLM_Moving_particle
#if TRANSLATION
	      , (*U).x, (*U).y
#elif TRANSLATION == 0
	      , 0., 0.
#endif
#if ROTATION
	      , (*w).z
#elif ROTATION == 0
	      , 0.
#endif
	      
#elif DLM_Moving_particle == 0
	      , 0., 0.
	      , 0.
#endif
	      );
      fflush(pdata[k]);
    }
#elif dimension == 3
    if(pid() == 0) {
      fprintf(pdata[k], "%20.18f\t %20.18f\t %20.18f\t %20.18f\t %20.18f\t %20.18f %20.18f\t %20.18f\t %20.18f\t %20.18f\n", t, (*GCi).center.x, (*GCi).center.y, (*GCi).center.z

#if DLM_Moving_particle
#if TRANSLATION
	      , (*U).x, (*U).y, (*U).z
#elif TRANSLATION == 0
	       , 0., 0., 0.
#endif	      
#if ROTATION
	      ,(*w).x, (*w).y, (*w).z
#elif ROTATION == 0
	       , 0., 0., 0.
#endif
#elif DLM_Moving_particle == 0
	      , 0., 0., 0.
	      , 0., 0., 0.
#endif
	      );
      fflush(pdata[k]);
    }
#endif
  }
}

coord sumLambda (particle * p, FILE ** sl, const double t, const double dt, scalar Flagfield, vector Lambda, vector Index_lambda, const double rho_f, vector pshift) {
  /* Compute hydrodynamic force & torque and write to a file */
  /* Fh = <lambda,V>_P + (rho_f/rho_s)MdU/dt */
  /* Mh = <lambda,xi^GM>_P + (rho_f/rho_s).( Idom/dt + om x (I.om)) */

  /* The inertia tensor is: */
  /* Ip[0] = Ixx */
  /* Ip[1] = Iyy */
  /* Ip[2] = Izz */
  /* Ip[3] = Ixy */
  /* Ip[4] = Ixz */
  /* Ip[5] = Iyz */  

  coord lambdasumint;
  coord lambdasumboundary;
  coord lambdasum;
  coord crossLambdaSumInt;
  coord crossLambdaSumBoundary;
  coord crossLambdaSum; 
  Cache * Interior[NPARTICLES];
  Cache * Boundary[NPARTICLES];
  GeomParameter * gci;
  SolidBodyBoundary * sbm;
  
  /* Loop over all particles */
  for (int k = 0; k < NPARTICLES; k++) {

    foreach_dimension() {
      lambdasumint.x = 0.;
      lambdasumboundary.x = 0.;
      lambdasum.x = 0;
      
      crossLambdaSumInt.x = 0.;
      crossLambdaSumBoundary.x = 0.;
      crossLambdaSum.x = 0.;
    }

#if dimension ==2 
    lambdasumint.z = 0.;
    lambdasumboundary.z = 0.;
    lambdasum.z = 0;
      
    crossLambdaSumInt.z = 0.;
    crossLambdaSumBoundary.z = 0.;
    crossLambdaSum.z = 0.;
#endif

    /* Interior domain of the particle */
    Interior[k] = &(p[k].Interior);

    /* Boundary domain of the particle */
    Boundary[k] = &(p[k].reduced_domain);
    gci = &(p[k].g);
    sbm = &(p[k].s);
    coord lambdapos;
#if DLM_Moving_particle
    double rho_s = p[k].rho_s;
#if TRANSLATION
    double M = p[k].M;
    coord Unm1 = p[k].Unm1;
    coord U = p[k].U;
#endif
#if ROTATION
    coord wnm1 = p[k].wnm1;
    coord w = p[k].w;
    coord Iom, omIom, Idomdt;
#endif
#endif
    
    // Particule's interior multipliers
    foreach_cache((*Interior[k])) {
      if (Flagfield[] < 1 && k == Index_lambda.y[]) {

	/* Compute Fh = <lambda,V>_P */
	/* For moving particles the additional term +
	   (rho_f/rho_s).M.dU/dt is added after the mpi_calls below */
	
	foreach_dimension() {
	  lambdasumint.x += Lambda.x[];
	}

	/* Compute Mh = <lambda,xi^GM>_P */  
	/* For moving particles, the additional term +
	   (rho_f/rho_s).(I.dom/dt + om ^ (I.om)) is added after mpi
	   calls below */
	
	/* Modify temporerly the particle center position for periodic boundary condition */
	foreach_dimension()
	  (*gci).center.x += pshift.x[];

	crossLambdaSumInt.x += (Lambda.z[]*(y - (*gci).center.y) - Lambda.y[]*(z - (*gci).center.z));
	crossLambdaSumInt.y += (Lambda.x[]*(z - (*gci).center.z) - Lambda.z[]*(x - (*gci).center.x));
	crossLambdaSumInt.z += (Lambda.y[]*(x - (*gci).center.x) - Lambda.x[]*(y - (*gci).center.y));
	foreach_dimension()
	  (*gci).center.x -= pshift.x[];

      }
    }//end foreach_cache()

    // Particle's boundary multipliers
    foreach_cache((*Boundary[k])) {
      if ((Index_lambda.x[] > -1) && (p[k].pnum == (int)Index_lambda.y[])) {
	lambdapos.x = (*sbm).x[(int)Index_lambda.x[]];
	lambdapos.y = (*sbm).y[(int)Index_lambda.x[]];
	lambdapos.z = 0.;
#if dimension == 3
	lambdapos.z = (*sbm).z[(int)Index_lambda.x[]];
#endif

	/* Compute Fh = <lambda,V>_P */
	foreach_dimension() {
	  lambdasumboundary.x += Lambda.x[];
	}


	/* Modify temporerly the particle center position for periodic boundary condition */
	foreach_dimension()
	  (*gci).center.x += pshift.x[];
	
	/* Compute Mh = <lambda,xi^GM>_P */
        crossLambdaSumBoundary.x += (Lambda.z[]*(lambdapos.y - (*gci).center.y) - Lambda.y[]*(lambdapos.z - (*gci).center.z));
        crossLambdaSumBoundary.y += (Lambda.x[]*(lambdapos.z - (*gci).center.z) - Lambda.z[]*(lambdapos.x - (*gci).center.x));
        crossLambdaSumBoundary.z += (Lambda.y[]*(lambdapos.x - (*gci).center.x) - Lambda.x[]*(lambdapos.y - (*gci).center.y));

	foreach_dimension()
	  (*gci).center.x -= pshift.x[];
      }
    }
 
#if _MPI
    foreach_dimension() {
      mpi_all_reduce (lambdasumint.x, MPI_DOUBLE, MPI_SUM);
      mpi_all_reduce (lambdasumboundary.x, MPI_DOUBLE, MPI_SUM);
      mpi_all_reduce (crossLambdaSumInt.x, MPI_DOUBLE, MPI_SUM);
      mpi_all_reduce (crossLambdaSumBoundary.x, MPI_DOUBLE, MPI_SUM);
    }
#if dimension == 2
    mpi_all_reduce (crossLambdaSumInt.z, MPI_DOUBLE, MPI_SUM);
    mpi_all_reduce (crossLambdaSumBoundary.z, MPI_DOUBLE, MPI_SUM);
#endif
#endif

    /* The Force term (rho_f/rho_s)MdU/dt is added here for
       translating particles */
    
#if DLM_Moving_particle
#if TRANSLATION
    foreach_dimension()
      lambdasumint.x += (rho_f/rho_s)*M*(U.x - Unm1.x)/dt;
#endif
    
    /* The Torque term (rho_f/rho_s).(Idom/dt + om ^ (I.om)) is added
	 here for rotating particles */
    
#if ROTATION
    /* I.om term */
    Iom.x = p[k].Ip[0]*w.x + p[k].Ip[3]*w.y + p[k].Ip[4]*w.z;
    Iom.y = p[k].Ip[3]*w.x + p[k].Ip[1]*w.y + p[k].Ip[5]*w.z;
    Iom.z = p[k].Ip[4]*w.x + p[k].Ip[5]*w.y + p[k].Ip[2]*w.z;

    /* om^(I.om) term */
    omIom.x = w.y*Iom.z - w.z*Iom.y;
    omIom.y = w.z*Iom.x - w.x*Iom.z;
    omIom.z = w.x*Iom.y - w.y*Iom.x;

    /* Idom/dt term; */
    Idomdt.x = p[k].Ip[0]*(w.x - wnm1.x) + p[k].Ip[3]*(w.y - wnm1.y) + p[k].Ip[4]*(w.z - wnm1.z);
    Idomdt.y = p[k].Ip[3]*(w.x - wnm1.x) + p[k].Ip[1]*(w.y - wnm1.y) + p[k].Ip[5]*(w.z - wnm1.z);
    Idomdt.z = p[k].Ip[4]*(w.x - wnm1.x) + p[k].Ip[5]*(w.y - wnm1.y) + p[k].Ip[2]*(w.z - wnm1.z);
	
    crossLambdaSumInt.x += (rho_f/rho_s)*((Idomdt.x)/dt + omIom.x);
    crossLambdaSumInt.y += (rho_f/rho_s)*((Idomdt.y)/dt + omIom.y);
    crossLambdaSumInt.z += (rho_f/rho_s)*((Idomdt.z)/dt + omIom.z);
#endif	  
#endif

    
    foreach_dimension() {
      lambdasum.x = lambdasumint.x + lambdasumboundary.x;
      crossLambdaSum.x = crossLambdaSumInt.x + crossLambdaSumBoundary.x;
    }
    

#if dimension == 2
    crossLambdaSum.z = crossLambdaSumInt.z + crossLambdaSumBoundary.z;
      
    if(pid() == 0) {
      fprintf(sl[k], "%f\t %20.18f\t %20.18f\t %20.18f\n", t, lambdasum.x,  lambdasum.y, crossLambdaSum.z);
      fflush(sl[k]);
    }
#elif dimension == 3
    if(pid() == 0) {
      fprintf(sl[k], "%20.18f\t %20.18f\t %20.18f\t %20.18f\t %20.18f\t %20.18f\t %20.18f\n", t, lambdasum.x,  lambdasum.y, lambdasum.z, 
              crossLambdaSum.x, crossLambdaSum.y, crossLambdaSum.z);
      fflush(sl[k]);
    }
#endif
  }
  
  return lambdasum;
}

void init_file_pointers (FILE ** p, FILE ** d, const size_t rflag) {
  char name[800];
  char name2[800];
#if _MPI
  if (pid() == 0)
#endif
    {
    for (int k = 0; k < NPARTICLES; k++) {

      sprintf (name, "particle-data-%d", k);
      sprintf (name2, "sum_lambda-%d", k);

      if (!rflag) {
	p[k] = fopen (name,  "w"); 
	d[k] = fopen (name2, "w");

	/* if (k == (NPARTICLES - 1)) { */
	/*   printf ("pointer sum lambda= %p\n", (void *)d[k]); */
	/*   printf ("pointer particle data = %p\n", (void *)p[k]); */
	/* } */
	/* Write headers in these files */
	writer_headers (p[k],  d[k]);
      }
      else {
	p[k] = fopen (name,  "a");
	d[k] = fopen (name2, "a");
      }
    }
  }
}

void vorticity_3D (const vector u, vector omega) {
  
  foreach()
    foreach_dimension()
      omega.x[] = ( (u.z[0,1,0] - u.z[0,-1,0]) - (u.y[0,0,1] - u.y[0,0,-1]) )/2.*Delta;

    
  
  boundary((scalar *){omega});
}

void compute_flowrate_right (FILE * pf, const vector u, const int level) {

  coord velo = {0., 0., 0.};

  foreach_dimension()
    velo.x = 0.;

#if !adaptive 
  Cache intDomain = {0};
  Point lpoint;
  foreach() {
    if (x > (L0 - Delta)) { 
	lpoint = locate(x, y, z);
	cache_append(&intDomain, lpoint, 0);
      }
  }
		
  cache_shrink(&intDomain);
    
  foreach_cache(intDomain){
    foreach_dimension()
      velo.x += sq(Delta)*u.x[];
  }
  
  free(intDomain.p);
#endif

#if adaptive
  double hh = L0/pow(2,level);
  double xi = 0., yj = 0., zval = 0.;
  int ii = 0, jj = 0;
  
  foreach_level(level) {
    zval = L0 - 0.5*Delta + Z0;
  }

  /* printf("zz = %f, h = %f\n", zval, hh); */
  
  for (ii = 0; ii < pow(2,level); ii++) {
    xi = 0.5*hh + ii*hh + X0;
    for (jj = 0; jj < pow(2,level); jj++) {
      yj = 0.5*hh + jj*hh + Y0;
      coord uu = {0., 0., 0.};
      foreach_dimension() {
	uu.x = interpolate (u.x, xi, yj, zval);
	if (uu.x < nodata)
	  velo.x += sq(hh)*uu.x;
      }
      /* printf("thread %d, u = (%f,%f,%f)\n",pid(), interpolate(u.x, xi, yj, zz), interpolate(u.y, xi, yj, zz), interpolate(u.z, xi, yj, zz)); */
    }
  }  
#endif

  /* printf("thread %d, velo = (%f,%f,%f)\n",pid(), velo.x, velo.y, velo.z); */


#if _MPI
  foreach_dimension ()
    mpi_all_reduce (velo.x, MPI_DOUBLE, MPI_SUM);
#endif
  
  fprintf(pf, "%f\t %20.18f\t %20.18f\t %20.18f\n", t, velo.x, velo.y, velo.z);
  fflush(pf);

}

void pressure_law (scalar pres, FILE * ppf, const coord hv, const double t) {

  double plaw = 0.;
  
  plaw = interpolate (pres, hv.x, hv.y, hv.z);

#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
  if (plaw == nodata)
    plaw = 0.;
  MPI_Barrier(MPI_COMM_WORLD);

  if (plaw != nodata)
    mpi_all_reduce(plaw, MPI_DOUBLE, MPI_SUM);
#endif
 
  if (pid() == 0) {
    fprintf(ppf, "%20.18f\t %20.18f\n", t, plaw);
    fflush(ppf);
  }
}

/** dump/resto particle's structure functions */
void dump_particle (particle * p) {

  FILE * fp;
  char * file = "dump_particle";

  if (file && (fp = fopen (file, "w")) == NULL) {
    perror (file);
    exit (1);
  }
  assert (fp);

  fwrite(p, sizeof(*p), NPARTICLES, fp);
  fclose(fp);
}

bool restore_particle (particle * p) {

  FILE * fp;
  char * file = "dump_particle";
  if (file && (fp = fopen (file, "r")) == NULL)
    return false;
  assert (fp);

  if (fread (p, sizeof(*p), NPARTICLES, fp) < 1) {
    fprintf (ferr, "restore_particle(): error reading particle data\n");
    exit (1);
  }
  fclose(fp);

  return true;
}

void inverse3by3matrix(double oldMatrix[3][3], double ** inversedMatrix) {

  double block1 = oldMatrix[1][1]*oldMatrix[2][2] - oldMatrix[1][2]*oldMatrix[2][1];
  double block2 = oldMatrix[1][2]*oldMatrix[2][0] - oldMatrix[1][0]*oldMatrix[2][2];
  double block3 = oldMatrix[1][0]*oldMatrix[2][1] - oldMatrix[1][1]*oldMatrix[2][0];

  double determinat = oldMatrix[0][0]*block1 + oldMatrix[0][1]*block2 + oldMatrix[0][2]*block3;

  /* m[i][j], *(*(m + i) + j) */
  
  *(*(inversedMatrix + 0) + 0) = block1/determinat;
  *(*(inversedMatrix + 0) + 1) = (oldMatrix[0][2]*oldMatrix[2][1] - oldMatrix[0][1]*oldMatrix[2][2])/determinat;
  *(*(inversedMatrix + 0) + 2) = (oldMatrix[0][1]*oldMatrix[1][2] - oldMatrix[0][2]*oldMatrix[1][1])/determinat; 
  *(*(inversedMatrix + 1) + 0) = block2 / determinat;
  
  *(*(inversedMatrix + 1) + 1) = (oldMatrix[0][0]*oldMatrix[2][2] - oldMatrix[0][2]*oldMatrix[2][0])/determinat;
  *(*(inversedMatrix + 1) + 2) = (oldMatrix[0][2]*oldMatrix[1][0] - oldMatrix[0][0]*oldMatrix[1][2])/determinat;
  *(*(inversedMatrix + 2) + 0) = block3/determinat;
  *(*(inversedMatrix + 2) + 1) = (oldMatrix[0][1]*oldMatrix[2][0] - oldMatrix[0][0]*oldMatrix[2][1])/determinat;
  *(*(inversedMatrix + 2) + 2) = (oldMatrix[0][0]*oldMatrix[1][1] - oldMatrix[0][1]*oldMatrix[1][0])/determinat ;
}

int totalcells() {
  int t = 0;
  foreach() {
    t++;
  }
#if _MPI
  mpi_all_reduce (t, MPI_INTEGER, MPI_SUM);
#endif
  
  return t;
}

int total_dlmfd_cells (particle * p, const int np) {
  int apts = 0;
  for (int k = 0; k < np; k++) {
#if debugInterior == 0
    apts += p[k].Interior.n;
#endif
#if debugBD == 0
    apts += p[k].reduced_domain.n;
#endif
  }
  
#if _MPI
  mpi_all_reduce (apts, MPI_INTEGER, MPI_SUM);
#endif

  return apts;
}

int total_dlmfd_multipliers (particle * p, const int np) {

  int apts = 0;
  for (int k = 0; k < np; k++) {
#if debugInterior == 0
    apts += p[k].Interior.n;
#endif
  }
  
#if _MPI
  mpi_all_reduce (apts, MPI_INTEGER, MPI_SUM);
#endif
  
  for (int k = 0; k < np; k++) {
#if debugBD == 0
    SolidBodyBoundary * bla;
    bla = &(p[k].s);
    apts += bla->m;
#endif
  }

  return apts;
}

double getDelta() {
  /* Return the smallest grid size  */
  /* ----------------------------- */
  int ii = 0;
  double dd = 0.;
  foreach_level(depth()) {
    ii++;
    dd = Delta;
    if (ii > 0)
      break;
  }
#if _MPI
  mpi_all_reduce (dd, MPI_DOUBLE, MPI_MIN);
#endif
  return dd;
}
