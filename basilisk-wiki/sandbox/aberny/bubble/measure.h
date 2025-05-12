/**
## Tools
We first define some tools for the drops (or the bubble) we will measure

The first think will be to declare a new type of structure: Phase. It contains
all the information we want for a droplet (or a bubble in the ROMEO case).*/

#include <stdlib.h>

/**
This library can be used for a 3D simulation or one with the retracting film. If
FILM has not be defined before, we set it to 0*/
#ifndef FILM
#define FILM 0
#endif

typedef struct Phase Phase;
struct Phase {
  int tag; // tag number
  double x; // position
  double velo; // velocity
  double veloAdim; // dimensionless velocity
  double volume; // volume of the drop (or bubble)
  double radius; // radius of the drop (or bubble)
  double reynolds; // Reynolds number of the drop (or bubble)
  double cellsRadius; // equivalent number of cell per radius.
  int cells; // number of cell in the drop (or bubble)
};

/**
We initialize the value in our structure. All is set to 0.*/

void initPhase (Phase * drop) {
  drop->tag = 0;
  drop->x = 0;
  drop->velo=0;
  drop->veloAdim=0;
  drop->volume=0;
  drop->radius=0;
  drop->reynolds=0;
  drop->cellsRadius=0;
  drop->cells=0;
}

/** We define the comparaison function, for the qsort function in stdlib.h. The
    function will sort the results from the biggest to the smallest one. We
    will use the position as a sort argument*/

static int comparePosition  (const void *a, const void *b) {
  const Phase * pa = a, * pb = b;
  return pb->x > pa->x ? 1 : -1;
}


/**
## Physical measures: drops

The post-process event measure 2 types of physical parameters. The first one is
a measure of the jet parameters. It's measuring the velocity and the position 
of the jet, at every time steps. The second one  is detecting the formation of
liquid droplets at every time steps.


The jet measures should only measure the velocity of the jet. It's not a 
problem as long as we have one phase. In the other case, we need to make sure 
that we are in the phase of the jet. For that, we will tag the liquid phase. 
If there is more than one phase, we will measure the phase volume. The jet is 
in the biggest phase.*/

double tprev = -1, xprev = -2, yprev = 0, rprev = 0; // initialisation of parameters to some "stupid" value
double tCross = 1., tcompt = 0;

int phasePrev = 1;
int compt = 0;
bool sortie = true;
double f_min = 1e-3;

event measure (i++) {

  /**   
  When we produce a drop, we output the simulation. For this reason, we compute
  the vorticity flow to have access to it in the dump file.*/

  scalar omega[];
  vorticity (u, omega);
  /**
  Initialisation of the field before the measure*/
  vector h[]; 
  double xMax = -2;
  double xMin = 5;
  double yMin = L0;
  double velo = 0., veloFilm = 0.;
  // heights (f, h);
  

  /**
  We tag the fluid, to detect the different phases. */

  scalar m[];

  foreach()
    m[] = f[] > (10*f_min);
  int n = tag(m);

  double v[n];
  double dropVelo[n];
  coord b[n];
  Phase drop[n-1];
  int nc[n];

  for (int j = 0; j < n; j++){
    dropVelo[j] = v[j] = b[j].x = b[j].y = b[j].z = 0.;
    nc[j] = 0;
  }

  for(int j =0; j< n; j++)
    initPhase(drop); // initialisation of the drop phase

  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;

      /**
      Calculation of the volume of each phase*/

      v[j] += dv()*f[];

      /**
      Calculation of the number of cells in each drop*/
      
      nc[j]++;

      /**
      Computation of each drop velocity*/

      dropVelo[j] += dv()*f[]*u.x[];

      coord p = {x,y,z};
      foreach_dimension()
        b[j].x += dv()*f[]*p.x;
    }

  #if _MPI // If we work on a supercomputer or simply with several core
    MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, dropVelo, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, nc, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    int world_size;
//     MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//     fprintf(stderr, "N proc %d\n", world_size); // for debug reason
  #endif

  /**
  We detect the main phase. This correspond to the largest one. We use that to
  measure correctly the position of the top of the jet*/
  int tagMainPhase = 1;
  double maxVol = 0;

  if (n>1){
    for (int j = 0; j<n; j++) {
      if (v[j]>maxVol){
        maxVol = v[j];
        tagMainPhase = (j+1);
      }
    }
  }

  scalar g[];
  foreach(){
    if (m[] == tagMainPhase)
      g[] = f[];
    else
      g[] = 0;
  }

  boundary ({g});

  heights (f, h);


  
  /**      
  We will take the maximum interface position for only a slice of the
  simulation. This slice begins at the axi of symmetry and has a height equal 
  to the size of the biggest cells (i.e. $L0/2^7$).

  We will also restrain our test on the main phase*/

  foreach(reduction(max:xMax))
    #if dimension>2 //dimension == 3
      if (fabs(y) <= Delta && fabs(z)<= Delta && h.x[] != nodata && g[]>(0.+f_min) && g[]<(1-f_min)) 
    #else //dimension == 2
      if (y< Delta && h.x[] != nodata && g[]>(0.+f_min) && g[]<(1.-f_min)) 
    #endif
    {
      double xi = x + height(h.x[])*Delta;
      if (xi>xMax)
        xMax = xi;
    }

  foreach(reduction(min:xMin)) 
    if (h.x[] != nodata && g[]>(0.+f_min) && g[]<(1.-f_min)) 
    {
      double xi = x + height(h.x[])*Delta;
      if (xi<xMin)
        xMin = xi;
          // fprintf(stderr, "y = %f, Delta = %f\n", y, Delta);
    }

  double deltaX = xMax-xprev;

  if (xMax == -2 || (fabs(deltaX)>0.5 && xMax<-1.5)){
    dump("error");
    return 0;
  }

  #if FILM // To track a liquid retracting film during a bursting bubble event
    foreach(reduction(min:yMin))
    if(x>0. && h.y[] != nodata) {
      double yi = y+height(h.y[])*Delta;
      if (yi<yMin) {
        yMin = yi;
      }
    }
  #endif

  scalar xpos[];
  scalar ypos[];
  scalar zpos[]; // 3D case

  /**
  Capturing the film position and it's velocity. This is only the closest point
  of the film to the axis of symmetry*/

#if dimension>2 // 3D case
    position(g, xpos, {1,0,0});
    position(g, ypos, {0,1,0});
    position(g, zpos, {0,0,1});

  boundary({xpos, ypos, zpos});
  #if FILM
    scalar r[];
    double epsilon = 0.02;
    scalar out[];
    scalar out2[];

    foreach(){
      if (xpos[] != nodata)
        r[] = sqrt(sq(ypos[])+sq(zpos[]));
      else 
        r[] = nodata;
    }

    boundary({r});

    double rMin = HUGE;
    double xHeight = -HUGE;

    foreach(reduction(min:rMin) reduction(max:xHeight)){
      if(x>0. && r[] != nodata){
        if (rMin > r[])
          rMin = r[];
      if (xpos[] != nodata && xHeight<xpos[])
        xHeight = xpos[];
      }
    }
    double xr = HUGE;
    double yr = HUGE;
    double zr = HUGE;
    double rh = HUGE;
    veloFilm = HUGE;

    /**
    Get point with minimal radius. Compute the velocity for this point*/

    foreach(reduction(min:xr) reduction(min:yr) reduction(min:zr) 
      reduction(min:veloFilm) reduction(min:rh)) {
      if (r[] == rMin){
        xr = xpos[];
        yr = ypos[];
        zr = zpos[];
        veloFilm = sqrt(sq(u.z[])+sq(u.y[]) + sq(u.x[]));
      }
      if (xpos[] == xHeight){
        rh = r[];
      }
    }
    fprintf(stderr, "rh %g xH %g xr %g\n", rh, xHeight, xr);

    /**
    tag cells closest to the minimum radius*/
    foreach() {
      if ((fabs(r[]-rh)/rh<=0.01 && x > xr) ){
        out2[] = 1;
        out[] = 0;
      }

      /**
      tag cells closest to the maximum height of the film*/
      else if ((fabs(x-xr)<= Delta/2. && xpos[] != nodata) && r[]< rh){
        out[] = 1;
        out2[] = 0;
      }
      else{
        out[] = 0;
        out2[] = 0;
      }
    }

    boundary({out});
    boundary({out2});


    /**
    Output of the data in two files: coordClose and coordHeight. coordClose
    correspond to the closest coordinate from the axis of symmetry. coordHeight
    corrspond to the hieghest one.

    If we are in MPI, output in the separated file (one per core)*/

    #if _MPI
    char fileFilmClose[80];
    sprintf(fileFilmClose, "filmCoordClose-%d", pid());
    FILE * fpFilmClose = fopen(fileFilmClose, "a");

    char fileFilmHeight[80];
    sprintf(fileFilmHeight, "filmCoordHeight-%d", pid());
    FILE * fpFilmHeight = fopen(fileFilmHeight, "a");
    #else //No MPI
    static FILE * fpFilmClose = fopen("filmCoordClose", "w");
    static FILE * fpFilmHeight = fopen("filmCoordHeight", "w");
    #endif //MPI

    foreach(){
      if (out[] > 0.5){
        fprintf(fpFilmClose, "%g %d %g %g %g %g %g %g %g\n" ,t, i,
         xpos[], ypos[], zpos[], u.x[], u.y[], u.z[], out[]);
      }
      if (out2[] > 0.5){
        fprintf(fpFilmHeight, "%g %d %g %g %g %g %g %g %g\n" ,t, i,
         xpos[], ypos[], zpos[], u.x[], u.y[], u.z[], out2[]);
      }
    }
    #if _MPI
    fclose(fpFilmHeight);
    fclose(fpFilmClose);
    #endif //MPI
    
    #endif //FILM

  #endif //3D

  if (xMax<-2)
    dump("toto");
  
  foreach(){
    if (x+ height(h.x[])*Delta == xMax )
      velo = u.x[];
  #if dimension < 3
    if (y+height(h.y[])*Delta == yMin)
      veloFilm = u.y[];
  #endif // 3D

  }

  
  /**   
  We also want to output the velocity in a dimensionless way. The
  operation to have a dimensionless velocity is:

  $$ V_{adim} = \frac{\rho\mu V_{tip}}{\sigma} $$
  */


  double deltaT = t - tprev;

  #if _MPI

    velo = (xMax-xprev)/deltaT;

      #if FILM
        #if dimension >2
        // veloFilm = (rMin - rprev)/deltaT;
        #endif //3D
      #endif //FILM

  #endif //MPI

  #if FILM
  double veloFilmAdim = veloFilm*sqrt(thick/2.);
  #endif
  double veloAdim = velo*sqrt(1./La);


  if (xprev<0. && xMax>0.) {
    tcompt++;
    fprintf (stderr, "crossing %g %f\n", t, veloAdim);

    if (tcompt==1){
      tCross = t;

      /**
      We output the simulation when the jet cross the axis $x=0$ .*/

      char name[80];
      sprintf(name, "jetX-0-t-%07f",t);
      dump (name);
    }
  }

  #if dimension>2
    #if FILM
      #if _MPI
        if (pid() == 0)
          printf("%d %f %g %f %f %f %f %f %f %f %f %f %d\n", 
          i, t, deltaT, xMax, velo, veloAdim, 
          rMin, xr, yr, zr, veloFilm, veloFilmAdim, n);
          // printf("%d %f %g %f %f %f %f %f %f %d\n", 
          // i, t, deltaT, xMax, velo, veloAdim, 
          // rMin, veloFilm, veloFilmAdim, n);
      #else //Single Core
        printf("%d %f %g %f %f %f %f %f %f %f %f %f %d\n", 
        i, t, deltaT, xMax, velo, veloAdim, 
        rMin, xr, yr, zr, veloFilm, veloFilmAdim, n);
      #endif //MPI
    #else //No film, but in 3D
      #if _MPI
        if (pid() == 0)
          printf("%d %f %g %f %f %f %f %d\n", 
            i, t, deltaT, xMax, velo, veloAdim, xMin, n);
      #else //No mpi
        printf("%d %f %g %f %f %f %f %d\n", 
          i, t, deltaT, xMax, velo, veloAdim, xMin, n);
      #endif //MPI
    #endif //Film
  #else //no 3D

    #if _MPI
    if (pid() == 0) {
      #if FILM
      printf("%d %f %g %f %f %f %f %f %f %d\n", 
        i, t, deltaT, xMax, velo, veloAdim, yMin, veloFilm, veloFilmAdim, n);
      #else //No film, but MPI, no 3D
      if (i==0)
        printf("i, t, deltaT, xMax, velo, veloAdim (Ca), xMin, n");
      printf("%d %f %g %f %f %f %d\n", i, t, deltaT, xMax, velo, veloAdim, n);
      #endif
    }
    #else //no mpi
      #if FILM 
        printf("%d %f %g %f %f %f %f %f %f %d\n", 
          i, t, deltaT, xMax, velo, veloAdim, yMin, veloFilm, veloFilmAdim, n);
      #else //No film, no mpi, no 3D
	if (i==0)
	  printf("i, t, deltaT, xMax, velo, veloAdim (Ca), xMin, n");

        printf("%d %f %g %f %f %f %f %d\n", i, t, deltaT, xMax, velo, veloAdim, xMin, n);
      #endif //film
    #endif //MPI
  #endif // 3D

  fflush (stdout);

  tprev = t;
  xprev = xMax;
  yprev = yMin;

  #if dimension>2
    #if FILM
      rprev = rMin;
    #endif //FILM
  #endif //3D

  /**   
  We will output those data in 2 external file, call "dropCorr" and
  "dropErr". DropCorr correspond to correct droplets. dropErr is where we put 
  the droplets that might be incorrect.*/

  static FILE * fp = fopen("dropCorr", "w");
  static FILE * fp2 = fopen("dropErr", "w");

  /**
  We first compute all the variables.

  We calculate:

   - the position of the drop
   - the velocity of the drop
   - the dimensionless velocity: $v\mu/\gamma$ 
   - the volume of the drop
   - the radius of the drop
   - the Reynolds number of the drop
   - the number of cells (at maximum level) in the radius
   - the number of cells in the drop*/

  for (int j =1; j < n; j++) {
    int k = j-1;
    drop[k].tag = j+1; 
    drop[k].x = b[j].x/v[j];
    drop[k].velo = dropVelo[j]/v[j];
    drop[k].veloAdim = dropVelo[j]/v[j]*sqrt(1/La);
    drop[k].volume = v[j];
    #if dimension > 2
      drop[k].radius = pow((3./(4*pi)*drop[k].volume),(1./3.));
    #endif
    #if AXI
      drop[k].radius = pow((3./2.*drop[k].volume),(1./3.));
    #endif
    drop[k].reynolds = drop[k].radius*fabs(drop[k].velo)*sqrt(La);
    drop[k].cellsRadius = drop[k].radius*(1 << LEVEL)/L0;
    drop[k].cells = nc[j];

    
  }

  /**
  We sort the variables. We want to have access to the bubble with respect
  to their order of creation.*/

  int size = (n-1);

  /**
  The first created drop is the highest one. So, we will sort all the
  variable by sorting the droplet from the highest one to the lowest one.*/

  qsort(drop, sizeof drop / sizeof *drop, sizeof *drop, comparePosition);

  /**
  There is actually (16/02/18) a bug in the tag function. It's tagging
  some very small area. To avoid that, we will print our output a second time,
  but this time, we will skip the line with a drop volume below $3.10^{-6}.
  Since we don't want to have a "bad" suppression, we are outputting this in a
  second file.

  We also add a count variable call phase. It will count the number of phases 
  to know if we have to generate a new dump file or not.*/

  int numPhase = 0;

  /**
  We sort the droplet. The one with less than a certain amount of cells will be
  outputed in the error file*/

  int minCells = 60;

  for (int j = 0; j< size; j++) {
    fprintf((drop[j].cells>minCells) ? fp : fp2 , "%d %g %d %g %g %g %g %g %g %g %d\n",
	    i, t,(j+1), drop[j].x, drop[j].velo, drop[j].veloAdim,
	    drop[j].volume, drop[j].radius, drop[j].reynolds,
	    drop[j].cellsRadius, drop[j].cells);

    if (j == (size-1)){
      fprintf(fp, "\n"); // dropCorr
      fprintf(fp2, "\n"); // dropErr
    }

    if (nc[j]> minCells)
      numPhase++;

    fflush(fp);
    fflush(fp2);
  }

  /**   
  If a drop is very small, we add the possibility to remove it. This part
  is still under developpement.œ

  We remove the very small droplet based on the 2 following conditions:
   - The droplet's cells number is below 60
   - The droplet is at least at 0.2 from the top of the jet*/

  // if (i>10)
  //   foreach()
  //     for (int j = 0; j < size; j++)
  //       if (drop[j].tag == m[])
  //         if (drop[j].cells<minCells && drop[j].x > (xMax+0.2)){
	 //          f[] = 0;
	 //          fprintf(stderr,"at i= %i t=%g a droplet has been removed\n", i, t);
	 //        }

  // boundary ({f});


  /**
  We can dump some simulation, to observe the evolution of the bursting
  bubble. To avoid too many output we restrain ourselves to 3 dump every 500 
  steps.

  We dump a simulation if there is a change in the number of tagged phases.*/

  #if DumpSimulation

  if (i%500==0)
    sortie = true;
  /**
  We dump a simulation if there is a new drop. But we avoid to have a new file
  at each time steps*/
  if ((numPhase > phasePrev || compt != 0)  && i>30 && sortie) {
    compt++;
    char rupt[80];
    sprintf( rupt, "rupture-%07d", i);
    dump (file = rupt);

    if (compt == 2){
      compt = 0;
      sortie = false;
    }

    if (compt == 1) {
      fprintf(stderr, "drop production, i = %i, t = %f\n", i, t);
    }

  }

  #endif

  /**
  At the beginning of the simulation, there is sometimes some small
  droplet that should not be here. We remove them. */

  // if (i >= 10 && i <= 30)
  //   remove_droplets_vol(f, 1, true);

  /**   
  When the simulation grid became too big, we stop it. This should avoid
  an overload of the RAM. We also add a dump file of this state, in order to 
  chek the simulation and see what happens. */


  // if (grid->tn > 80000){
  //   char dumpFile[80];
  //   sprintf (dumpFile, "dump");
  //   dump (file = dumpFile);
  //   fprintf (stderr, "grid overpass 80,000 cells. Abort.\n");
  //   return 1;
  // }

  phasePrev = numPhase;

}

/** 
## Time detection

We want to output the facets for a precise time. However, using a time event may
change the results (since we force the timestep). We will just test the sign of
a variable. When the sign change, we output the state of the simulation.

We will perform an output at t = 0.3, t = 0.35, t = 0.4 and t = 0.45.

This event can be disabled by defining DetectTime to 0

We will use the variable tPrev.

This second event will occur after the event measure, so the value of tPrev will
already be the "good" one*/

#define DetectTime 0


double tPrev = -1;

#if DetectTime

void timeOutput (double t, double tPrevious, double target) {
  if (timeDetect>0){
    scalar omega[];
    vorticity(u, omega);
    char name[80];
    sprintf(name, "interface-t-%g", target);
    FILE * fp  = fopen(name,"w");
    output_facets(f, fp);
    char dumpName[80];
    sprintf(dumpName, "simulation-t-%g", target);
    dump(file = dumpName);
  }
}

event detectTime (i++) {
  double tTarget1 = 0.3;
  double tTarget2 = 0.35;
  double tTarget3 = 0.4;
  double tTarget4 = 0.45;

  #if AXI
    #if ROMEO
      for (int i = 0; i<250; i++)
        timeOutput(t, tPrev, 0.35+i/1000.);
    #else
      for (int i = 0; i <80; i++)
        timeOutput(t, tPrev, 0.46+i/200.);
    #endif
  #endif
  
  timeOutput(t, tPrev, tTarget1);
  timeOutput(t, tPrev, tTarget2);
  timeOutput(t, tPrev, tTarget3);
  timeOutput(t, tPrev, tTarget4);

  tPrev = t;
}

#endif

#if ROMEO

/**
## Physical measurement: bubble

Dedicated event for the measure of the bubble under the cavity. The behavior is
similar to the event for the drops. We just replace the drops by the bubble*/

event detectBubble(i++){
  scalar m[];

  /**
  Tag of the gaz phase*/

  foreach()
    m[] = f[] < (1-10*f_min);
  int n = tag(m);
  fprintf(stderr, "n = %d\n", n);

  if (n>1) {  
    /**
    Initialisation of the variable*/
    double v[n];
    double bubbleVelo[n];
    coord b[n];
    Phase bubble[n-1];
    initPhase(bubble);
    int nc[n];

    for (int j = 0; j < n; j++){
      bubbleVelo[j] = v[j] = b[j].x = b[j].y = b[j].z = 0.;
      nc[j] = 0;
    }

    foreach_leaf()
      if (m[] > 0) {
        int j = m[] - 1;

        /**
        Calculation of the volume of each phase*/

        v[j] += dv()*(1-f[]);

        /**
        Calculation of the number of cells in each bubble*/
        
        nc[j]++;

        /**
        Computation of each bubble velocity*/

        bubbleVelo[j] += dv()*(1-f[])*u.x[];

        coord p = {x,y,z};
        foreach_dimension()
          b[j].x += dv()*(1-f[])*p.x;
      }

    #if _MPI
      MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce (MPI_IN_PLACE, bubbleVelo, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce (MPI_IN_PLACE, nc, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      int world_size;
      MPI_Comm_size(MPI_COMM_WORLD, &world_size);
      fprintf(stderr, "N proc %d\n", world_size);
    #endif

    /**
    We detect the main gaz phase (the biggest one). We only */

    int tagMainPhase = 0;
    double maxVol = v[0];

    for (int j = 0; j<n; j++) {
      if (v[j]>maxVol){
        maxVol = v[j];
        tagMainPhase = (j+1);
      }
    }

    fprintf(stderr, "tag air %d\n", tagMainPhase);
    
    /**     We will output those data in 2 external file, call "bubbleCorr"
and"bubbleErr". buubbleCorr correspond to correct bubble. bubbleErr is where we
put the bubble that might be incorrect.*/

    static FILE * fp = fopen("bubbleCorr", "w");
    static FILE * fp2 = fopen("bubbleErr", "w");

    /**
    We first compute all the variables.

    We calculate:

     - the position of the bubble
     - the velocity of the bubble
     - the dimensionless velocity: $v\mu/\gamma$ 
     - the volume of the bubble
     - the radius of the bubble
     - the Reynolds number of the bubble
     - the number of cells (at maximum level) in the radius
     - the number of cells in the bubble*/

    for (int j =0; j < n; j++) {
      bubble[j].tag = j+1; 
      bubble[j].x = b[j].x/v[j];
      bubble[j].velo = bubbleVelo[j]/v[j];
      bubble[j].veloAdim = bubbleVelo[j]/v[j]*sqrt(1/La);
      bubble[j].volume = v[j];
      #if dimension > 2
        bubble[j].radius = pow((3./(4*pi)*bubble[j].volume),(1./3.));
      #endif
      #if AXI
        bubble[j].radius = pow((3./2.*bubble[j].volume),(1./3.));
      #endif
      bubble[j].reynolds = 1./998.*bubble[j].radius*fabs(bubble[j].velo)*sqrt(La);
      bubble[j].cellsRadius = bubble[j].radius*(1 << LEVEL)/10.;
      bubble[j].cells = nc[j];
    }
        
    int numPhase = 0;

    int minCells = 10;

    for (int j = 0; j< n; j++) {
      if (bubble[j].tag!=tagMainPhase){
        fprintf((bubble[j].cells>minCells) ? fp : fp2 , "%d %g %d %g %g %g %g %g %g %g %d\n",
          i, t,(j+1), bubble[j].x, bubble[j].velo, bubble[j].veloAdim,
          bubble[j].volume, bubble[j].radius, bubble[j].reynolds,
          bubble[j].cellsRadius, bubble[j].cells);
        }
    
      if (nc[j]> minCells)
        numPhase++;

    }

    fprintf(fp, "\n");
    fprintf(fp2, "\n");

    fflush(fp);
    fflush(fp2);

    /**
    Dump of the simulation if there is production of a small bubble. We output the steps to have a chance to see it at formation*/

    #if DumpSimulation

    if (i%500==0)
      sortie = true;

    if ((numPhase > phasePrev)  && i>30 && sortie) {
      compt++;
      char rupt[80];
      sprintf( rupt, "bubble-%07d", i);
      dump (file = rupt);

    }

    #endif
    phasePrev = numPhase;
  }
}
#endif
/**
## Tracking the divergence field*/
#if DIV>2


event divTracking(i++; i <= Iend) {
  scalar divU[];
  foreach(){
    divU[] = 0;
    if (f[]>1e-6){
      foreach_dimension()
        divU[] += (uf.x[1] - uf.x[-1])/(2.*Delta);
    }
  }
  fprintf(stderr, "%d div max %f div mean %f\n",
   i, normf(divU).max, normf(divU).avg);
  {
    static FILE * fpNoise = popen ("ppm2mp4 divField.mp4", "w");
    clear();
    view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
    draw_vof("f");
    squares("divU", min = -5, max = 5);
    mirror({0,1,0}, 0) {
      draw_vof("f");
      squares("divU", min = -50, max = 50);
    }
    save(fp = fpNoise);
  }
}
#endif
/**
## Interface tracking in the noise case*/
#if IntTracking

event interfaceTracking(i+=5; i<=150) {
  char name[80];
  
  #if (BRUIT || InterfaceNoise) 
    sprintf(name, "interface-bruit-%03d",i);
  #else
    sprintf(name, "interface-%03d",i);
  #endif

  FILE * fp = fopen(name, "w");
  output_facets(f, fp);
}

#endif
