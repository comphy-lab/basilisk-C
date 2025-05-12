/**
[Return to my homepage](http://basilisk.fr/sandbox/nlemoine/README)

# Hydrodynamic modeling of the Meuvaines / Ver-sur-Mer coastal wetland (Calvados, France)

![Picture : [Conservatoire du Littoral](https://www.conservatoire-du-littoral.fr/siteLittoral/252/28-marais-de-ver-14_calvados.htm)](https://www.conservatoire-du-littoral.fr/uploads/Image/b3/6811_261__77H6444.jpg){width="500px"}
*/

#include "grid/quadtree.h"
#include "run.h"
#include "diffusion.h"
#include "input.h"
#include "output.h"
#include "terrain.h"
#include "nlemoine/vector-geom.h"
#include "view.h"
#include "nlemoine/view-utils.h"
#include "nlemoine/groundwater/Meuvaines/colormaps.h"

/**
## Definitions and global variables
*/
#define LEVEL 8
#define MINLEVEL 5
#define HE 1.

/**
We first define some space / time constants 
*/

// Chart datum and mean sea level (RAM_PACK data from SHOM)
#define ZH_ref_Arromanches -4.019
#define NM_Arromanches 4.43
#define ZH_ref_Courseulles -3.99
#define NM_Courseulles 4.59

/** 1-hour runtime limit on the server does not allow to simulate more than about 18 month with `LEVEL = 8`. Just change the `tfin` value (in days) to simulate a longer period, e.g. `#define tfin 1462.` to simulate from 1963-01-01 until 1966-12-31 and validate against observed streamflow in 1965-1966 (a two-year warm-up period covering 1963-1964 is necessary to remove the effet of the (unknown) initial head distribution.
*/
#define tfin 518.
#define DAYS_PER_YEAR 365.25
#define DAYS_PER_MONTH 30.4375

// Some reference points (RGF Lambert 93, EPSG:2154)
double XlimW,XlimE,YlimS,YlimN;
coord SegNeumannW = {430330.,6912160.};
coord SegNeumannE = {435380.,6911570.};
coord n_int;
double ZMER;
double datenum0 = 1614.; // initial date expressed as elapsed days since 1958-08-01, hence 1963-01-01

/** A small hack to convert elapsed simulation time into year & month
   (avoids bringing out the big guns of "time.h")
   Assumes that `t+datenum0` is the number of (fractional) days elapsed since
   the begining of SAFRAN forcings on 1958-08-01:
*/

int hack_date(double tt, int * _year, int * _month)
{
  int elapsed_years_since_1958 = (int) floor((tt+datenum0+212)/DAYS_PER_YEAR);
  int elapsed_months_since_1958 = (int) floor((tt+datenum0+212)/DAYS_PER_MONTH);
  *_month = (elapsed_months_since_1958 % 12) + 1; 
  *_year = 1958 + elapsed_years_since_1958;
  return(0);
}

typedef struct {
  coord Site;
  char BSS_ID[12];
  int layer;
} Piezo ;

struct TimeSeries Forcings;
double CumP,CumPE,P_dt,PE_dt;
double SMA0,SMA1,SMA2,Excess0,Excess1,Excess2;
double maxStorage0; // marls soil storage
double * Forc_t;
double oneday = 1. [0,1];
double sum_dt_month;

double PARAM[11] = {200., 200., 2.e-4, 1.e-5, 0.08, 3.e-4, 150., 1.5e-5, 1.e-4, 0.025, 1.e-5};
int nparam;

scalar isNeumann[], isDirichlet[];
scalar zb[];
scalar msk1[],msk2[];

struct Polygon PolyDirichlet;
struct Polygon BV_Gronde;
struct Polygon BV_Provence;
double QM_Gronde, QM_Provence;
Piezo * PIEZO;
int npiezo;
double * HPIEZO;

scalar maliere_mur[], marnesPB_mur[], bathonien_mur[];

/** All scalars are indexed (1) for the semi-confined aquifer in Bajocian limestone, and (2) for the unconfined aquifer in mid-Bathonian limestone. 
*/
scalar zbottom1[],zbottom2[]; // elevation of the base of each aquifer layer
scalar ztop1[],ztop2[];
scalar h1[],h2[]; // groundwater head
face vector T1[], T2[]; // transmissivity
scalar S1[],S2[]; // storage
scalar r1[],r2[]; // recharge
scalar beta1[],beta2[];
scalar seep1[],seep2[];
scalar isConfined1[],isConfined2[]; 
scalar hasThickness1[],hasThickness2[];
scalar recharge_layer[];

double dt;
mgstats mgd;

// Definition of hydrodynamic attributes (uniform for the moment, but could be scalar [])
// They will be attached to the scalars zbottomN[]

attribute {  
  double maxStorage;
  double Ksat; // Saturated hydraulic conductivity
  double Tf; // fixed transmissivity component
  double Ss; // Spectific storage (confined case)
  double omega_d; // porosité efficace (cas libre)                           
} 
 
FILE * resfile, * flowfile;

/**
## Main
*/
#include "nlemoine/groundwater/Meuvaines/aquifer_ML.h"

int main ()
{
  if(pid()==0){
    
    system("wget https://dropsu.sorbonne-universite.fr/s/Zg24YkSA3w2e7wk/download && mv download data.zip");
    system("unzip data.zip");
    
    resfile = fopen("results.csv","wt");
    flowfile = fopen("streamflow.dat","wt");
  }

  run(); 

  fclose(resfile);
  fclose(flowfile);
  return(0);
}

/**
## Initialization
*/
event init (i = 0)
{
  // Mean sea level

  ZMER = 0.5*((NM_Arromanches+ZH_ref_Arromanches)+
              (NM_Courseulles+ZH_ref_Courseulles));
    
  // Max resolution: (2^12) * (2^12) cells of 5 m 

  L0 = 5. * (double)(1 << 12);

  // Initialize grid

  XlimW = 429162.5; // East of Port-en-Bessin
  XlimE = 449102.5; // Courseulles
  YlimS = 6911432.5;
  YlimN = YlimS+L0;

  size (L0);
  X0 = (XlimW+XlimE)/2. - L0/2;
  Y0 = YlimS;
  origin (X0, Y0);
  N = 1 << LEVEL;
  init_grid (1 << LEVEL);

  // Inward normal to the domain at the south border,
  // defined by a segment approximately following the N13 road south of Bayeux
  // between Aure river (SegNeumannW) and Seulles river (SegNeumannE)

  n_int = (coord){-(SegNeumannE.y-SegNeumannW.y),+(SegNeumannE.x-SegNeumannW.x)};
  double norm = sqrt(sq(n_int.x) + sq(n_int.y));
  n_int.x /= norm;
  n_int.y /= norm;

  // Read polygon specifying Dirichlet conditions (Aure and Seulles rivers)

  (void) ReadPolygon ( "PolyDirichlet_v2.dat", & PolyDirichlet);
  printf("%d vertices read.\n",PolyDirichlet.nv);

  // Read contours of gauged catchments

  (void) ReadPolygon ( "Gronde_SGR_simplif.dat", & BV_Gronde);
  printf("%d vertices read.\n",BV_Gronde.nv);

  (void) ReadPolygon ( "Provence_SGR_simplif.dat", & BV_Provence);
  printf("%d vertices read.\n",BV_Provence.nv);

  // Read locations of observation boreholes
    
  (void) read_piezo_sites ( & PIEZO, "BSS_Basilisk.dat", & npiezo);

  HPIEZO = (double *) malloc(npiezo*sizeof(double));

  if(pid()==0)
  {
    fprintf(resfile,"# PARAM :");  
    for(int kp=0;kp<nparam;kp++)
      fprintf(resfile," %g",PARAM[kp]);
    fprintf(resfile,"\n");
          
    fprintf(resfile,"t,P,PE,Q_Gronde,Q_Provence");
    for(int kp=0;kp<npiezo;kp++)
      fprintf(resfile,",%s",PIEZO[kp].BSS_ID);
    fprintf(resfile,"\n");
  }
    
  // Load bathy-topographic data

  terrain (zb, "bessin_topo", NULL);
  boundary({zb});

  // Load interface elevations

  terrain (maliere_mur, "maliere_mur", NULL);
  terrain (marnesPB_mur, "marnesPB_mur", NULL);
  terrain (bathonien_mur, "bathonien_mur", NULL);

  (void) update_geom_stack();
    
  /** We set the [hydrodynamic properties](https://books.gw-project.org/hydrogeologic-properties-of-earth-materials-and-principles-of-groundwater-flow/chapter/properties-of-aquifers-and-confining-units/) of each aquifer layer according to the content of the `PARAM` vector. Specific storage (for the confined case) is a volume of water released by a unit volume of porous medium in response of a unit head drop, hence in $\textrm{L}^3\cdot\textrm{L}^{-3}\cdot\textrm{L}^{-1}$ i.e. $\textrm{L}^{-1}$.*/
    
    // Soil storage over marls (0)

    maxStorage0 = PARAM[0]; // in mm
    
    // Bajocian semi-confined aquifer (1)
    
    zbottom1.maxStorage = PARAM[1];
    zbottom1.Ksat = PARAM[2]*86400.; // convert m/s => m/day
    zbottom1.Tf = PARAM[3]*86400.;  // convert m2/s => m2/day   
    zbottom1.omega_d = PARAM[4];
    zbottom1.Ss = PARAM[5]; // specific storage in m^-1
    
    // mid-Bathonian unconfined aquifer (2)
    
    zbottom2.maxStorage = PARAM[6];
    zbottom2.Ksat = PARAM[7]*86400.; // convert m/s => m/day
    zbottom2.Tf = PARAM[8]*86400.; // convert m/s => m/day
    zbottom2.omega_d = PARAM[9];
    zbottom2.Ss = PARAM[10]; // specific storage in m^-1
    
  // Initialized isDirichlet[] and isNeumann[] scalars

  (void) update_BC_cells();

  // Initialize groundwater heads

  scalar zbmin = zb.dmin;
  scalar msk[];
  scalar msk_base_maliere[], msk_base_marnesPB[], msk_base_bathonien[];
  scalar msk_outcrop_maliere[], msk_outcrop_marnesPB[], msk_outcrop_bathonien[];

  foreach()
  {
    msk[] = (isDirichlet[]>0.) | (isNeumann[]>0.) ? -1. : 1. ;
    double drainage_level = zbmin[]<=zb[] ? fmax(zbmin[],ZMER) : fmax(zb[],ZMER);
      
    h1[] = fmax(zbottom1[],drainage_level);
    seep1[] = 0.;
    
    h2[] = fmax(zbottom2[],drainage_level);
    seep2[] = 0.;

    msk_base_maliere[]    = (maliere_mur[]  <=zb[])      && (msk[]>0.) ? 1. : -1.;
    msk_base_marnesPB[]   = (marnesPB_mur[] <=(zb[]-3.)) && (msk[]>0.) ? 1. : -1.;
    msk_base_bathonien[]  = (bathonien_mur[]<=zb[])      && (msk[]>0.) ? 1. : -1.;

    msk_outcrop_maliere[]    = (maliere_mur[]  <=zb[])      && (marnesPB_mur[] >=(zb[]-3.)) && (msk[]>0.) ? 1. : -1.;
    msk_outcrop_marnesPB[]   = (marnesPB_mur[] <=(zb[]-3.)) && (bathonien_mur[]>=zb[])      && (msk[]>0.) ? 1. : -1.;
    msk_outcrop_bathonien[]  = msk_base_bathonien[];
  }

  // Illustrate interfaces
  float qview_x[4],qview_z[4],qview[4];
  (void) gl_axis_to_quat ((float[]){1,0,0}, 0.42*PI, qview_x); // angle w.r.t vertical
  (void) gl_axis_to_quat ((float[]){0,0,1}, 0.72*PI, qview_z); // rotation in horizontal plane
  gl_add_quats(qview_x, qview_z, qview);
  
  char str[120];

  view (fov = 13.5, quat = {qview[0],qview[1],qview[2],qview[3]},
        sx = 1., sy = 1., sz = 40.,
        width = 1200, height = 768);

  translate(x = -X0-L0/2.+0.12*L0,y = -Y0-L0/2.+0.12*L0,z=+40.){
    NCCLASS = 25;
    masked_squares("maliere_mur", linear = true, z = "maliere_mur", mask = msk_base_maliere, min = -170., max = 80., map = discrete_blues);
    masked_squares("marnesPB_mur", linear = true, z = "marnesPB_mur", mask = msk_base_marnesPB, min = -140., max = 110., map = discrete_grays);
    masked_squares("bathonien_mur", linear = true, z = "bathonien_mur", mask = msk_base_bathonien, min = -120., max = 130., map = discrete_reds);
    surf_cells(zb,mask=msk_outcrop_maliere,lc={0.,0.,0.7});
    surf_cells(zb,mask=msk_outcrop_marnesPB,lc={0.4,0.4,0.4});
    surf_cells(zb,mask=msk_outcrop_bathonien,lc={0.7,0.,0.});

    sprintf (str, "                             Blue: Base (solid) and outcrops (wireframe) of Bajocian limestone");
    draw_string (str, 1, size = 100, lc = {0,0,1}, lw = 1);  // from top-left

    sprintf (str, "\n                             Gray: Base (solid) and outcrops (wireframe) of lower-Bathonian marls");
    draw_string (str, 1, size = 100, lc = {0,0,0}, lw = 1);  // from top-left

    sprintf (str, "\n\n                             Red:  Base (solid) and outcrops (wireframe) of mid-Bathonian limestone");
    draw_string (str, 1, size = 100, lc = {1,0,0}, lw = 1);  // from top-left

    sprintf (str, "Z-stretch x 40 ");
    draw_string (str, 3, size = 64, lc = {0,0,0}, lw = 1); // bottom-right
  }
  save("interfaces.png");
  clear();

/**
## 3D view of the geological structure
The coastal wetland is at the edge of the outcrops of both Bajocian limestone and mid-Bathonian limestone. Geological layers have a tilt oriented towards the North-East direction.

![3D geological structure](meuvaines/interfaces.png){width="700px"}
*/
  
  // Initialize T & S maps

  (void) update_TS(zbottom1,ztop1,h1,T1,S1,isConfined1,hasThickness1,seep1);
  (void) update_TS(zbottom2,ztop2,h2,T2,S2,isConfined2,hasThickness2,seep2);
  boundary(all);

  // Load cumulative P & PE time series
  const char * separators = " ,;|";
  (void) load_timedata(& Forcings,"safran.csv",1,separators,datenum0); 
    
  // Initialize cumulative P & PE
  Forc_t = (double *) malloc(2*sizeof(double));
  (void) interpolate_timedata(& Forcings,0.,&Forc_t);
  CumP = Forc_t[0];
  CumPE = Forc_t[1]; 
   
  SMA0 = 0.;
  SMA1 = 0.;
  SMA2 = 0.;
  sum_dt_month = 0;
  QM_Gronde = 0.;
  QM_Provence = 0.;
}

/**
## Final event (assemble movies)
*/
event stop (t = tfin)
{  
  system ("for f in seepage-*.png; do convert $f ppm:- && rm -f $f; done | "
	  "ppm2mp4 seepage.mp4");
  system ("for f in gw_levels-*.png; do convert $f ppm:- && rm -f $f; done | "
	  "ppm2mp4 gw_levels.mp4");
  fprintf(stderr,"Done.\n");
  return 1;
}

/**
## Time integration
*/
event integration (i++)
{
  scalar S1back[],S2back[];
  scalar r1back[],r2back[];
  
  dt = dtnext(2.);  // fixed time step = 2 days
    
  (void) update_geom_stack();  
  (void) update_TS(zbottom1,ztop1,h1,T1,S1,isConfined1,hasThickness1,seep1);
  (void) update_TS(zbottom2,ztop2,h2,T2,S2,isConfined2,hasThickness2,seep2);
  (void) apply_head_conditions();

  (void) interpolate_timedata(& Forcings,(t+dt)/oneday,&Forc_t);
  P_dt = Forc_t[0]-CumP;
  PE_dt = Forc_t[1]-CumPE;  
  CumP = Forc_t[0];
  CumPE = Forc_t[1];
  
  /** The soil moisture accounting (SMA) procedure is taken from the [GR4J model](https://webgr.inrae.fr/webgr-eng/tools/hydrological-models/daily-hydrological-model-gr4j). It is computed in the three zones defined, namely over impervious marls (0), over outcropping Bajocien limestone (1), and over outcropping mid-Bathonian limestone (2):
  */
  (void) SMA_GR(& SMA0, maxStorage0,P_dt,PE_dt,(dt/oneday),false,& Excess0); 
  double excess_rate0 = (1.e-3)*Excess0/dt;

  (void) SMA_GR(& SMA1, zbottom1.maxStorage,P_dt,PE_dt,(dt/oneday),false,& Excess1); 
  double excess_rate1 = (1.e-3)*Excess1/dt;

  (void) SMA_GR(& SMA2, zbottom2.maxStorage,P_dt,PE_dt,(dt/oneday),false,& Excess2); 
  double excess_rate2 = (1.e-3)*Excess2/dt;
      
  foreach()
  {
    r1[] = recharge_layer[]==1. ? excess_rate1 : 0.; 
    r2[] = recharge_layer[]==2. ? excess_rate2 : 0.; 
    beta1[] = 0.;
    beta2[] = 0.;
    S1back[] = S1[];
    r1back[] = r1[];
    S2back[] = S2[];
    r2back[] = r2[];
  }

/**
Then we solve the diffusion equation over the timestep. Note that the $r$, $\beta$ and $\theta$ fields will be modified by the solver.
*/
  mgd = diffusion(h1,dt,T1,r1,beta1,S1);
  mgd = diffusion(h2,dt,T2,r2,beta2,S2);
  
  foreach()
  {
      S1[] = S1back[];
      r1[] = r1back[];
      S2[] = S2back[];
      r2[] = r2back[];
  }
    
  (void) update_seepage();
  
  // Sum seepage fluxes in catchments to estimate discharge
    
  double Q_Gronde = 0.;
  double Q_Provence = 0.;
    
  foreach(reduction(+:Q_Gronde) reduction(+:Q_Provence))
  {
     coord P = (coord){x,y};  
     int inGronde = isInPolygon(P, & BV_Gronde );
     int inProvence = isInPolygon(P, & BV_Provence );
 
     if(inGronde>0.)
     {
        if(seep1[]>0.) Q_Gronde += seep1[]*sq(Delta)*oneday/86400.;
        if(seep2[]>0.) Q_Gronde += seep2[]*sq(Delta)*oneday/86400.;
        if(recharge_layer[]==0.) Q_Gronde += excess_rate0*sq(Delta)/86400.;   
     }       

     if(inProvence>0.)
     {
        if(seep1[]>0.) Q_Provence += seep1[]*sq(Delta)*oneday/86400.;
        if(seep2[]>0.) Q_Provence += seep2[]*sq(Delta)*oneday/86400.;
        if(recharge_layer[]==0.) Q_Provence += excess_rate0*sq(Delta)/86400.;   
     }
  }

  QM_Gronde += Q_Gronde * dt/oneday;
  QM_Provence += Q_Provence * dt/oneday;
  sum_dt_month += dt/oneday;

  // Groundwater head at boreholes

  for(int kp=0;kp<npiezo;kp++)
  {
     Point point = locate ( PIEZO[kp].Site.x, PIEZO[kp].Site.y );
     if (point.level < 0 ) // borehole is not in the process subdomain
       HPIEZO[kp] = -9999.;
     else
     {
        if( PIEZO[kp].layer == 1 ) HPIEZO[kp] =  interpolate_linear(point,h1,PIEZO[kp].Site.x, PIEZO[kp].Site.y);   
        if( PIEZO[kp].layer == 2 ) HPIEZO[kp] =  interpolate_linear(point,h2,PIEZO[kp].Site.x, PIEZO[kp].Site.y);   
      }
  }
    
@if _MPI
    MPI_Allreduce (MPI_IN_PLACE, HPIEZO, npiezo, MPI_DOUBLE, MPI_MAX,MPI_COMM_WORLD);
@endif    

  if (pid() == 0) {
    fprintf(resfile,"%f,%.2f,%.2f,%.1f,%.1f",t+oneday*datenum0,P_dt*oneday/dt,PE_dt*oneday/dt,1000.*Q_Gronde,1000.*Q_Provence);  
  
    for(int kp=0;kp<npiezo;kp++)
      fprintf(resfile,",%.2f",HPIEZO[kp]);
    fprintf(resfile,"\n");
    fflush(resfile);
  }

  boundary(all);
}

/**
## Grid adaptivity settings
*/
int adapt() {

#if TREE
  scalar h1w[],h2w[],maxseep[],indic_MPB[],indic_cote[];
  
  foreach()
  {
    h1w[] = (isDirichlet[]>0.) | (isNeumann[]>0.) ? ZMER : h1[];
    h2w[] = (isDirichlet[]>0.) | (isNeumann[]>0.) ? ZMER : h2[];
    maxseep[] = (isDirichlet[]>0.) | (isNeumann[]>0.) ? 0. : fmax(seep1[],seep2[]);
    indic_MPB[] = (zb[]>=marnesPB_mur[]) && (zb[]<=bathonien_mur[]) ? 1. : 0.;
    indic_cote[] = zb[]>=ZMER ? 1. : 0.;
  }
  
  boundary ({h1w,h2w,maxseep});
  double SE = 5.e-4;
  double IE = 0.1;
    
  astats s = adapt_wavelet ({h1w,h2w,maxseep,indic_MPB,indic_cote}, (double[]){HE,HE,SE,IE,IE},
			    LEVEL, MINLEVEL);

  // Updating isNeumann[] and isDirichlet[] to newly adapted grid
  (void) update_BC_cells();

  boundary(all);
    
//  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}

event do_adapt(t+=8.)
{
  (void) adapt();
}

/**
  Movie output
*/

event movie (t+=2.)
{
  scalar seepmsk1[],seepmsk2[],msk[],l[];
  scalar h1draw[], h2draw[];
  scalar zbpos[],colorcode[],seasurface[],msksea[];
  
  foreach()
  {
//    msk[] = (isDirichlet[]>0.) || (isNeumann[]>0.) || (zb[]<-4.) ? -1. : 1.;
    msk[] = (isDirichlet[]>0.) || (isNeumann[]>0.) ? -1. : 1.;
    msk1[] = (msk[] < 0.) || (hasThickness1[] < 1.) ? -1. : 1. ;
    msk2[] = (msk[] < 0.) || (hasThickness2[] < 1.) ? -1. : 1. ;
    seepmsk1[] = 0.5*(1.+msk1[]) * seep1[] ;
    seepmsk2[] = 0.5*(1.+msk2[]) * seep2[] ;
    l[] = level;
    h1draw[] = h1[];
    h2draw[] = msk2[] > 0 ? fmax(h1[]+0.3,h2[]) : h2[];
    seasurface[] = ZMER;    
    msksea[] = (zb[] >= (ZMER+2.)) ||  (msk[]<0.) ? -1. : 1.; 

    colorcode[] = 0.3;
    if(marnesPB_mur[]<(zb[]-3.))colorcode[] = 0.5;
    if(bathonien_mur[]<zb[])colorcode[] = 0.7;
    if(seep1[]>0. && zb[]>= (ZMER-5.)) colorcode[] = 0.1;
    if(seep2[]>0. && zb[]>= (ZMER-5.)) colorcode[] = 0.9;
  }

  boundary({zbpos});
    
  char str[100];
  int year, month;
  (void) hack_date(t/oneday,& year, & month);

  float qview_x[4],qview_z[4],qview[4];
  (void) gl_axis_to_quat ((float[]){1,0,0}, 0.36*PI, qview_x); // angle w.r.t vertical
  (void) gl_axis_to_quat ((float[]){0,0,1}, 0.7*PI, qview_z); // rotation in horizontal plane
  gl_add_quats(qview_x, qview_z, qview);

  view (fov = 11., quat = {qview[0],qview[1],qview[2],qview[3]},
        sx = 1., sy = 1., sz = 25.,
        width = 1200, height = 768);

  translate(x = -X0-L0/2.+0.05*L0,y = -Y0-L0/2.+0.12*L0){
    masked_squares_checkerboard("colorcode", linear = true, z = "zb", mask = msk, min = 0., max = 1., map = topocolormap);
    surf_cells (zb,mask = msk);
    masked_squares("seasurface", linear = true, z = "seasurface", mask = msksea, min = 0., max = 1., map = uniform_gray);
    surf_cells (seasurface,mask = msksea);

    draw_polygon_on_surface ( & BV_Gronde, zb, lc = {0,1,0}, lw = 3.5, z_offset = 2.5);
    draw_polygon_on_surface ( & BV_Provence, zb, lc = {0,1,0}, lw = 3.5, z_offset = 2.5);

    sprintf (str, "%4.4d-%2.2d",year,month);
    draw_string (str, 0, size = 32, lc = {0,0,0}, lw = 4); // bottom-left

    sprintf (str, "                                             Blue:  Bajocian limestone w/ semi-confined aquifer");
    draw_string (str, 1, size = 98, lc = {0,0,1}, lw = 1);  // from top-left

    sprintf (str, "\n                                             White: Lower-Bathonian marls (impervious layer)");
    draw_string (str, 1, size = 98, lc = {0,0,0}, lw = 1);  // from top-left

    sprintf (str, "\n\n                                             Red:   Mid-Bathonian limestone w/ unconfined aquifer");
    draw_string (str, 1, size = 98, lc = {1,0,0}, lw = 1);  // from top-left

    sprintf (str, "\n\n\n\n\n   Satured color = seepage from underlying aquifer");
    draw_string (str, 2, size = 120, lc = {0,0,0}, lw = 1);  // from top-left

    sprintf (str, " Green outlines:");
    draw_string (str, 1, size = 98, lc = {0,1,0}, lw = 2); // top-left
    sprintf (str, "\n Gronde & Provence catchments");
    draw_string (str, 1, size = 98, lc = {0,1,0}, lw = 2); // top-left
   
    sprintf (str, "Z-stretch x 25 ");
    draw_string (str, 3, size = 64, lc = {0,0,0}, lw = 1); // bottom-right
  }
  char fname[200];
  sprintf(fname,"seepage-%4.4d.png",(int)t);  
  save (fname);
  clear();


  view (fov = 11., quat = {qview[0],qview[1],qview[2],qview[3]},
        sx = 1., sy = 1., sz = 25.,
        width = 1200, height = 768);

  translate(x = -X0-L0/2.+0.05*L0,y = -Y0-L0/2.+0.12*L0){
    NCCLASS = 7;
    masked_squares("h1", linear = true, z = "h1draw", mask = msk1, min = 0., max = 70., map = discrete_blues);
    surf_cells (h1draw,mask = msk1,lc={0.,0.,0.2});
    masked_squares("h2", linear = true, z = "h2draw", mask = msk2, min = 0., max = 70., map = discrete_reds);
    surf_cells (h2draw,mask = msk2,lc={0.2,0.,0.});

    sprintf (str, "%4.4d-%2.2d",year,month);
    draw_string (str, 0, size = 32, lc = {0,0,0}, lw = 4); // bottom-left

    sprintf (str, "                                           Blue: Head in Bajocian semi-confined aquifer");
    draw_string (str, 1, size = 90, lc = {0,0,1}, lw = 1);  // from top-left

    sprintf (str, "\n                                           Red:  Head in mid-Bathonian unconfined aquifer");
    draw_string (str, 1, size = 90, lc = {1,0,0}, lw = 1);  // from top-left

    sprintf (str, "Z-stretch x 25 ");
    draw_string (str, 3, size = 64, lc = {0,0,0}, lw = 1); 
  }
  sprintf(fname,"gw_levels-%4.4d.png",(int)t);  
  save (fname);
  clear();

}

/** ## Function to compute Gronde & Provence average monthly flow 'on the fly'
This would be painfull with gnuplot. The function also writes observed data in face of each simulated value for years 1965 & 1966.
*/

event monthly_discharge(t=DAYS_PER_MONTH;t+=DAYS_PER_MONTH)
{
  // Monthly flow observations for 1965 & 1966, in L/s (see )
  double Gronde_obs[24] = {580.,689.,113.,91.,44.,14.,2.9,0.5,1.7,1.6,25.,408.,455.,178.,105.,85.,89.,32.,11.,3.,0.7,90.,201.,224.};
  double Provence_obs[24] = {64.,68.,24.,32.,34.,31.,20.,21.,29.,23.,27.,79.,83.,73.,57.,53.,53.,37.,33.,31.,20.,38.,68.,86};
  double S_Gronde = 20.30; // in km2
  double S_Provence = 9.48; // in km2

  int year,month;
  QM_Gronde *= 1000./sum_dt_month; // average in L/s
  QM_Provence *= 1000./sum_dt_month; // average in L/s

  hack_date(t/oneday-DAYS_PER_MONTH/2.,& year, & month);
  if( pid()==0 && year>=1965 && year<=1966)
  {
    int ix = 12*(year-1965)+month-1;
    fprintf(flowfile,"%d %d %g %g %g %g %g %g %g %g\n",year,month,
      QM_Gronde,Gronde_obs[ix],
      QM_Provence,Provence_obs[ix],
      DAYS_PER_MONTH * 0.0864 * QM_Gronde/S_Gronde,DAYS_PER_MONTH * 0.0864 * Gronde_obs[ix]/S_Gronde,
      DAYS_PER_MONTH * 0.0864 * QM_Provence/S_Provence,DAYS_PER_MONTH * 0.0864 * Provence_obs[ix]/S_Provence);
    fflush(flowfile);
  }

  sum_dt_month = 0.;
  QM_Gronde = 0.;
  QM_Provence = 0.;
}

/**
## Animation of the solution (simulation timestep: 2 days)
![Saturated areas](./meuvaines/seepage.mp4)(width=75% )
![Groundwater head in both aquifers. Filled contours have 10-meter interval. The head in both aquifers is initially set to the topographic level, so that a warm-up period is needed in order for it to drop to a deeper, more realistic level](./meuvaines/gw_levels.mp4)(width=75% )
*/