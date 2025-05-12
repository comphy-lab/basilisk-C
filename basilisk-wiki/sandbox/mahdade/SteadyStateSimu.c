/** This is an example of a steady state simulation on the Garonne River using parallel computing and multiple inlets*/
#include "saint-venant.h"
#include "input.h"
#include "output.h"
#include "esri-binary-grids.h"
#include "terrain.h"
#include "sourceterm/manning-tilt.h"
#include "inflow-tilt.h"

#define MINLEVEL 6
#define ETAE     1e-2 // error on free surface elevation (1 cm)
#define HE       1e-2 // error on water depth (1 cm)
#define ZE       1e-2 // error on bed elevation (1 cm)
#define UE       1e-2 // 0.01 m/s
#define tmax     86400.
#define sec_per_day 86400.
#define NbInlets 6
#define Qverdun 90.
#define etaini 7.5

int LEVEL;

struct Inlet INLETS[NbInlets];
double Kmin = 25., Kmaj = 15.;
double rotation = -0.43798, xc0 = 565e3, yc0 = 6294e3;
FILE * fpcontrol;
FILE * results ;
double Qlast[NbInlets];

scalar zte[];

coord LambertToLocal(coord P)
{
   return (coord) {cos(rotation)*(P.x-xc0)-sin(rotation)*(P.y-yc0) , sin(rotation)*(P.x-xc0)+cos(rotation)*(P.y-yc0) };
}
coord LocalToLambert(coord P)
{
   return (coord) {xc0+cos(-rotation)*P.x-sin(-rotation)*P.y, yc0+sin(-rotation)*P.x+cos(-rotation)*P.y };
}

int output_flt_lambert (struct OutputFLT p); // prototype, see definition at the end of the file

int main (int argc, char * argv[])
{
  G = 9.81;

/** The boundary conditions*/
  u.n[top] = neumann(0.);
  u.t[top] = dirichlet(0.);
  h[top] = neumann(0.);
  zb[top] = neumann(0.);
  eta[top] = neumann(0.);
  nmanning[top] = neumann(0.);
  zte[top] = neumann(0.);

  run();

  for(int in=0;in<NbInlets;in++)
   (void) Deallocate_Inlet(INLETS+in);
  
  if(pid()==0.)
   fclose(fpcontrol);
 
  return(0);
}

/** Initial conditions*/

event init (i = 0)
{

  LEVEL = 12; 
  L0 = 32768.;

  tilt.x = 0.;
  tilt.y = 1.103e-3;

/** Init grid /!\ local coordinates, not Lambert 93*/

  size (L0);
  X0 = -L0/2.;
  Y0 = -17068.;
  origin (X0, Y0);
  N = 1 << LEVEL;
  init_grid (1 << LEVEL);

  terrain (zb, "./garonne_2m_local", NULL);
  boundary({zb});

  terrain (zte, "./enveloppe_4m_local", NULL);
  boundary({zte});

  char fltfile[200];
  sprintf(fltfile,"./topo_ini.flt");
  output_flt(zb,fltfile);

  coord LeftBank,RightBank;

/** Discharge inlets, numbered from upstream to downstream from 0 to (NbInlets-1)
    The coordinates of the cross-sections boundaries are defined in Lambert 93
    hen transformed into local coordinates  */

   // Inlet #0 : Garonne in Pont des Catalans
   int inj=0;
   LeftBank = LambertToLocal((coord) {573025.,6279345.});   
   RightBank = LambertToLocal((coord) {573085.,6279655.});
   (void) Init_Inlet(&INLETS[inj],LeftBank,RightBank);
   INLETS[inj].ID = inj;
   INLETS[inj].with_velocity = true;
   printf("Left bank : (%g,%g) ; Right bank : (%g,%g)\n",INLETS[inj].Segment[0].x,INLETS[inj].Segment[0].y,INLETS[inj].Segment[1].x,INLETS[inj].Segment[1].y);

   // Inlet #1 : Touch (528 km2)
   inj++;
   LeftBank = LambertToLocal((coord){570317.,6281695.});   
   RightBank = LambertToLocal((coord){570352.,6281296.});
   (void) Init_Inlet(&INLETS[inj],LeftBank,RightBank);
   INLETS[inj].ID = inj;
   INLETS[inj].with_velocity = false;
   printf("Left bank : (%g,%g) ; Right bank : (%g,%g)\n",INLETS[inj].Segment[0].x,INLETS[inj].Segment[0].y,INLETS[inj].Segment[1].x,INLETS[inj].Segment[1].y);

   // Inlet #2 : Aussonnelle (190 km2)
   inj++;
   LeftBank = LambertToLocal((coord){567088.,6289477.});   
   RightBank = LambertToLocal((coord){567355.,6289178.});
   (void) Init_Inlet(&INLETS[inj],LeftBank,RightBank);
   INLETS[inj].ID = inj;
   INLETS[inj].with_velocity = false;
   printf("Left bank : (%g,%g) ; Right bank : (%g,%g)\n",INLETS[inj].Segment[0].x,INLETS[inj].Segment[0].y,INLETS[inj].Segment[1].x,INLETS[inj].Segment[1].y);

   // Inlet #3 : Hers-Mort (1524 km2)
   inj++;
   LeftBank = LambertToLocal((coord){567779.,6297539.});   
   RightBank = LambertToLocal((coord){567937.,6297906.});
   (void) Init_Inlet(&INLETS[inj],LeftBank,RightBank);
   INLETS[inj].ID = inj;
   INLETS[inj].with_velocity = false;
   printf("Left bank : (%g,%g) ; Right bank : (%g,%g)\n",INLETS[inj].Segment[0].x,INLETS[inj].Segment[0].y,INLETS[inj].Segment[1].x,INLETS[inj].Segment[1].y);

   // Inlet #4 : Save (1146 km2)
   inj++;
   LeftBank = LambertToLocal((coord){562061.,6298072.});   
   RightBank = LambertToLocal((coord){562447.,6298177.});
   (void) Init_Inlet(&INLETS[inj],LeftBank,RightBank);
   INLETS[inj].ID = inj;
   INLETS[inj].with_velocity = false;
   printf("Left bank : (%g,%g) ; Right bank : (%g,%g)\n",INLETS[inj].Segment[0].x,INLETS[inj].Segment[0].y,INLETS[inj].Segment[1].x,INLETS[inj].Segment[1].y);

   // Inlet #5 : Saint-Pierre (135 km2)
   inj++;
   LeftBank = LambertToLocal((coord){559482.,6302451.});   
   RightBank = LambertToLocal((coord){559759.,6302163.});
   (void) Init_Inlet(&INLETS[inj],LeftBank,RightBank);
   INLETS[inj].ID = inj;
   INLETS[inj].with_velocity = false;
   printf("Left bank : (%g,%g) ; Right bank : (%g,%g)\n",INLETS[inj].Segment[0].x,INLETS[inj].Segment[0].y,INLETS[inj].Segment[1].x,INLETS[inj].Segment[1].y);

   DT = 1.0;

  scalar zmin = zb.dmin;
  foreach(){
    nmanning[] = zmin[]<(zte[]+3.) ? 1./Kmin : 1./Kmaj;
    eta[] = fmax(etaini,zte[]+0.4);
    h[] = eta[] > zb[] ? eta[]-zb[] : 0.;
    eta[] = zb[] + h[];
    u.x[] = 0.;
    u.y[] = 0.;
  }

  boundary ((scalar *){u});
  boundary({h,eta,nmanning});

  for(int inj=0;inj<NbInlets;inj++)
    Qlast[inj] = 0.;

  if(pid()==0)
   fpcontrol = fopen("./inflow_control_6inj.txt","w");
    results = fopen("./outflow_results.txt","w");
}

double etainj;
double Qinj[NbInlets] = {0.925*Qverdun,0.15*0.075*Qverdun,0.05*0.075*Qverdun,0.43*0.075*Qverdun,0.33*0.075*Qverdun,0.04*0.075*Qverdun};
double Tinj[NbInlets] = {0.,3600.,3600.,3600.,3600.,3600.};

event inflows (i++)
{
  for(int inj=0;inj<NbInlets;inj++){
    if(t>=Tinj[inj]){
      double Qtarg = Qlast[inj]+(1.-exp((Tinj[inj]-t)/180.))*(Qinj[inj]-Qlast[inj]);
      etainj = eta_b ( Qtarg*INLETS[inj].Controller.fQ, &INLETS[inj]);
      if(etainj>18.0)etainj=18.0;
  
      (void) set_inlet_fields(&INLETS[inj],etainj);
    }
  }

  foreach()
   eta[] = zb[]+h[];

  boundary ((scalar *){u});
  boundary({h,eta});

}

double etas;

event inflow_control(i++)
{
   coord Section[2];
   double Qmes[NbInlets];

   for(int j=0;j<NbInlets;j++){
     
     // The discharge is controlled slightly downstream of the inlet sections (distance dmes)
     double dmes = sqrt(2.)*(L0/N);
     Section[0] = (coord) { INLETS[j].Segment[0].x+dmes*INLETS[j].vec_n.x , INLETS[j].Segment[0].y+dmes*INLETS[j].vec_n.y };
     Section[1] = (coord) { INLETS[j].Segment[1].x+dmes*INLETS[j].vec_n.x , INLETS[j].Segment[1].y+dmes*INLETS[j].vec_n.y };
     // "Measure" discharge
     Qmes[j] = segment_flux (Section,h,u);
     if(t<=Tinj[j])Qlast[j] = Qmes[j];

     // Set PID controllers parameters
     // https://en.wikipedia.org/wiki/PID_controller
     double Ku = 0.5, Tu = dmes*1.0; 
     INLETS[j].Controller.Kp = 0.2*Ku;
     INLETS[j].Controller.Ti = 30.;
     INLETS[j].Controller.Td = 5.;
     double Qtarg = t>Tinj[j] ? Qlast[j]+(1.-exp((Tinj[j]-t)/180.))*(Qinj[j]-Qlast[j]) : 0.;

     if(Qmes[j]>0. & Qtarg>0. & fabs(log(INLETS[j].Controller.fQ))<log(10.) )        
        (void) update_PID_error(Qmes[j],Qtarg, & (INLETS[j].Controller) );
     else
	(void) Reset_PID_Controller( & (INLETS[j].Controller) );
   }

   double Vol;
   Vol = 0.;

   foreach(reduction(+:Vol))
      Vol+=sq(Delta)*h[];

   if(pid()==0){
    fprintf(fpcontrol,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",t,dt,Vol,
                    Qmes[0],INLETS[0].Controller.e,INLETS[0].Controller.E,INLETS[0].Controller.fQ , 
                    Qmes[1],INLETS[1].Controller.e,INLETS[1].Controller.E,INLETS[1].Controller.fQ , 
                    Qmes[2],INLETS[2].Controller.e,INLETS[2].Controller.E,INLETS[2].Controller.fQ , 
                    Qmes[3],INLETS[3].Controller.e,INLETS[3].Controller.E,INLETS[3].Controller.fQ , 
                    Qmes[4],INLETS[4].Controller.e,INLETS[4].Controller.E,INLETS[4].Controller.fQ , 
                    Qmes[5],INLETS[5].Controller.e,INLETS[5].Controller.E,INLETS[5].Controller.fQ );
    fflush(fpcontrol);
   }
}

event outflow(i++)
{
   
   coord Avalcoord[2];
   double Qaval, etaaval;
   
   
   Avalcoord[0] = LambertToLocal((coord) {558575.,6308295.}); 
   Avalcoord[1] = LambertToLocal((coord) {558733.,6308374.}); 
   
   
   
   Qaval = segment_flux (Avalcoord,h,u);
   etaaval = segment_etamean (Avalcoord,h,zb);
  

   if(pid()==0){
    fprintf(results,"%g %g %g %g\n",t,dt,Qaval,etaaval);
    fflush(results);
   }}


scalar l[];

event snapshot (t+=600. ; t<=tmax) {

  char fltfile[200];
  sprintf(fltfile,"./h_l%d_Q_90_Kmin_%d_Kmaj_%d_%5.5d.flt",LEVEL,(int)Kmin,(int)Kmaj,(int)(t/60.));
  output_flt(h,fltfile);

  sprintf(fltfile,"./ux_l%d_Q_90_Kmin_%d_Kmaj_%d_%5.5d.flt",LEVEL,(int)Kmin,(int)Kmaj,(int)(t/60.));
  output_flt(u.x,fltfile);

  sprintf(fltfile,"./uy_l%d_Q_90_Kmin_%d_Kmaj_%d_%5.5d.flt",LEVEL,(int)Kmin,(int)Kmaj,(int)(t/60.));
  output_flt(u.y,fltfile);

  foreach()
    l[] = level;

  sprintf(fltfile,"./level_l%d_Q_90_Kmin_%d_Kmaj_%d_%5.5d.flt",LEVEL,(int)Kmin,(int)Kmaj,(int)(t/60.));
  output_flt(l,fltfile);
}

// Adaptivity

scalar absu[];

int adapt() {
#if TREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  foreach()
    absu[]= h[] > dry ? sqrt(u.y[]*u.y[]+u.x[]*u.x[]) : 0.;

  /**
  We can now use wavelet adaptation on the list of scalars *{Î·,hmax}*
  with thresholds *{ETAE,HMAXE}*. The compiler is not clever enough yet
  and needs to be told explicitly that this is a list of *double*s,
  hence the *(double[])*
  [type casting](http://en.wikipedia.org/wiki/Type_conversion). 
  
  The function then returns the number of cells refined.*/

  astats s = adapt_wavelet ({h,absu}, (double[]){HE,UE},
			    LEVEL, MINLEVEL);

  // Refine inlets

  for(int in=0;in<NbInlets;in++){

    double xA = INLETS[in].Segment[0].x , xB = INLETS[in].Segment[1].x, yA = INLETS[in].Segment[0].y , yB = INLETS[in].Segment[1].y;
    coord vec_n = {yA - yB, xB - xA};
    coord vec_t = {xB - xA, yB - yA};
    double norm = sqrt(sq(vec_n.x) + sq(vec_n.y));
    assert (norm > 0.);
    vec_n.x = vec_n.x/norm + 1e-6, vec_n.y = vec_n.y/norm - 1.5e-6;
    vec_t.x = vec_t.x/norm, vec_t.y = vec_t.y/norm;

//    printf(" Refining inlet #%d, A = (%g,%g), vec_n = (%g,%g)\n",in,xA,yA,vec_n.x,vec_n.y);

     refine (level < LEVEL
            && fabs((x-xA)*vec_n.x+(y-yA)*vec_n.y) <= 2.83*(L0/(1 << LEVEL))
            && ((x-xA)*vec_t.x+(y-yA)*vec_t.y) >= 0.
            && ((x-xA)*vec_t.x+(y-yA)*vec_t.y) <= norm
          );
  }

  scalar zmin = zb.dmin;
  foreach(){
    nmanning[] = zmin[]<(zte[]+3.) ? 1./Kmin : 1./Kmaj;
    l[] = level;    
  }

  boundary(all);
  
//  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;

#else // Cartesian
  return 0;
#endif
}

// And finally we set the event that will apply our *adapt()* function at every timestep.

event do_adapt (i++) adapt();

////////////////////////////////////////////////////////////////////////////////////

int output_flt_lambert (struct OutputFLT p)
{
  FILE * flt , * hdr ;
  scalar input = p.f;
  char * pch;
  char buffer[100];
  double * data_double;
  coord current_point;
  char hdrfile[200]; 

  // default values
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;

  data_double = (double *)malloc(nx*ny*sizeof(double)); 
  for(int i=0;i<(nx*ny);i++)
    data_double[i] = nodata;

  // data
  for (int j = ny-1; j >= 0; j--) {
    double yp93 = Delta*j + p.box[0][1] + Delta/2.;    // center of pixel
    for (int i = 0; i < nx; i++) {
      double xp93 = Delta*i + p.box[0][0] + Delta/2., v;
      float vf;

      current_point.x = xp93;
      current_point.y = yp93;
      current_point = LambertToLocal(current_point);
      double xp = current_point.x; // local x
      double yp = current_point.y; // local y

      if (p.mask.i) { // masking
    if (p.linear) {
      double m = interpolate (p.mask, xp, yp);
      if (m < 0.)
        v = nodata;
      else
            v = interpolate (p.f, xp, yp);
    }
    else {
      Point point = locate (xp, yp);
      if (point.level < 0 || val(p.mask) < 0.)
        v = nodata;
      else
        v = val(p.f);
    }
      }
      else if (p.linear){
    Point point = locate (xp, yp);
        v = point.level >= 0 ? interpolate (p.f, xp, yp) : nodata;
      }
      else {
    Point point = locate (xp, yp);
    v = point.level >= 0 ? val(p.f) : nodata;
      }
      if (v == nodata)
    data_double[i+(ny-1-j)*nx] = -9999.;
      else
    data_double[i+(ny-1-j)*nx] = v;

    }
  }

    if (pid() == 0) { // master writes .FLT and .HDR files
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, data_double, nx*ny, MPI_DOUBLE, MPI_MAX, 0,
        MPI_COMM_WORLD);
@endif

// Set path to header file
//  hdrfile = (char *)malloc(sizeof(char)*strlen(p.file));
  strcpy(hdrfile,p.file);
 
  pch = hdrfile + strlen(hdrfile)-3;
//  strncpy (pch,"hdr",3);
  memcpy (pch,"hdr",3);
 
  if(!(hdr = fopen (hdrfile, "w")))
  {
    printf("Failed to open header file %s\n",hdrfile);
    return -1;
  }

  // header
  fprintf (hdr, "ncols          %d\n", nx);
  fprintf (hdr, "nrows          %d\n", ny);
  fprintf (hdr, "xllcorner      %.8g\n", p.box[0][0]);
  fprintf (hdr, "yllcorner      %.8g\n", p.box[0][1]);
  fprintf (hdr, "cellsize       %.8g\n", Delta);
  fprintf (hdr, "nodata_value   -9999\n");
  fprintf (hdr, "byteorder   LSBFIRST\n");

  fflush(hdr);
  fclose(hdr);
//  free(hdrfile);

  // open flt file

  bool opened = false;
  if( !(flt = fopen (p.file, "wb")) ) {
      perror (p.file);
      exit(1);
  }
  else
      opened = true;

  for(int i=0;i<(nx*ny);i++)
  {
    float vf = (float)data_double[i];
    fwrite(&vf,1,sizeof(float),flt);
  }
  fclose(flt);
}
@if _MPI
  else // slave does not write anything
    MPI_Reduce (data_double, NULL, nx*ny, MPI_DOUBLE, MPI_MAX, 0,
        MPI_COMM_WORLD);
@endif

  free(data_double);
  return(0);

}