
#define USE_CONJUGATE_HEAT 2   //0:no conjugate heat-transfer 1:explicit coupling 2:implicit coupling

#define INT_USE_UF
#define CONSISTENTPHASE2
#define SHIFT_TO_GAS
#define INIT_TEMP


#include "axi.h"
#include "mypoisson.h"   //I don't understant why, but I has to add this line here to compile the code
#include "centered-evaporation.h"
#include "centered-doubled.h"
#include "two-phase.h"
#include "mytension.h"
#include "evaporation.h"
#include "temperature-gradient.h"
#include "view.h"
#if USE_MY_SOLID
#include "mysolid.h"
#endif //USE_MY_SOLID

#if USE_CONTACT_ANGLE
vector h[];
double theta0 = 50.;
h.t[left] = contact_angle (theta0*pi/180.);
#endif

double lambda1, lambda2, cp1, cp2, dhev;
double TL0, TG0, TIntVal, Tsat, Tbulk;
#if USE_CONJUGATE_HEAT
double rhos, lambdas, cps;
#endif //USE_CONJUGATE_HEAT

double rhos2, lambdas2, cps2;
double thickness_heater;
scalar vof_heater[]; //volume fraction of heater


int maxlevel = 11, minlevel = 6;
int dump_num = 0;

const double uemax = 1e-2;
const double Temax = 1e-2;
const int num_refine = 1;

const double tshift = 0.0;
double Lsize = 7.0e-3;
const double LsizeLarge = 7.0e-3;
const double LsizeSmall = 2.3e-3;
double R0 = 10e-6;

void solidSetup1()
{

  double delta = LsizeSmall / (1024.0);
  int ncs = ceil(1.0e-3 / delta);
  double lenx = ncs * delta;
  ncs = ceil(1.0e-3 / delta);
  double leny = ncs * delta;
  if(pid() == 0.0)
  {
    printf("solid_len_x: %g\n",lenx);
    printf("solid_len_y: %g\n",leny);
  }
  SOLID_LEN_x = lenx;
  SOLID_LEN_y = leny;
}

void solidSetup2()
{
  double delta = LsizeLarge / (1024.0);
  int ncs = ceil(1.0e-3 / delta);
  double lenx = ncs * delta;
  if(pid() == 0.0)
  {
    printf("solid_len_x: %g\n",lenx);
  }
  SOLID_LEN_x = lenx;
  SOLID_LEN_y = 0.0;//0.00099951171875;
}

#include "myfunctions.h"

int Ncell = 64;
int step = 0;
int main (void) {
  for (maxlevel = 13; maxlevel <= 13; maxlevel++)
  {
    step = 0;
    printf("maxlevel:%d begins\n", maxlevel);
    // first step: output grid
    Ncell = 1 << minlevel;
    Lsize = LsizeSmall;
    init_grid(Ncell);
    solidSetup1();
    size(Lsize);
    origin(-SOLID_LEN_x, -SOLID_LEN_y);
    run();

    // second step: read grid and do the interpolation
    step++;
    Lsize = LsizeLarge;
    solidSetup2();
    size(Lsize);
    origin(-SOLID_LEN_x, -SOLID_LEN_y);
    run();

    // thrid step: do the initialization and output the dumpfile
    step++;
    Ncell = 1 << minlevel;
    Lsize = LsizeSmall;
    init_grid(Ncell);
    solidSetup1();
    size(Lsize);
    origin(-SOLID_LEN_x, -SOLID_LEN_y);
    run();

    printf("maxlevel:%d ends\n", maxlevel);
  }
}

event defaults(i = 0)
{
  printf("setflags\n");
  setSolidFlag();
  TS.boundarySolid_y = boundarySolidNeumman_y;
}

#define circle(x, xc, y, yc, R) (sq(R) - sq(x - xc) - sq(y - yc))

void writeGridPosition()
{
  char filename[30];
  sprintf(filename, "./data/gridT-%d.dat", maxlevel);
  FILE * fp = fopen (filename, "w");
  int counter = 0;
  foreach()
  {
    counter++;
    fprintf (fp, "%.13lf %.13lf\n", x, y);
    fflush(fp);
  }
  printf("%d cells output\n", counter);
  fclose(fp);
  return;
}

void interpolateTemp()
{
  char name_grid[30];
  char name_ini[30];

  sprintf(name_grid, "./data/gridT-%d.dat", maxlevel);
  sprintf(name_ini, "./data/iniT-%d.dat", maxlevel);

  FILE *input = fopen(name_grid, "r");
  FILE *output = fopen(name_ini, "w");
  double x, y, interp_val;

  if (input == NULL || output == NULL) {
      printf("Error opening file!\n");
      return;
  }

  scalar T1[], T2[];
  foreach()
  {
    T1[] = TS[];
    T2[] = TS[];
  }

  foreach()
  {
    if(is_solid[] == 1. && is_solid[1] == 0.)
    {
      T2[] = TS[1];
    }

    if(is_solid[] == 0. && is_solid[-1] == 1.)
    {
      T1[] = TS[-1];
    }
  }

  int counter = 0;
  while (fscanf(input, "%lf %lf", &x, &y) == 2) {
    counter ++;
      //printf("%dth \n", counter++);
    Point tmp = locate(x, y);
    scalar Ti = x >= 0.0 ? T2 : T1;
    interp_val = interpolate(Ti, x, y); 
    interp_val = y < 0.0 ? 0.0 : interp_val;
    if(y > 0.0 && interp_val < 0.0)
    {
      printf("wrong cell!\n");
    }
    fprintf(output, "%.11lf %.11lf %g\n", x, y, interp_val); 
    fflush(output);
  }
  printf("number of cells: %d\n", counter);
  fclose(input);
  fclose(output);

  scalar Tout[];

  foreach()
  {
    Tout[] = TS[];
  }

  boundarySolidNeummanNoauto(Tout);
  char Tname[80];
  sprintf(Tname, "./data/T-wall-%d.dat", maxlevel);
  outputWallTemperature(Tname, Tout);
  return;
}

void outputDumpTemp()
{
  char name_ini[30];
  sprintf(name_ini, "./data/iniT-%d.dat", maxlevel);

  FILE *input = fopen(name_ini, "r");
  int counter = 0;
  foreach ()
  {
    counter++;
    double xx, yy, TT;
    if(fscanf(input, "%lf %lf %lf", &xx, &yy, &TT) != 3)
    {
      printf("read error for %dth \n", counter);
    }
    if(TT != 0.0)
    {
      int pause = 0;
    }
    TS[] = TT;
  }
  printf("%d cells initialized\n", counter);
  boundary({TS});
  char dumpname[80];
  sprintf(dumpname, "./data/dumpTini-%d", maxlevel);
  dump(dumpname, list = {TS, is_solid});
  scalar Tout[];

  foreach()
  {
    Tout[] = TS[];
  }

  boundarySolidNeummanNoauto(Tout);
  char Tname[80];
  sprintf(Tname, "./data/T-wall_inter-%d.dat", maxlevel);
  outputWallTemperature(Tname, Tout);
  return;
}

void refineRegions()
{
  scalar solid_refine[];
  scalar solid_refiney[];
  scalar f_refine[];
  foreach ()
  {
    solid_refine[] = is_solid[];
    solid_refiney[] = is_solid_y[];
    f_refine[] = 0.0;
  }

  foreach ()
  {
    if (is_solid[] == 1.0 && is_solid[1] == 0.0)
    {
      solid_refine[] = 0.5;
    }
    if (is_solid[] == 0.0 && is_solid[-1] == 1.0)
    {
      solid_refine[] = 0.5;
    }
    if (is_solid_y[] == 0.0 && is_solid_y[0, -1] == 1.0)
    {
      solid_refiney[] = 0.5;
    }
    if (is_solid_y[] == 1.0 && is_solid_y[0, 1] == 0.0)
    {
      solid_refiney[] = 0.5;
    }
  }

  setJumpVarAdaptation(solid_refine, num_refine);
  adapt_wavelet({solid_refine, solid_refiney}, (double[]){uemax, uemax}, maxlevel, minlevel);
}

event init (i = 0) {
  printf("init\n");
  if (step == 0)
  {
    for(int ilevel = 6; ilevel < maxlevel; ilevel ++)
      refineRegions();
    
    writeGridPosition();
    return 1;
  }
  else if(step == 1)
  {
    restore (file = "T_input", list = {TS});
    setSolidFlag();
    interpolateTemp();
    return 1;
  }
  else 
  {
    for(int ilevel = 6; ilevel < maxlevel; ilevel ++)
      refineRegions();
    
    outputDumpTemp();
    return 1;
  }
}