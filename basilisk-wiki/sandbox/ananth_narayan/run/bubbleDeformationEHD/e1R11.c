/**
Compilation: ```CC99='mpicc -std=c99 -D_GNU_SOURCE' qcc -Wall -O2 -D_MPI=1 e1R11.c -o e1R11 -lm```

Run: ```mpirun -np 8 ./e1R11```

** OR **

```CC='mpicc -std=c99 -D_GNU_SOURCE -D_MPI=8' make e1R11.tst```
*/

/**
# Introduction

When a liquid droplet suspended in an immiscible fluid is subjected to an electric field, the mismatch of the physical properties between the two fluids gives rise to the electrical stress on the interface, which is likely to deform the droplet. This was first intuitively understood by Taylor. Later Taylor and Melcher formulated Leaky Dielectric Model (LDM). Realistically, the liquids in the aforementioned  applications are considered to be leaky dielectrics with finite electrical conductivity, and the electrohydrodynamics (EHD) of the droplet can be described by the leaky dielectric model.

However, the same conclusion cannot be derived when it comes to the transient deformation of the droplets and temporal evolution of the fluid flow toward the steady state. Since such understanding could be of great relevance in microfluidic applications or electric field-driven enhancement of heat and mass transfer, several studies have been reported to reveal the various characteristics of the transient electrohydrodynamics of liquid droplets.


# Code setup
Including header for axisymmetric, Navier-Stokes, EHD leaky dielectric, stress terms, VoF and surface tension
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "ehd/implicit.h"
#include "ehd/stress.h"
#include "vof.h"
#include "tension.h"

/**
Initialization of variables for calculation and post processing. 
*/

face vector muv[];
scalar f[], * interfaces = {f};
#define LEVEL0 7
#define LEVELMAX 10
int LEVEL;
double Radius = 1., R = 0.1, Q = 0.1, Lambda = 1, CaE, refD;
int cases;
char csvName[200];
FILE * fpD;

/**
For performing parametric study, the data's has been setup in structure and looped in `int main()`. For details to these test data, refer [Lac2007](#Lac2007) Fig 5, 7 and 10. The `refD` is used to determine the error (%) for each case against `caseNo`. Fig 10 refers to test data similar to [Jiang2020](#Jiang2020) which is the article commonly used for validation of droplet deformation under electric field.
*/
typedef struct {
  int caseNo;
  double R;
  double Q;
  double Lambda;
  double CaE;
  double refD;
} CaseData;

static CaseData casesList[] = {
{1,0.1,0.1,1,0.05,0.0235},
{2,0.1,0.1,1,0.1,0.0512},
{3,0.1,0.1,1,0.15,0.0833},
{4,0.1,0.1,1,0.2,0.122},
{5,0.1,0.1,1,0.25,0.171},
#if 0 // set to 1 to run all the cases
{6,0.1,0.1,1,0.28,0.211},
{7,0.1,0.1,1,0.3,0.246},
{8,0.1,0.1,1,0.31,0.267},
{9,0.1,0.1,1,0.32,0.294},
{10,0.1,0.1,1,0.33,0.33},
{11,0.1,0.1,1,0.34,0.412},
{12,0.1,5,1,0.05,0.0189},
{12,0.1,5,1,0.05,0.0189},
{13,0.1,5,1,0.1,0.0413},
{14,0.1,5,1,0.15,0.0671},
{15,0.1,5,1,0.2,0.0964},
{16,0.1,5,1,0.25,0.129},
{17,0.1,5,1,0.3,0.17},
{18,0.1,5,1,0.35,0.225},
{19,0.1,5,1,0.4,0.306},
{20,0.1,5,1,0.45,0.506},
{21,0.1,5,1,0.5,0.745},
{22,0.1,5,1,0.55,0.818},
{23,0.1,5,1,0.6,0.862},
{24,0.1,5,1,0.65,0.892},
{25,0.1,5,1,0.7,0.912},
{26,0.1,1.37,1,0.025,0.0115},
{27,0.1,1.37,1,0.05,0.023},
{28,0.1,1.37,1,0.1,0.0509},
{29,0.1,1.37,1,0.15,0.0788},
{30,0.1,1.37,1,0.2,0.115},
{31,0.1,1.37,1,0.25,0.159},
{32,0.1,1.37,1,0.3,0.22},
{33,0.1,1.37,1,0.32,0.255},
{34,0.1,1.37,1,0.35,0.337},
{35,0.1,1.37,1,0.36,0.399},
{36,0.1,1.37,1,0.38,0.883}
#endif
};
#define NCASES (sizeof(casesList)/sizeof(casesList[0]))

/**
`deformationFactor(f)` computes the deformation factor, $D = \dfrac{L-B}{L+B}$ which is modified version of `output_facets(f)`. As this is axisymmetric case, `y` coordinate is doubled.
 */
double deformationFactor (scalar c, face vector s = {{-1}})
{
  double xmin = HUGE, xmax = -HUGE;
  double ymax = 0.;

  foreach(reduction(min:xmin)
          reduction(max:xmax)
          reduction(max:ymax)) {

    if (c[] > 1e-6 && c[] < 1. - 1e-6) {

      coord n = facet_normal(point, c, s);
      double alpha = plane_alpha(c[], n);

      coord segment[2];

      if (facets(n, alpha, segment) == 2) {

        for (int k = 0; k < 2; k++) {

          double xf = x + segment[k].x*Delta;
          double yf = y + segment[k].y*Delta;

          xmin = min(xmin, xf);
          xmax = max(xmax, xf);
          ymax = max(ymax, yf);
        }
      }
    }
  }

  double L = xmax - xmin;
  double B = 2.*ymax;

  return (L - B)/(L + B);
}

/**
`int main ()` contains the loop for running cases for different `CaE` numbers. The domain is $16R \times 16R$ as suggested in [Jiang2020](#Jiang2020) with the droplet located at $(0,0)$. We export the consolidated data of deformation factor to a separate csv file `bubbleDeformationEHD.csv`
*/
int main ()
{
  system ("rm -r Validation");
  system ("mkdir Validation");
  X0 = (-8.)[0];
  L0 = (16.)[0];
  DT = 1e-1[0];
  TOLERANCE = 1e-7[0];
  for (int k = 0; k < NCASES; k++) 
  {
    LEVEL = LEVEL0;
    init_grid(1 << LEVEL);
    f.sigma = 1.; 
    mu = muv;  
    R = casesList[k].R;
    Q = casesList[k].Q;
    Lambda = casesList[k].Lambda;
    cases = casesList[k].caseNo;
    CaE   = casesList[k].CaE;
    refD = casesList[k].refD;
    fprintf(stderr,"\nRunning Case %d with CaE = %g\n",cases, CaE);
    run();
  } 
  system(
    "cd Validation/ && "
    "echo 'CaE,D,refD' > bubbleDeformationEHD.csv; "
    "for file in $(ls e1R11v*.csv | sort -V); do "
    "CaE=$(awk -F',' 'NR==2 {print $7}' \"$file\"); "
    "refD=$(awk -F',' 'NR==2 {print $8}' \"$file\"); "
    "D=$(tail -n 1 \"$file\" | awk -F',' '{print $2}'); "
    "echo \"$CaE,$D,$refD\"; "
    "done >> bubbleDeformationEHD.csv"
  );
}
/**
Initializing the domain and the setup with variable `phi[]` where the electric field vector points along $-x$. The whole setup is non-dimensionalized such that $E_\infty \sim \sqrt{\text{Ca}_E}$.
 */

scalar uN[];
event init (t = 0)
{
  fraction (f, sq(Radius)-(sq(x)+sq(y))); 
  foreach()
    phi[] = sqrt(CaE)*(x+L0/2.);
  if (pid() == 0) {
      sprintf(csvName, "Validation/e1R11v%d.csv",cases);
      fpD = fopen(csvName, "w");
      fprintf(fpD, "Time,D,Err %%,R,Q,Lambda,CaE,refD\n");
      fprintf(fpD, "0,0,0,%g, %g, %g, %g, %g\n", R, Q, Lambda, CaE, refD);
      fflush(fpD);
  }  
}
phi[top]   = dirichlet(sqrt(CaE)*(x+L0/2.));
phi[left]  = dirichlet(sqrt(CaE)*(x+L0/2.));
phi[right] = dirichlet(sqrt(CaE)*(x+L0/2.));
uf.n[left] = 0.;
uf.n[right] = 0.;
uf.n[top] = 0.;
uf.n[bottom] = 0.;

/**
The properties of ion conductivity `K`, permittivity `epsilon` and viscosity `mu = muv` is written as arithmetic interpolation as there are two phases. Some studies uses harmonic interpolation which may show better results. But I found that arithmetic interpolation is sufficient to capture the deformation factor.
 */
event properties (i++)
{
    foreach_face()
    {
        double Tx = (f[] + f[-1])/2.;
        epsilon.x[] = (Tx*Q + (1.-Tx))*fm.x[];
        K.x[] = (Tx + (1.-Tx)*R)*fm.x[];
        muv.x[] = (Tx*Lambda + (1.-Tx))*fm.x[];
    }
}
/**
The whole system is being monitored for steady state convergence through velocity and deformation factor. Thus the adpative mesh transistions from coarse `LEVEL0` to max of LEVELMAX as in `return (LEVEL++ == LEVELMAX)`

 */

event initErr (i = 0)
{
  foreach()
  {
    uN[] = u.x[];
  }
  fprintf(stderr, "Max potential difference is %g \n",statsf(phi).max);
}
double Dprev = 0.;
event monitor (i+=20; t<=500.)
{
    double du = change (u.x, uN);
    double D = deformationFactor(f);
    double err = fabs(D - Dprev);
    Dprev = D;
    if (pid() == 0) {
        fprintf(fpD, "%g,%g,%g\n", t, D, (D - refD)*100/refD);
        fflush(fpD);
    } 
    if (i > 0 && du < 5e-6 && err < 1e-6)
    {
      return (LEVEL++ == LEVELMAX); 
    }   
}
/**
Further adpative grid was utilized to capture the bubble dynamics
 */
event adapt (i += 5) 
{
  adapt_wavelet ({f, rhoe, phi, u.x, u.y}, {1e-6,1e-5,1e-4,1e-5,1e-5}, maxlevel = LEVEL);
}

event stop (t = end)
{
  if (pid() == 0) 
    fclose(fpD);
}

/**

~~~gnuplot Bubble deformation under electric field
reset
set datafile separator ","

set xlabel "Ca_E"
set ylabel "D"
set size square
set key top left
set grid

plot \
    "Validation/bubbleDeformationEHD.csv" using 1:3 with points pt 6 ps 1 title "Benchmark (Lac2007)", \
    "Validation/bubbleDeformationEHD.csv" using 1:2 with points pt 5 title "Present"
~~~

~~~bib
@Article{Lac2007,
  author    = {Lac, Etienne and Homsy, G. M.},
  journal   = {Journal of Fluid Mechanics},
  title     = {Axisymmetric deformation and stability of a viscous drop in a steady electric field},
  year      = {2007},
  issn      = {1469-7645},
  month     = oct,
  pages     = {239--264},
  volume    = {590},
  comment   = {Fig 5 can be well validated with numerical},
  doi       = {10.1017/s0022112007007999},
  file      = {:PDF/Lac2007.pdf:PDF},
  keywords  = {BubbleDynamics,Validation},
  priority  = {prio1},
  publisher = {Cambridge University Press (CUP)},
}
@Article{Jiang2020,
  author     = {Jiang, Zhengwei and Gan, Yunhua and Luo, Yanlai},
  journal    = {Physics of Fluids},
  title      = {Effect of viscosity ratio on the dynamic response of droplet deformation under a steady electric field},
  year       = {2020},
  issn       = {1089-7666},
  month      = may,
  number     = {5},
  volume     = {32},
  comment    = {Validation for bubble under EHD},
  doi        = {10.1063/5.0003449},
  file       = {:PDF/Jiang2020.pdf:PDF},
  keywords   = {Validation,BubbleDynamics},
  priority   = {prio1},
  publisher  = {AIP Publishing},
  readstatus = {skimmed},
}

~~~
*/