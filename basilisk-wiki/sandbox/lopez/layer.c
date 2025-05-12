/**
# Singularities in the thin 2D film equation

In a coming work of Dallaston et al.,
2021 (under revision in JFM) it is
studied the motion of a thin film of liquid under surface forces
with origin in the solid substrate.

The nonlinear lubrication model used in the work is formed by the
following equations:

$$ h_t = \nabla \cdot \left[h^n \nabla p - h^{n-m-1} \nabla h \right] $$ 
$$ p = -\nabla^2 h $$ 

where $h$ is the height of the film and $p$ is the 2D pressure
distribution. $n$ and $m$ are arbitrary exponents.

## Numerical method

The discretization is of second order in space and of first order in
time (although the code is prepared to increase it to second
order). The resulting nonlinear set of equations is solved iteratively
with a Newton-Raphson (N-R) scheme, 

$$
p = p^* + \delta p \quad h = h^* + \delta h
$$

where $(p^*, h^*), (\delta p, \delta h)$ and $(p, h)$ are the trial,
correction and converged solution, respectively.  
*/

#include "poisson.h"
#include "run.h"
#include "view.h"

/**
An struct is created to facilitate the back and forth of data through
the different functions.*/

struct Layer {
  scalar f;
  double dt;
  double m;
  double n;
  face vector D;
  face vector beta;
  face vector gamma;
  scalar lambda;
  double tolerance;
};

/** 
The linearization of the equations leads to the following set, 

$$
p + h_{xx} + h_{yy} = 0
$$
$$
(\alpha p_x)_x + (\alpha p_y)_y + (\beta h_x)_x + (\beta h_y)_y +
(\gamma_x h)_x + (\gamma_y h)_y = h_t + (\gamma_x h^*)_x +
(\gamma_y h^*)_y
$$

with
$$
\alpha = (h^*)^n, \: \beta = -(h^*)^{n-m-1} \: \text{and} \:
\gamma_{x,y} = n(h^*)^{n-1} p^*_{x,y} -(n-m-1)(h^*)^{n-m-2}  h^*_{x,y}.
$$

which is solved with the iterative Multigrid (MG) scheme integrated in
Basilisk. The MG scheme requires the computation of the residual of the
linear equation... */

static double residual_layer (scalar * al, scalar * bl,
			      scalar * resl, void * data)
{
  scalar h = al[0], b = bl[0], res = resl[0];
  struct Layer * p = (struct Layer *) data;
  (const) face vector alpha = p->D;
  (const) face vector beta = p->beta;
  (const) face vector gamma = p->gamma;
  (const) scalar lambda = p->lambda;
  double maxres = 0.;
  
  face vector g[], hf[];

  foreach_face()
    hf.x[] = face_gradient_x (h, 0);
  boundary_flux ({hf});

  scalar pres[];
  foreach() {
    pres[] = 0.;
    foreach_dimension()
      pres[] -= (hf.x[1]- hf.x[0])/Delta;
  }
  boundary({pres});
  
  foreach_face() {
    g.x[] = alpha.x[]*face_gradient_x (pres, 0);
    hf.x[] = beta.x[]*face_gradient_x (h, 0) +
      gamma.x[]*(h[]+h[-1])/2.;
  }
  boundary_flux ({g, hf});
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*h[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta + (hf.x[1] - hf.x[])/Delta;

    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  boundary (resl);
  return maxres;
}

/** 
...together with the subsequent relaxation toward the solution of the
linear system. */

static void relax_layer (scalar * al, scalar * bl, int l, void * data)
{
  scalar h = al[0], b = bl[0];
  struct Layer * p = (struct Layer *) data;
  (const) face vector alpha = p->D;
  (const) face vector beta = p->beta;
  (const) face vector gamma = p->gamma;
  (const) scalar lambda = p->lambda;

  
#if JACOBI
  scalar c[];
#else
  scalar c = h;
#endif

  foreach_level_or_leaf (l) {
    double n =  sq(sq(Delta))*b[], d =  lambda[]*sq(sq(Delta));
    d -= 5*(alpha.x[] + alpha.x[1] + alpha.y[] + alpha.y[1]);
    d += 0.5*Delta*sq(Delta)*(gamma.x[1]-gamma.x[] + gamma.y[1]-gamma.y[]);
    foreach_dimension() {
      n -= (beta.x[1]*h[1] + beta.x[]*h[-1])*sq(Delta);
      d -= (beta.x[1] + beta.x[])*sq(Delta);
    }

    n -= alpha.x[]*(-h[-2] - h[-1,-1] + 5*h[-1] - h[-1,1] + h[0,-1] + h[0,1] + h[1]) +
      alpha.x[1]*(h[-1] + h[0,-1] + h[0,1] - h[1,-1] + 5*h[1] - h[1,1] - h[2]) +
      alpha.y[]*(-h[-1,-1] + h[-1] - h[0,-2] + 5*h[0,-1] + h[0,1] - h[1,-1] + h[1]) +
      alpha.y[1]*(h[-1] - h[-1,1] + h[0,-1] + 5*h[0,1] - h[0,2] + h[1] - h[1,1]);

    n -= 0.5*Delta*sq(Delta)*(h[1]*gamma.x[1] - h[-1]*gamma.x[] +
			      h[0,1]*gamma.y[1] - h[0,-1]*gamma.y[]); 
   
    c[] = n/d;
  }
  
#if JACOBI
  foreach_level_or_leaf (l)
    h[] = (h[] + 2.*c[])/3.;
#endif
  
#if TRASH
  scalar a1[];
  foreach_level_or_leaf (l)
    a1[] = h[];
  trash ({h});
  foreach_level_or_leaf (l)
    h[] = a1[];
#endif
}

/**
$h$ and $p$ are solved in each timestep with the interface below. */

trace
mgstats layer (struct Layer p)
{
  if (p.dt == 0.) {
    mgstats s = {0};
    return s;
  }

  scalar f = p.f;
  scalar b[], fold[], fa[];
  const scalar lambda[] = - 1./p.dt;
  double m = p.m, n = p.n;
  face vector D[], beta[], gamma[];

  /**
  The first trial solution is that of the previous time step,
  $h^*=h(t- \delta t)$. */
  
  foreach() 
    fold[] = f[];
  //  boundary({fa});

  double maxres;
  mgstats s;
  int niter = 0, nitermax = 3;

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;
  
  /**
  In each iteration step of the N-R scheme, the auxiliary functions
  $\alpha(h^*)$, $\beta(h^*)$ and $\gamma_{x,y}(h^*)$ are first
  computed. */

  do {
    foreach_face()
      beta.x[] = face_gradient_x (f, 0);
    boundary_flux ({beta});
    
    scalar pres[];
    foreach() {
      pres[] = 0.;
      foreach_dimension()
	pres[] -= (beta.x[1]- beta.x[0])/Delta;
    }
    boundary({pres});
    
    face vector g[];
    foreach_face () {
      double ff = (f[]+f[-1])/2.;
      D.x[] = pow(ff, n);
      beta.x[] = -pow(ff, n-m-1.);
      gamma.x[] = n*pow(ff, n-1.)*(pres[]-pres[-1])/Delta
	-(n-m-1.)*pow(ff, n-m-2.)*(f[]-f[-1])/Delta;
      g.x[] = gamma.x[]*ff;
    }
    boundary_flux({g});
  
    foreach() {
      b[] = -fold[]/p.dt;
      fa[] = f[];
      foreach_dimension()
	b[] += (g.x[1] - g.x[])/Delta;
    }
    
    restriction ({beta, D, lambda, gamma});
    p.lambda = lambda;
    p.D = D;
    p.beta = beta;
    p.gamma = gamma;
  
    /**
    The MG linear solver is launched next, */

    s = mg_solve ({f}, {b}, residual_layer, relax_layer, &p);

    /**
    The convergence of the N-R is checked below. If the convergence is
    not reached and the number iteration is below a threshold value, a
    new iteration proceed.  */

    maxres = 0; 
    foreach (reduction(max:maxres)) {
      double res = fabs(f[]-fa[]);
      if (res > maxres)
	maxres = res;
    }
    boundary({fa});
    niter++;
    //    printf("maxres = %g, niter = %d\n", maxres, niter);
  } while(maxres > 1e-8 && niter < nitermax);

  /**
  We restore the default tolerance. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}

/**
## Time integration 

The variables we wish follow are: the height of the film *h*, the
inverse of the heigth *hi* and the level of refinement *l*. The problem
is fully symmetric and the the domain is a [-1.5:1.5]x[-1.5:1.5] square. */

#define T2OR 0 // active second order in time integration
#define pert1 0.05
#define pert2 0.03
#define HMIN A*(1-pert1)*(1-pert2)

scalar h[], hi[], l[];
double A = 0.05, n = 3., m = 1., DTMAX = 1e-3;
double dtold , hmin, hmin1, xmin = 0.5, ymin = 0.5;
int levelm;
int LEVEL= 12, lrejct = 0, LEVEL_MIN = 7;
double dtfactor = 0.995;

#if T2OR
scalar h2[], h1[], err[];
double errorstep = 0., DTMIN = 1e-6;
#endif

/**
initial timestep can be given by screen. */

int main(int argc, char * argv[]) {
  if (argc > 1)
    DTMAX = atof (argv[1]);

  periodic(right);
  periodic(top);
  origin(0.,0.);
  init_grid (1 << LEVEL_MIN);
  L0 = 3;
  run();
}

/**
Initially, a small sinusoidal pertubation is seeded. */

event init (i = 0) {
#if T2OR
  dtold = DTMIN;
  h2.nodump = true;
  h1.nodump = true;
  err.nodump = true;
#endif
  hmin = HMIN;
  hmin1 = HMIN + 1e-9;
  foreach() 
    h[] = A*(1.-pert1*cos(2*M_PI*(x/L0-0.5)))*
    (1.-pert2*cos(2*M_PI*(y/L0-0.5)));
  boundary ({h});
}

/**
Several timestep control procedures have been implemented. For example,
a time step control based in the error in time integration is computed
below. Has only sense if second order in time is used. */

#if T2OR
double ckstp (double dtmin, double dtmax, double dt, double error)
{
         
double dtgrow= 1.189207115;     
    
if ( error<1e-4 )
  {
    lrejct = lrejct + 1;
    if( (lrejct > 10) & (error < 2.5e-5))
      dt = dt * dtgrow;
  }
 else  // optional
   {
     dt=0.5*dt;
     lrejct = 0;
   }
 return ( min( dtmax, max(dtmin, dt))) ;
}
#endif


/**
Or the time step is proportional  to  the minimum heigth.  */

double mastep ( double hmin) {
  return 1.e-4*hmin;
}

/**
Or the time step depends on the maximum level of refinement. This criteria
although implemented has not been really used. */

double steplevel ( double levelmax) {
return 1.e-4/pow(2.,levelmax);
}


/**
In this case, the time step come after a linear interpolation. (Method
proposed by M.A. Herrada) */

double lininterp (double h1, double h2, double dt,
		  double alpha, double dtmax) {
  printf("%g %g %g %g %g \n",
	 h1,h2,dt,alpha, (alpha*h2-h1)/(h2-h1)*dt -dt);
  double step = min(((alpha*h2-h1)/(h2-h1)*dt -dt), dtmax);
    return (step <= 0. ? dtmax : step); 
}


event integration (i++) {
#if T2OR  
    dt = dtnext (ckstp (DTMIN, DTMAX, dtold, errorstep));
    dtold = dt;
    foreach() {
      h2[] = h[];
      h1[] = h[];
    }
  boundary({h2, h1});
  layer (h1, dt, m, n);  
  layer (h2, dt/2., m, n);  
  layer (h2, dt/2., m, n);

  double errormax = 0;
  foreach() {
    h[] = 2.*h2[] -h1[];
    err[] = fabs(h2[]-h1[]);
    if (err[] > errormax)
      errormax = err[];
  }
  boundary({h2, h1});
  stats s = statsf(err);
  errorstep = sqrt(s.sum);
#else
  dtold = dt;
  dt = dtnext (lininterp (hmin1, hmin, dtold, dtfactor, DTMAX));
  layer (h, dt, m, n);
#endif

  /** 
  Once *h[]* is computed, we calculate its inverse and the spatial
  distribution of the level of refinement. Beside, position and value
  of the minimum *hmin* is found. */
  
  double minres = 1000; 
  levelm = -1;
  hmin1 = hmin;
  foreach (reduction(min:minres)) {
    hi[] = 1/h[];
    l[] = level;
    double res = fabs(h[]);
      if(levelm < level)
	levelm = level;
    if (res < minres) {
      xmin = x;
      ymin = y;
      hmin = res;
      minres = res;
    }
  }
}

/**
In each timestep the grid is adapted accordingly to the error
estimation of the inverse of *h*. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({hi}, (double[]){0.2}, LEVEL, LEVEL_MIN);
}
#endif

/** 
# Outputs

![initial](hi_n3_m1_ho0.05_t0.png)
![initial](hi_n3_m1_ho0.05_t3.8944.png)

The plots of figure 12 of the paper are generated. */
  
event safe_figure (t = {0, 3.8944}) {
  view (fov = 20., tx = -0.5, ty = -0.5,
	width = 600, height = 600, samples = 1);
 squares("hi", linear =1, spread = -1);
 char name[80];
 sprintf(name,"hi_n%g_m%g_ho%g_t%g.png", n, m, A, t);
 save (name);
}

/** 
We could save also the time evolution of the minimum h. */  

#if 0
event output (i += 5) {
  static FILE * fp = fopen("hmin", "w");
  if (i==0)
    fprintf (fp, "#m = %g n = %g ho = %g LEVEL = %d LEVEL_MIN = %d dt = %g alpha = %g  \n",
	     m, n, A, LEVEL, LEVEL_MIN, DTMAX, dtfactor);
  fprintf (fp, "%.15g %.15g %.15g %d\n", t, hmin, dt, levelm);
  fflush(fp);
}
#endif
