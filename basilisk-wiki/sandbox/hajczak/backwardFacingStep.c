/** 
# Backward-facing step flow at *Re=800*.

This is the classical backward-facing step test case at *Re=800* in a bounded channel with aspect ratio of 2. 
Fluid is injected to the left, all the walls are no-slip. We use the centered Navier-Stokes solver, with embedded boundaries.

Reference data for this configuration are available in [Gartling (1990)](https://doi.org/10.1002/fld.1650110704) and [Erturk (2008)](https://doi.org/10.1016/j.compfluid.2007.09.003)

![Streamwise velocity contours](backwardFacingStep/Streamwise-velocity.png)

 */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "utils.h"

double Reynolds = 800.;
int maxlevel = 10;
double xmin = -5, xmax = 40.; // limits for region of interest
int lvl_coarse = 7; // coarsest level in channel
double tmax = 500.; // maximum simulation time 
double ev = 1e-5;   // tolerance on velocity to assess convergence
face vector muv[];

double h = 1.;
double U0 = 3./2.;
scalar un[];

/** 
This is a function for extracting aerodynamic quantities profiles at a given position in the channel. Probably not the most elegant way to do it in (Basilisk-)C.
*/

void extract_profile(double xp, vector u, vector gu, vector gv, double h)
{
  double uu[21], vv[21], dudx[21], dudy[21], dvdx[21], dvdy[21];
  char fname[20] ;
  snprintf(fname, 20, "profile_xh%g.dat",xp);
  FILE * fp = fopen (fname,"w");
  double yp = 1.;
  int j = 0;
  while (yp < 3.*h) {
    uu[j] =   interpolate (u.x ,xp,yp);
    vv[j] =   interpolate (u.y ,xp,yp);
    dudx[j] = interpolate (gu.x ,xp,yp);
    dudy[j] = interpolate (gu.y ,xp,yp);
    dvdx[j] = interpolate (gv.x ,xp,yp);
    dvdy[j] = interpolate (gv.y ,xp,yp);
    
    fprintf (fp,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t\n", (yp-1.)/2., uu[j], vv[j],
	     dvdx[j]-2*dudy[j], dudx[j], dudy[j]*2.,
	     dvdx[j], dvdy[j]*2.); /* y-profiles are rescaled between 0 and 1 for coherence with Erturk (beware of the typo in first column of tables) and Gartling */
    j++;
    yp = yp+2.*h/20.;
  }
  fflush (fp);
}

/**
You can pass arguments to the program in order to run a parametric study (*eg* on Reynolds number) in a shell
 */
int main(int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    Reynolds = atof (argv[2]);
  if (argc > 3)
    xmin = atof (argv[3]);
  if (argc > 4)
    xmax = atof (argv[4]);
  if (argc > 5)
    lvl_coarse = atoi (argv[5]);
  if (argc > 6)
    tmax = atof (argv[6]);
  if (argc > 7)
    ev = atof (argv[7]);

/**
The domain is 100 steps long, centered vertically. */
  
  L0 = 100. [1];
  origin (-20., -L0/2.);
  N = 32;
  mu = muv;
  
  run();
}

/**
We set a constant viscosity based on the Reynolds number, twice the
  step height $2h$ (hydraulic radius) and the inlet velocity
  mean $2/3 \times U_0 = 1$. */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*(2.*h*(2./3.)*U0)/Reynolds;
}


/**
We impose a Poiseuille flow with a centerline velocity of 3/2 to get a unit mean velocity 
flow inside the upstream channel. An outflow
condition is used on the right boundary. */

u.n[left]  = (fabs(y-2.5*h) < h/2.) ? dirichlet(4*U0*((y-2*h)/h)*(1-(y-2*h)/h)) : dirichlet(0.);
u.t[left]  = dirichlet(0.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The channel walls are no-slip. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

		       
event init (t = 0)
{
  refine (fabs(y-2*h) <= 1.5*h && level < maxlevel); // get a nice resolution in the channel
  unrefine(x>xmax && level > lvl_coarse); // coarsen elsewhere
  unrefine(x<xmin && level > lvl_coarse);

  /** Seems that [embedded boundaries need to not coincide with the grid](https://groups.google.com/g/basilisk-fr/c/SLnKb8U0BXs). I have struggled for a while with 
      this issue */
  
    solid (cs, fs, union(intersection(-x, intersection (3. - y, y-2)),
			 intersection(x, intersection (3. - y, y-1.))));
    
    /** Somehow seems to work. The channel is $2<y<3$ for $x<0$ and $1<y<3$ for $x>0$.*/
    
    /**
       We set the initial velocity field as Poiseuille in both upstream and downstream parts of the 
       channel.*/

    foreach()
      {
	u.x[] = (cs[] && x < 0) ?
	  4*U0*((y-2*h)/h)*(1-(y-2*h)/h) : (cs[] && x > 0) ?
	  4*U0*((y-h)/(2*h))*(1-(y-h)/(2*h)) : 0.;
	un[] = u.y[]; // allows to track convergence
      }
}



/**
We check the number of iterations of the Poisson and viscous
problems and compute the recirculation bubbles positions. The info are written in the log file */

event logfile (i++, t <= tmax)
{
  double du = change (u.y, un);
  if (i > 1 && du < ev)
    return 1; /* stop */

  /** Finding the bubbles positions. We scan $y=1.1$ and $y=2.9$ lines which is not very nice. Possible fix using   [embed_vorticity](http://basilisk.fr/src/embed.h#surface-force-and-vorticity) ?
  */
  
  double X1 = 0.5;
  while (interpolate (u.x ,X1, 1.1) < 0 && X1 < xmax)
    X1=X1+0.1;		       
  			       
  double X4 = X1;	       
  while (interpolate (u.x ,X4, 1.1) > 0 && X4 < xmax)
    X4=X4+0.1;		       
  			       
  double X2 = 0.;	       
  while (interpolate (u.x ,X2, 2.9) > 0 && X2 < xmax)
    X2=X2+0.1;		       
  			       
  double X3 = X2;	       
  while (interpolate (u.x ,X3, 2.9) < 0 && X3 < xmax)
    X3=X3+0.1;
    
  fprintf (stderr, "%d %g %d %d %g %g %g %g %g\n", i, t, mgp.i, mgu.i, du, X1, X2, X3, X4);
  
}

/** Outputs for comparison with reference data */
event profile (t=end)
{
  scalar omega[], l[], m[];
  vector gu[], gv[];

  foreach()
    {
      l[] = level;
      m[] = cs[]-0.1;
    }

  gradients ({u.x,u.y}, {gu,gv});  

  output_ppm (u.x, file="Streamwise-velocity.png",
	      n=512, box={{-5.,1.},{20.,3.}},
	      min=-1, max=1, linear=true, mask=m);


  extract_profile(6.,  u, gu, gv, h);
  extract_profile(14., u, gu, gv, h);
  extract_profile(30., u, gu, gv, h);

  dump();

}


/**
Script for plotting data in python
~~~pythonplot Profiles of aerodynamic quantities
import numpy as np
import matplotlib.pyplot as plt

varnames=['u','v','omega','dudx','dudy','dvdx','dvdy']
xlabels=['$u$','$v$','$\omega$', '$\partial u / \partial x$', '$\partial u / \partial y$', '$\partial v / \partial x$', '$\partial v / \partial y$']
pos=[6,14,30];
for p in pos:
    ref_e = np.loadtxt('dataErturk/flowVariables_xh'+str(p)+'_Re800.dat')
    ref_g = np.loadtxt('dataGartling/flowVariables_xh'+str(p)+'_Re800.dat')

    res1 = np.loadtxt('profile_xh'+str(p)+'.dat')
    res2 = np.loadtxt('maxlevel11_Re800/profile_xh'+str(p)+'.dat')
    res3 = np.loadtxt('maxlevel12_Re800/profile_xh'+str(p)+'.dat')

    f, axes = plt.subplots(3,3, sharex=False, sharey=False, figsize=(10,10))
    ax= np.ravel(axes)
    for i in range(1,8):
        ax[i-1].plot( ref_e[:,i],ref_e[:,0] , 'ok', label='Erturk (2008)')
        if (p != 6):
            ax[i-1].plot( ref_g[:,i],ref_g[:,0] , 'sr', label='Gartling (1990)')
        ax[i-1].plot(res1[:,i],res1[:,0], label='Basilisk (Max. level = 10)')
        ax[i-1].plot(res2[:,i],res2[:,0], label='Basilisk (Max. level = 11)')
        ax[i-1].plot(res3[:,i],res3[:,0], label='Basilisk (Max. level = 12)')
        ax[i-1].set_xlabel(xlabels[i-1])
        ax[i-1].set_ylabel('$y/2h$')
    f.delaxes(ax[7])
    f.delaxes(ax[8])
    #ax[6].legend(bbox_to_anchor=(0.5, -0.1))
        #ax[i].title('Flow variables at $x/h = $'+ str(p)+', $Re = 800$')
    plt.tight_layout()
    plt.savefig('profiles_xh'+str(p)+'.png')
    plt.close()

~~~
![Comparison with data at $x/h=6$. $\textcolor{red}{\blacksquare}$ : Gartling (1990), $\bullet$:
Erturk (2008), ($\textbf{\textcolor{blue}{-}}$): Basilisk (max level=10), 
($\textbf{\textcolor{orange}{-}}$): Basilisk (max level=11), ($\textbf{\textcolor{green}{-}}$): Basilisk (max level=12)](backwardFacingStep/profiles_xh6.png)

![Comparison with data at $x/h=14$ $\textcolor{red}{\blacksquare}$ : Gartling (1990), $\bullet$:
Erturk (2008), ($\textbf{\textcolor{blue}{-}}$): Basilisk (max level=10), 
($\textbf{\textcolor{orange}{-}}$): Basilisk (max level=11), ($\textbf{\textcolor{green}{-}}$): Basilisk (max level=12)](backwardFacingStep/profiles_xh14.png)

![Comparison with data at $x/h=30$ $\textcolor{red}{\blacksquare}$ : Gartling (1990), $\bullet$:
Erturk (2008), ($\textbf{\textcolor{blue}{-}}$): Basilisk (max level=10), 
($\textbf{\textcolor{orange}{-}}$): Basilisk (max level=11), ($\textbf{\textcolor{green}{-}}$): Basilisk (max level=12)](backwardFacingStep/profiles_xh30.png)
 */


/**
* NB : There seems to be a typo in Erturk -> $y$-gradients are wrong if the first column is $0<y<2$ (should be rescaled to $0<y<1$ for consistency with Gartling who considers the downstream channel to be of height $h$ and therefore a step of height $h/2$).
* $x$-gradients are not very well captured
*/

/**
# Recirculation bubbles position with Reynolds

The code has been run with increasing Reynolds for the same target residual (not reached for the highest $Re$, $\sim 10^{-1}$ at most) with the same grid settings.

~~~pythonplot 
import numpy as np
import matplotlib.pyplot as plt

Reynolds=[100*i for i in range(1,16)]

X=np.loadtxt('positionsFlowReattachmentComputed.dat')
b = np.loadtxt('dataErturk/positionsFlowReattachment.dat')

plt.plot(Reynolds,X[:,0], 'bo', label='X1 (Basilisk)')
plt.plot(b[:,0],b[:,3],'ro'   , label='X1 (Erturk, 2008)'  )

plt.plot(Reynolds,X[:,1], 'bs', label='X2 (Basilisk)')
plt.plot(b[:,0],b[:,4],'rs'   , label='X2 (Erturk, 2008)'  )

plt.plot(Reynolds,X[:,2], 'bd', label='X3 (Basilisk)')
plt.plot(b[:,0],b[:,5],'rd'   , label='X3 (Erturk, 2008)'  )

plt.xlim(50.,1550.)
plt.xlabel('$Re$')
plt.ylim(0,35.)
plt.ylabel('$X_n$')

plt.title('Recirculation bubbles positions')
plt.legend()
plt.savefig('Xn.png')

~~~
*/

/**
## References
* Erturk, E. (2008). Numerical solutions of 2-D steady incompressible flow over a backward-facing step, Part I: High Reynolds number solutions. Computers & Fluids, 37(6), 633-655.
* Gartling, D. K. (1990). A test problem for outflow boundary conditions—flow over a backward‐facing step. International Journal for Numerical Methods in Fluids, 11(7), 953-967.
*/
	
