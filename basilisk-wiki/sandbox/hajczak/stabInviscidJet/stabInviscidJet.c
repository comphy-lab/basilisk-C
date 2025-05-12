/** 
## Sinuous instability of $Re=2 \times 10^3$ Bickley jet
In flute-like instruments, the source of sound results from the interaction between an air jet and an edge. The so-called Jet-drive model is based on a description of the jet instability waves that develop from the channel outflow down to the edge. 

The amplification and phase speed of the disturbances are user-prescribed and correspond to the results that can be obtained by studying the linear stability of an inviscid and incompressible Bickley jet using the Rayleigh equation.

Theoretical results on this case are to be found in Drazin & Howard for instance (see p.40)

Experimental results have been published by Mattingly & Criminale, and later by Nolle and show a very nice agreement with theory. Nolle, on the basis of earlier claims (end of page 3693), rejected the possibility of a yet theory-compatible decaying wave near the jet exit and considered an exponentially growing only disturbance. It is important to note that the experimental nature of Nolle's work makes it difficult to measure this evanescent mode. Furthermore, Nolle is bound to measuring the jet centerline position and not the velocity. In our numerical study, using the jet centerline instead of the vertical velocity on the centerline led to important precision issues due to the necessity to interpolate the velocity profile in the first case.

Some numerical results have been obtained by Nagy & Paal who have also investigated viscous effects by using the more general Orr-Sommerfeld equation. Nagy et al also used CFD to obtain the jet-drive parameters in another study on tenor recorder geometry. Nagy & Paal report the apparition of a so-called "strange mode" at the jet exit.

We show that :

* There is no need for additional terms in the conservation equations to simulate this situation, nor is there for compressibility to be considered

* Basilisk simulations can reproduce theoretical results to a certain degree of accuracy. (Comparison with other authors ????)

* The "strange" mode reported by Nagy & Paal is nothing but the evanescent mode that Nolle rejected. We show that the amplification coefficient and phase speed of this numerical mode is consistent with theory.

We use the centered Navier-Stokes solver and set the Reynolds number to a value of 2000 for comparison with the inviscid theory. We record the vertical velocity perturbations $v_y$ on the jet centerline $y=0$. The jet is excited at the domain left boundary at a pulsation $\omega_r$. After a transient time, we record the perturbations for a total duration of 3 excitation periods. The jet has a squared hyperbolic secant profile with unit velocity on the centerline. The vertical coordinate is made non-dimensional with a length $b$ such that the jet width is $2b$. The amplitude of the imposed disturbance is 1% of the jet centerline velocity so that the linear stability analysis framework is respected.


 */

#include "navier-stokes/centered.h"
#include "utils.h"

double Reynolds = 2e3;
int maxlevel = 9;
int Nper = 3;  
double omega_r = 0.25 [0,-1];
double b = 1. [1,0];
double U0 = 1. [1,-1];
double xmin =0.;
double xmax=10.;
double transient_time=10.;
double  u_ac = 1e-2 [1,-1];

face vector muv[];

/**
A function for extracting (by interpolation) the vertical velocity on the $y=0$ line.

Note that this is very different from Nolle's experimental approach which cannot have 
easy access to this quantity and therefore relies on measuring the streamwise velocity 
profile, interpolating it to finally deduce the jet centerline at any $x$ position. 
Indeed, working on $v_y$ directly saves us the sech2 interpolation and thus some error.

*/

void extractPerturb(vector u)
{
  FILE *fp;
  double xp = xmin;
  char fname[50] ;
  snprintf(fname, 50, "omega_%.2f_vy.dat",omega_r); 
  fp = fopen (fname,"a");
  while (xp < xmax) {
    fprintf (fp,"%g\t", interpolate (u.y ,xp, 0.));
    xp=xp+0.2;
  }
  fprintf (fp,"\n");
  fclose(fp);
}

/**
The sech and sech2 functions are defined for clarity in the BC
 */

double sech2(double y)
{
  return 1./(sq(cosh(y)));
}
double sech(double y)
{
  return 1./cosh(y);
}


int main(int argc, char * argv[])
{

/**
The domain is 100 jet half-widths long. This is to avoid the vorticity generated when the jet breaks to interact back with the jet as large recirculation can appear in the simulation domain. The domain is also centered vertically.

TODO : study the dependancy of results if a (small) streamwise velocity is imposed to "wash" out the eddies.

The excitation pulsations are taken as those considered by Nolle for comparison
*/


  int i;
  double omegaList[13] = {0.04,0.05,0.07,0.10,0.14,0.20,0.25,0.40,0.50,0.60,0.80,1.00,1.20};
  for (i = 0; i < 13; i++)
    {
      L0 = 100. [1]; 
      origin (0., -L0/2.);
      N = 128;
      mu = muv;
      omega_r = omegaList[i];
      run();
    }
}

/**
The viscosity is everywhere imposed in accordance to the Reynolds number

TODO : Set the viscosity to zero and see what happens to the results
*/

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*(b*U0)/Reynolds; 
}


/**
Boundary conditions : we impose a sech2 profile plus a perturbation, of which analytical expression is taken from Drazin & Howard. It is important to note that the velocity profile must be imposed on the full $y$ range so that the jet is as receptive as in the theory.

TODO : study the varicose mode. So far I have only managed to get the sinuous mode (or do I have sinuous + varicose altogether ? This would explain some funny - ie not double-exponential like - figures at certain frequencies)

*/

u.n[left]  = dirichlet(U0*sech2(y/b) + u_ac*cos(omega_r*t)*(-2*tanh(y/b)*sech2(y/b))); //SINUOUS
u.t[left]  = dirichlet(u_ac*sech2(y/b)*cos(omega_r*t)) ; // SINUOUS

//u.n[left]  = dirichlet(U0*sech2(y/h) + u_ac*cos(omega_r*t)*sech(y/h)*(2*sech2(y/h)-1)); //VARICOSE
//u.t[left]  = dirichlet(u_ac*sech(y/h)*tanh(y/h)*cos(omega_r*t)) ; // VARICOSE 


p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
Initialise the velocity fied with the sech2 profile
*/

event init (t = 0)
{
  foreach()
    u.x[] =  U0*sech2(y/b);
}

/**
Use wavelet mesh adapation with a tolerance of 1e-3 on velocity field. How does this tolerance affect the results ?

It is not very satisfying to resolve the whole domain, which is long, with the same accuracy, while we are only interested in the linear region. Can we use some sort of region-specific adaptation ? This would look like https://groups.google.com/g/basilisk-fr/c/EMg6USbSVq0/m/R79XWddDBAAJ 
*/ 
 
 
event adapt (i++) { 
  adapt_wavelet ({u}, (double[]){1e-3,1e-3}, maxlevel, 4); 
}


/** 
Check the progress of the simulation
 */

event logfile (i++)
{    
  fprintf (stderr, "%d %g %d %d\n", i, t/(transient_time+(Nper)*2*pi/omega_r), mgp.i, mgu.i);  
}

/** Outputs for comparison with reference data */

event profile (t=transient_time; t+=(1./128.)*(2*pi/omega_r) ; t<transient_time+(Nper*2*pi/omega_r)) // With unit velocity on the centerline, t=10 means the center has travelled to x=10
{
  extractPerturb(u); 
}

/**
Script for plotting data in python
~~~pythonplot aaa

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import signal

def dble_expo(x,eta0,alpha0,eta1,alpha1): # We look for a linear sum of amplified and evanescent modes. The constraint on alpha0,1 signs is applied in the curve_fit procedure
    return eta0*np.exp(alpha0*x) + eta1*np.exp(alpha1*x) 

def aff(x,a,b):
    return a*x+b

def bestfit(x,dph): # Look for the "best" linear fit starting at the beginning/end of the x interval
    res1 = 1e3
    res2 = 1e3

    for i in range(3,len(x)-2):
        popt = np.zeros(3)
        pcov=np.ones((3,3))
        try:
            popt, pcov = curve_fit(aff,x[-i:],np.abs(dph[-i:]))
        except RuntimeError:
            pass
        if  np.sqrt(np.diag(pcov))[1] < res1:
            res1 = np.sqrt(np.diag(pcov))[1]
            i1=i
        try:
            popt, pcov = curve_fit(aff,x[:i],np.abs(dph[:i]))
        except RuntimeError:
            pass
        if  np.sqrt(np.diag(pcov))[1] < res2:
            res2 = np.sqrt(np.diag(pcov))[1]
            i2=i
    return i1,i2

Ns=128
dx=0.2
x1=0.
#x2 is a sensitive parameter. The jet breaks at lower x when the amplification is higher
Nper=3

results=np.ones((13,5))
k=0
for omega in [0.04,0.05,0.07,0.10,0.14,0.20,0.25,0.40,0.50,0.60,0.80,1.00,1.20]:
    if omega < 0.2:
        x2=10.
    else:
        x2=7.
    Nx=int(np.round((x2-x1)/dx + 1))
    x = np.linspace(x1,x2,Nx)
    dph=np.zeros(Nx)
    alpha=np.zeros(Nx)
    vy=np.transpose(np.loadtxt('omega_'+format(omega,".2f")+'_vy.dat'))
    for i in range(Nx):
        freq = np.fft.fftfreq(len(vy[i,:]),d=1/Ns)
        sp = np.fft.fft(vy[i,:])
        alpha[i]=(2*(1./(Nper*Ns))*np.abs(sp[Nper])) # Amplitude taken as Fourier coeff at excitation freq
        phi0 = np.angle(np.fft.fft(vy[0,:]))[Nper] 
        phi = np.angle(np.fft.fft(vy[i,:]))[Nper] 
        dph[i] = phi-phi0
    dph=np.unwrap(dph) # Allows having the phase in linear form
        
    i1,i2 = bestfit(x,np.abs(dph))
    popt_cp1, pcov_cp1 = curve_fit(aff,x[-i1:],np.abs(dph)[-i1:])
    popt_cp2, pcov_cp2 = curve_fit(aff,x[:i2],np.abs(dph)[:i2])        
    popt,pcov = curve_fit(dble_expo,x,alpha, bounds=[(0.,0.,0.,-2.),(1e-2,1.,1e-2,0.)])

    if (omega==0.20):
        plt.subplot(121)
        plt.plot(x,np.abs(dph),'b-')
        plt.plot(x[int(Nx/2.):],aff(x[int(Nx/2):],*popt_cp1),'k-',linestyle='dashed')
        plt.plot(x[:int(Nx/2)],aff(x[:int(Nx/2)],*popt_cp2),'k-',linestyle='dashed')
        plt.title('$\omega_r=0.2$')
        plt.ylabel('Phase difference')
        plt.xlabel('$x/b$')
        
        plt.subplot(122)
        plt.plot(x,alpha,'b-')
        plt.plot(x,dble_expo(x,*popt),'k-',linestyle='dashed')
        plt.title('$\omega_r=0.2$')
        plt.ylabel('Amplitude')
        plt.xlabel('$x/b$')
        plt.savefig('omega02.png')
    results[k,:]=(omega, popt[3], popt[1], omega/popt_cp1[0], omega/popt_cp2[0])
    k=k+1

# Figures generation
eig_sp = np.loadtxt('../eigenvalues_spatialTheory.dat',skiprows=1)

Nolle_cp_05= np.loadtxt('../expNolle_cp_flueWidth0.5mm.dat')
Nolle_cp_1= np.loadtxt('../expNolle_cp_flueWidth1mm.dat')
Nolle_mub_05= np.loadtxt('../expNolle_mub_flueWidth0.5mm.dat')
Nolle_mub_1= np.loadtxt('../expNolle_mub_flueWidth1mm.dat')

plt.subplot(221)
plt.semilogx(eig_sp[:,0],-eig_sp[:,2], color='k', label='Theory')

plt.semilogx(Nolle_mub_05[:,0],Nolle_mub_05[:,1],marker='o',linewidth=0, color='r', label='Nolle (1998) - Exp.')
plt.semilogx(Nolle_mub_1[:,0],Nolle_mub_1[:,1],marker='o',linewidth=0  , color='r')

plt.semilogx(results[:,0],results[:,2],marker='s',linewidth=0., color='b', label='Basilisk - Num.')

plt.legend()
plt.xlabel('$\omega_r$')
plt.ylabel('$-α_i$')
plt.title('Amplified mode')

plt.subplot(223)
plt.semilogx(eig_sp[:,0],eig_sp[:,0]/eig_sp[:,1], color='k', label='Theory')

plt.semilogx(Nolle_cp_05[:,0],Nolle_cp_05[:,1],marker='o',linewidth=0, color='r', label='Nolle (1998) - Exp.')
plt.semilogx(Nolle_cp_1[:,0],Nolle_cp_1[:,1],marker='o',linewidth=0  , color='r')

plt.semilogx(results[:,0],results[:,3],marker='s',linewidth=0., color='b', label='Basilisk - Num.')

plt.legend()
plt.ylabel('$c_p$')
plt.xlabel('$\omega_r$')

plt.subplot(222)
plt.semilogx(eig_sp[:,0],eig_sp[:,1], color='k', label='Theory')

plt.semilogx(results[:,0],-results[:,1],marker='s',linewidth=0., color='b', label='Basilisk - Num.')

plt.legend()
plt.xlabel('$\omega_r$')
plt.ylabel('$α_r$')
plt.title('Evanescent mode')

plt.subplot(224)
plt.semilogx(eig_sp[:,0],-eig_sp[:,0]/eig_sp[:,2], color='k', label='Theory')

plt.semilogx(results[:,0],results[:,4],marker='s',linewidth=0., color='b', label='Basilisk - Num.')
plt.ylim(0,10)

plt.legend()
plt.ylabel('$c_p$')
plt.xlabel('$\omega_r$')
plt.savefig('results.png')

~~~
*/

/**
![Streamwise evolution of disturbance phase (left) and amplitude (right). Dotted lines : fitted piecewise linear and double exponential models](stabInviscidJet/omega02.png)


 * TODO : perform $k-\omega$ transform

 * NB : poor predictions at lowest frequencies and highest frequencies for the evanescent mode
*/

/**
## References
To be completed
*/
	
