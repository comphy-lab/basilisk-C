# -*- coding: utf-8 -*-
"""

Stability of the plane Poiseuille flow in primitive variables

Created on Dec 3 2025

@author: D. Fabre (adapted from a previous Matlab program)


 *This document belongs to the [easyPYstab](/sandbox/easystab/easyPYstab/README.md) project, please consult the main page of the project for explanations.
 
 *The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapter 7](/sandbox/easystab/M2DET/Instabilities.md#lecture-7-shear-flow-instabilities-ii-viscous-instabilities)*


 In this code, we solve the linear eigenvalue problem for stability analysis of the Plane Poiseuille flow. 
 This is the flow in a flat channel between two infinite plates, driven by a pressure gradient. 

  We perform the analysis using the primitive variables u,v and p the horizontal and vertical components
  of the velocity and the pressure. Thus we use the Navier-Stokes equations, linearized 
  about a parabolic velocity profile, plus as well the continuity equation. 
  Since there is no time derivative in the continuity equation for incompressible flow this 
  system becomes a *descriptor system*, also known as *differential algebraic equation*. 
  This is a nice formulation because there are no fourth order derivative as for the Orr-Sommerfeld equation.
 
 The parameters are the wave number $k=2\pi/\lambda$ where $\lambda$ is the wavelength in 
 the direction of the flow, the height of the domain $L$ which must be $2$, the Reynolds number 
 and the number of grid points.
"""

import EasyStab as ES                     # coeur de la toolbox easystab

import matplotlib.pyplot as plt           # Pour tracer
import numpy as np                        # Pour calculer
from numpy import*                        # Pour construire les tableaux
from matplotlib.pyplot import*            # Pour construire le graphe
from scipy.linalg import eig              # Pour resoudre pb generalise
from scipy.sparse.linalg import eigs      # Pour resoudre pb generalise avec methode iterative
from scipy.optimize import root_scalar

plt.close('all')

"""
# Definition of functions for this program
"""

"""
## DEFINITION OF FUNCTION solve_eigenproblem
"""
def solve_eigenproblem(Re,k):

# %{
# # Construction of the matrices

# The problem is written as follows:

# $$
# \lambda  {\mathcal B} \, \hat{q} = {\mathcal A} \, \hat{q}
# $$

# with 
# $$
# {\mathcal B} = 
# \left[
# \begin{array}{ccc} 
# 1 & 0 & 0 \\ 
# 0 & 1 & 0 \\
# 0 & 0 & 0  
# \end{array} 
# \right] 
# $$

# $$
# {\mathcal A} = 
# \left[
# \begin{array}{ccc} 
# -i k \bar{U} + Re^{-1} ( \partial_y^2 - k^2) & - \partial_y \bar{U} & - i k \\ 
# 0 & -i k \bar{U} + Re^{-1} ( \partial_y^2 - k^2) & - \partial_y \\
# i k  & \partial_y  & 0 
# \end{array} 
# \right] 
# $$

# Here we build the dynamics matrices $A$ and $B$. 
# %}

  Z=np.zeros((Nx+1,Nx+1)); 
  I=np.eye(Nx+1); 
  dx=1j*k*I; dxx=-k**2*I;
  Delta=dxx+dyy;
  S=-diag(U)*dx+Delta/Re;

  A=np.block( 
            [ [S,   -diag(Uy), -dx] ,
              [Z,    S,        -dy] ,
              [dx,   dy,        Z ] ]
           )

 
  B=np.block( 
            [ [I,Z,Z] ,
              [Z,I,Z] ,
              [Z,Z,Z] ]
           );


  #%{
  # Boundary conditions
  #%}
  #% We have to modify the lines corresponding to first and last point for u and v.
  #% We first determine the corresponding indices. Here are the numbers of the four lines to modify in the matrices   
  loc=(0, n-1, n, 2*n-1 );  

  # Now we rewrite lines number "loc" in A with 1 on the diagonal, 0 elsewwhere
  III=np.eye(3*n);C=III[loc,:];A[loc,:]=C; 

  # Finally rewrite lines number "loc" in B with zeros 
  B[loc,:]=0; 


  #%{
  # Computing the eigenmodes
  #We can use the *eig* function for DAE just as for dynamic systems. the difference is that we will get infinite eigenvalues. Thinking of the incompressible system as the limit of the compressible one, you realize that these infinite eigenvalues are the limit of the sound modes. We remove them, we sort the eigenvalues and eigenvectors.
  #%}

  #% computing eigenmodes 
  s,em = eigs(A,M=B,k=10,sigma=1+1j);
  
  #% sort and remove spurious eigenvalues
  S = s[abs(s)<1e8];
  EM = em[:,abs(s)<1e8]; 
  order = np.argsort(real(-S))
  S = S[order];
  EM = EM[:,order];
  return S,EM


def leading_eigenvalue(Re,k):
   S,EM = solve_eigenproblem(Re,k)
   return S[0]
   
"""
## Definition of two functions for postprocessing (plot mode structure)
"""

def plotmode(event,mode,S,k):
    # This function is used to represent the structure of an eigenmode
    global y,dy,dyy
    n = y.shape[0];
    u = mode[0:n];
    v = mode[n:2*n];
    p = mode[2*n:3*n];
    MM =max(abs(mode))
    vorticity = (dy@u)-1j*k*v;
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 2]})
    ax1.plot(real(u),y,'b-',imag(u),y,'b--');
    ax1.plot(real(v),y,'g-',imag(v),y,'g--');
    ax1.plot(real(p),y,'k-',imag(p),y,'k--');
    ax1.set_ylabel('y');ax1.set_xlim(-MM*1.2,MM*1.2); 
    ax1.set_ylim(-1,1); 
    ax1.legend({'Re(u)','Im(u)','Re(v)','Im(v)','Re(p)','Im(p)'},loc='lower center')
    if (-imag(S)/k)<1 and (-imag(S)/k)>0:
      ycrit = sqrt(1+imag(S)/k);
      ax1.plot([-MM,MM],[ycrit,ycrit],'k:',[-MM,MM],[-ycrit,-ycrit],'k:')    
    Lx=2*pi/k; Nx =30;  x=linspace(-Lx/2,Lx/2,Nx);
    yy = y;
    pp = 2*real(np.outer(p,exp(1j*k*x)));
    uu=2*real(np.outer(u,exp(1j*k*x)));
    vv=2*real(np.outer(v,exp(1j*k*x)));
    vorticityvorticity=2*real(np.outer(vorticity,exp(1j*k*x)))
    ax2.contourf(x,yy,vorticityvorticity,10); 
    #ax2.quiver(x,yy,uu,vv,'k'); 
    ax2.set_xlabel('x'); 
    z_str = f"{S.real:.3f} + {S.imag:.3f}j"
    plt.title("Structure of the eigenmode (vorticity)\n" +" lambda = "+z_str);


def on_pick(event):
    # This function is defined to allow plotting of eigenmodes by clicking on corresponding eigenvalueq
    global EM
    ind = event.ind[0]  # Index du point cliqué
    EMM =EM[:,ind]
    z_str = f"{S[ind].real:.3f} + {S[ind].imag:.3f}j"
    print("Clic sur valeur propre "+str(ind)+" : Eigenvalue = "+z_str+ " ind " +str(ind))
    plotmode(event, EMM, S[ind], k)
    plt.show()


"""
# Main program
"""

"""
## Discretization and base flow
"""

print("\n #### Selection of parameters for this program ####  \n")

Nx = ES.myinput("==> Nx, number of points ? ",60); # The order of the Chebychev expansion
n = Nx+1;    # The number of collocation points (pseudo-grid-points) 
DiscetizationType = ES.myinput("==> Discretization type ? ","Cheb")

# %{
# # Differentiation matrices
# We need to compute the derivatives in $x$ and also in $y$. 
# But in fact as is usually done for stability of parallel flows, 
# we can do a Fourier transform in the direction where the system does not change, 
# so here the numerical differentiation is done in $y$, and the differentiation in $x$ 
# simply ammounts to multiplication with $ik$.
# %}

print("   (building differentiation matrices and constructing base flow)")

# Differentiation matrices
y,dy,dyy,w = ES.dif1D(DiscetizationType,parameters={"a":-1,"b":1},N=Nx); 

# Definition of the base flow : PLANE POISEUILLE FLOW. We also need its y-derivative.
U=1-y**2; Uy=-2*y;



"""
## Compute a spectrum and plot eigenvalues for one set of (k,Re)
"""

rep = ES.myinput("\n   Do you want to compute a spectrum and plot eigenmodes ?  ","yes")
if rep=="yes":
    
    # Selecting the parameters
    k  = float(ES.myinput("==> Wavenumber k ?    ",1));    # the wave number
    Re = float(ES.myinput("==> Reynolds number ? ",10000));    # the Reynolds number
    print("   (computing spectrum)")
    S,EM = solve_eigenproblem(Re,k)
    
    # We plot the spectrum
    print("   (plotting spectrum in figure 1)")
    fig, ax = plt.subplots()
    h, = plt.plot(-np.imag(S), np.real(S), 'r+',picker=5)
    plt.xlabel('\omega_r');plt.ylabel('\omega_i');plt.ylim([-.5,.05])
    fig.canvas.mpl_connect('pick_event', on_pick) # to associate the previously defined function to click events
    ax.set_title('Spectrum for Re = '+str(Re)+' ; k = '+str(k)+'\n (Click on eigenvalues to display corresponding modes)')
    plt.show()
    
    # We plot the leading eigenmode
    print("   (plotting leading eigenmode in figure 1)")
    plotmode(1,EM[:,0],S[0],k)
    plt.show()



"""
## Draw curves omega_i(k) and omega_r(k) for several values of Re
"""

rep = ES.myinput("\n   Do you want to draw curves omega(k) for several values of Re ?","yes")
if rep=="yes":
    alphatab = np.linspace(.1,1.5,101);
    fig, (ax1, ax2) = plt.subplots(2, 1)
    for Re in [5772, 6000, 8000, 10000, 20000, 50000, 100000]:
       print("   ( Computing branch for Re = "+str(Re)+" )");
       stab = zeros(len(alphatab),dtype=complex);
       for j in range(len(alphatab)):
           k= alphatab[j]   
           S,EM = solve_eigenproblem(Re,k)
           stab[j] = S[0]
       # plot the curves        
       ax1.plot(alphatab,real(stab),label="Re = "+str(Re));
       ax2.plot(alphatab,-imag(stab)/alphatab,label="Re = "+str(Re)); 
       plt.show()
    fig.suptitle("amplification rate and phase velocity vs. k")
    ax1.set_xlabel("k");ax1.set_ylabel(r"$\omega_i$");
    ax1.legend();ax1.grid();ax1.set_ylim([-.01,.01])
    ax2.set_xlabel("k");ax2.set_ylabel(r"$\omega_r/k$");
    ax2.legend();ax2.grid();ax2.set_ylim([.05,.3])



"""
## Finally : reconstuct the marginal stability curve in the (Re,k) plane
"""

rep = ES.myinput("\n   Do you want to Reconstuct the marginal stability curve in the (Re,k) plane ? ","yes")
if rep=="yes":
    
    Retab = np.array([5772, 5800, 6000, 6500, 7000, 8000,  10000, 15000, 20000, 35000, 50000, 75000, 100000])
    # point de depart : 57772 / 1.02
    kinftab = [1.02]
    ksuptab = [1.02]
    # 2e point : guess sur Retab[1]
    Re = Retab[1]
    sol = root_scalar(lambda a: np.real(leading_eigenvalue(Re, a)),
                      bracket = [1,1.01]);  
    kinftab.append(sol.root)
    sol2 = root_scalar(lambda a: np.real(leading_eigenvalue(Re, a)),
                      bracket = [1.03,1.04]); 
    ksuptab.append(sol2.root)    
    # boucle pour les suivants
    for Re in Retab[2:]:
        print("calcul points sup et inf pour Re =", Re)
        sol = root_scalar(lambda a: np.real(leading_eigsenvalue(Re, a)),
                          x0=kinftab[-1]*.9, method='newton')
        kinftab.append(sol.root)
        sol2 = root_scalar(lambda a: np.real(leading_eigenvalue(Re, a)),
                           x0=ksuptab[-1], method='newton')
        ksuptab.append(sol2.root)
        print("   ... "+str(sol)+" ; "+str(sol2))
    plt.figure()
    plt.plot(Retab, kinftab, 'b+-', label='kinf')
    plt.plot(Retab, ksuptab, 'bx-', label='ksup')
    plt.xlabel('Re')
    plt.ylabel('k')
    plt.title('Marginal stability curve of the plane Poiseuille flow')
    plt.legend()
    plt.show()

