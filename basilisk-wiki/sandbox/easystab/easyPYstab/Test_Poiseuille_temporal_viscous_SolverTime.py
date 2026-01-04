"""

**Stability of 2D poiseuille flow : test of time needed by solving eigenvalue problem**

 *This document belongs to the [easystab](/sandbox/easystab/README) project, please consult the main page of the project for explanations.*

 *The present program was specifically designed as a support for the course ["Introductions to hydrodynamical instabilities"](/sandbox/easystab/M2DET/Instabilities.md) for the "DET" Master's cursus in Toulouse. For the underlying theory please see [Lecture notes for chapter 7](/sandbox/easystab/M2DET/Instabilities.md#lecture-7-shear-flow-instabilities-ii-viscous-instabilities)*

"""

import matplotlib.pyplot as plt           # Pour tracer
import numpy as np                        # Pour calculer
from numpy import*                        # Pour construire les tableaux
from matplotlib.pyplot import*            # Pour construire le graphe
from scipy.linalg import eig              # Pour resoudre pb generalise
from scipy.sparse.linalg import eigs      # Pour resoudre pb generalise avec methode iterative
from EasyStab_TestPython import dif1D
import time


#
# Numerical Parameters
#


Nx=100;      # The order of the Chebychev expansion
n = Nx+1;    # The number of collocation points (pseudo-grid-points) 


# %{
# # Differentiation matrices
# We need to compute the derivatives in $x$ and also in $y$. 
# But in fact as is usually done for stability of parallel flows, 
# we can do a Fourier transform in the direction where the system does not change, 
# so here the numerical differentiation is done in $y$, and the differentiation in $x$ 
# simply ammounts to multiplication with $ik$.
# %}


y,dy,dyy,w = dif1D('Cheb',a=-1,b=1,N=Nx); 


#
# Definition of the base flow : PLANE POISEUILLE FLOW. We also need its y-derivative.
#

U=1-y**2; Uy=-2*y;




#
# DEFINITION OF FUNCTION solve_eigenproblem
#


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
  start_time = time.time()
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

  epsilon = 1e-8
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
  B[loc,:]=-epsilon*III[loc,:]; 


  #%{
  # Computing the eigenmodes
  #We can use the *eig* function for DAE just as for dynamic systems. the difference is that we will get infinite eigenvalues. Thinking of the incompressible system as the limit of the compressible one, you realize that these infinite eigenvalues are the limit of the sound modes. We remove them, we sort the eigenvalues and eigenvectors.
  #%}

  ass_time = time.time()
  elapsed_time = ass_time - start_time
  print(f"Temps d'assemblage : {elapsed_time:.6f} secondes")

  #% computing eigenmodes 
  s,em = eigs(A,M=B,k=10,sigma=1+1j);
  #print(s)
  end_time = time.time()
  elapsed_time = end_time - ass_time
  print(f"Temps d'inversion methode iterative (k=10) : {elapsed_time:.6f} secondes")

  s,em = eig(A,B);
  #print(s)
  end_time = time.time()
  elapsed_time = end_time - ass_time
  print(f"Temps d'inversion methode directe : {elapsed_time:.6f} secondes")

  
  #% sort and remove spurious eigenvalues
  S = s[abs(s)<1e8];
  EM = em[:,abs(s)<1e8]; 
  order = np.argsort(real(-S))
  S = S[order];
  EM = EM[:,order];
  sort_time = time.time()
  elapsed_time = sort_time-end_time 
  print(f"Temps pour tri des valeurs propres  : {elapsed_time:.6f} secondes")
  
  return S,EM



#
# Definition of two functions for postprocessing (plot mode structure)
#

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



#
# We compute a spectrum for Re = 10,000 ; k = 1
#


k=1;    # the wave number
Re=10000;    # the Reynolds number




S,EM = solve_eigenproblem(Re,k)

#
# We plot the spectrum
#


fig, ax = plt.subplots()
h, = plt.plot(-np.imag(S), np.real(S), 'r+',picker=5)
plt.xlabel('\omega_r');plt.ylabel('\omega_i');plt.ylim([-.5,.05])
fig.canvas.mpl_connect('pick_event', on_pick) # to associate the previously defined function to click events
ax.set_title('Spectrum for Re = '+str(Re)+' ; k = '+str(k)+'\n (Click on eigenvalues to display corresponding modes)')
plt.show()


# We plot the leading eigenmode
plt.figure(1);
plotmode(1,EM[:,0],S[0],k)
plt.show()


# We do a loop over k
for k in range(10):
  S,EM = solve_eigenproblem(Re,k)
  

