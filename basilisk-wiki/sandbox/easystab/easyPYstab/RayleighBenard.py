"""
# Rayleigh-Benard instability

*This document belongs to the [easystab](http://basilisk.fr/sandbox/easystab/README) project, 
please consult the main page of the project for explanations.*

*The present program was specifically designed as a support for the course 
["Introductions to hydrodynamical instabilities"](http://basilisk.fr/sandbox/easystab/M2DET/Instabilities.md) 
for the "DET" Master's cursus in Toulouse. For the underlying theory please see 
[Lecture notes for chapter 5](http://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md) *

*To use this program: click "raw page source" in the left column, copy-paste in the Matlab/Octave editor, and run.*


This program solves the linear stability problem for the Rayleigh-Bénard instability 
of a horizontal layer of fluid heated from below. This program is adapted from 
[/sandbox/easystab/rayleigh_benard.m](), which solved the problem in the case of slip 
conditions ("Free-Free boundaries") and validated the approach by comparing with analytical 
results which exist in this case. 
The present program considers the more physical case of no-slip conditions at the walls
("Rigid-Rigid boundaries").

The theory is presented in [lecture 5](http://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md) 
of the M2 DET course. The  problem is set in generalised eigenvalue form 
$\lambda B X = A X$. The construction of the matrices and resolution is done by [function RB(k,Ra,Pr)](#function-rb) 

@author: David Fabre 

This program is protected by Creative Commons Licence CC-BY-NC
"""                                                                                                       

import EasyStab as ES                     # coeur de la toolbox easystab
import matplotlib.pyplot as plt           # Pour tracer
import numpy as np                        # Pour calculer
from numpy import*                        # Pour construire les tableaux
from matplotlib.pyplot import*            # Pour construire le graphe
from scipy.linalg import eig              # Pour resoudre pb generalise
from scipy.sparse.linalg import eigs      # Pour resoudre pb generalise avec methode iterative
from scipy.optimize import root_scalar
from scipy.linalg import block_diag
plt.close('all')

#%% Definition of functions

"""

# FUNCTION RB

Here is the definition of the function RB which performs the eigenvalue
computation. Note that this function is designed so that it can be used
in two ways :
    
- [s,U] = RB(k,Ra,P) will return a vector s containing the 10 leading
eigenvalues and an array U containing (in column) the corresponding eigenvectors.
    
    
- lambdamax = RB(k,Ra,P) will return only one value corresponding to the leading eigenvalue. 
This is useful to use this function with fzero.
    
"""

def RB(k,Ra,Pr):

    #% renaming the differentiation matrices
    Z=np.zeros((n,n)); 
    I=np.eye(n); 
    dx=1j*k*I; dxx=-k**2*I;
    Delta=dxx+dyy;
    
    """
    ## System matrices
    
    
    As explained in  the [lecture notes](http://basilisk.fr/sandbox/easystab/LectureNotes_RayleighTaylor.md#case-of-a-horizontal-cell-of-large-dimension), 
    after nondimensionalization the system of equations can be written in a matrix form
    $$
    \lambda B q =A q
    $$
    with the matrices
    $$
    q=\left(\begin{array}{c}
    u \\ v\\ p\\ \theta
    \end{array}\right)
    , \quad
    A=\left(\begin{array}{cccc}
    Pr \Delta&0&-\partial_x&0\\
    0&Pr \Delta&-\partial_y& Pr \\
    \partial_x&\partial_y&0&0\\
    0& Ra &0&\Delta
    \end{array}\right)
    , \quad
    B=\left(\begin{array}{cccc}1&0&0&0\\0&1&0&0\\0&0&0&0\\0&0&0&1\end{array}\right)
    $$
    
    NB in this program we use a slightly different nodimensionalization 
    in which the velocities are resclaed by $\sqrt(Re)$; this is equivalent.
    
    
    """
    
    #% system matrices
    A=np.block( 
      [ [ Pr*Delta, Z,          -sqrt(Ra)*dx, Z] ,
        [ Z,        Pr*Delta,   -sqrt(Ra)*dy, Pr*sqrt(Ra)*I] ,
        [ dx,       dy,         Z,            Z] ,
        [ Z,        sqrt(Ra)*I, Z,            Delta] ] )
    
    B = block_diag(I,I,Z,I);
    
    """
    ## Boundary conditions
    
    The natural conditions are that $u$ and $v$ are zero at the walls, 
    and the temperature perturbation $\theta$ is also zero at the walls 
    (the temperature is imposed at the wall, without perturbations). 
    
    $$
    \begin{array}{l}
    u|_0=0\\
    u|_L=0\\
    v|_0=0\\
    v|_L=0\\
    \theta|_0=0\\
    \theta|_L=0\\
    \end{array}
    $$
    thus the boundary conditions are expressed $Cq=0$, with the constraint matrix
    
    $$
    C=\left(\begin{array}{cccc}
    I|_0&0&0&0\\
    I|_L&0&0&0\\
    0&I|_L&0&0\\
    0&I|_0&0&0\\
    0&0&0&I|_L\\
    0&0&0&I|_0\\
    \end{array}\right)
    $$
    
    """
    
    #% boundary conditions
    II=eye(4*n); 
    u0=0; uL=n-1; v0=n; vL=2*n-1; T0=3*n; TL=4*n-1;
    loc=(u0,uL,v0,vL,T0,TL);
    C=II[loc,:];  A[loc,:]=C;
    B[loc,:]=0; 
    
    """
    ## Computing eigenmodes
    
    In the way the linear system $\lambda B q=Aq$ hzs been constucted, the matrix $B$ is not invertible, so we solve a generalized eigenvalue problem which will have infinite eigenvalues. We remove them form the results.  
    """
    
    #% computing eigenmodes
    s,em = eig(A,B);
    
    #% sort eigenalues
    S = s[abs(s)<1e8];
    EM = em[:,abs(s)<1e8]; 
    order = np.argsort(real(-S))
    S = S[order];
    EM = EM[:,order];
    return S,EM    

# function returning only leading eigenvalue (useful for building neutral curve)
def leading_eigenvalue(k,Ra,Pr):
   S,EM = RB(k,Ra,Pr)
   return S[0]

# Function to plot eigenmodes
def plotmode(event,mode,S,k):
    # This function is used to represent the structure of an eigenmode
    global y,dy,dyy
    n = y.shape[0];
    u = mode[0:n];
    v = mode[n:2*n];
    p = mode[2*n:3*n];
    T = mode[3*n:4*n];
    MM =max(abs(mode))
    vorticity = -(dy@u)+1j*k*v;
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 2]})
    ax1.plot(real(u),y,'b-',label='Re(u)');ax1.plot(imag(u),y,'b--',label='Im(u)');
    ax1.plot(real(v),y,'g-',label='Re(v)');ax1.plot(imag(v),y,'g--',label='Im(v)');
    ax1.plot(real(p),y,'k-',label='Re(p)');ax1.plot(imag(p),y,'k--',label='Im(p)');
    ax1.plot(real(T),y,'r-',label='Re(T)');ax1.plot(imag(T),y,'r--',label='Im(T)');
    ax1.set_ylabel('y');ax1.set_xlim(-MM*1.2,MM*1.2); 
    ax1.set_ylim(0,1); 
    ax1.legend(loc='lower center',fontsize=9,ncol=2)
    if (-imag(S)/k)<1 and (-imag(S)/k)>0:
      ycrit = sqrt(1+imag(S)/k);
      ax1.plot([-MM,MM],[ycrit,ycrit],'k:',[-MM,MM],[-ycrit,-ycrit],'k:')    
    Lx=2*pi/k; Nx =30;  x=linspace(-Lx/2,Lx/2,Nx);
    yy = y;
    pp = 2*real(np.outer(p,exp(1j*k*x)));
    uu=2*real(np.outer(u,exp(1j*k*x)));
    uu
    vv=2*real(np.outer(v,exp(1j*k*x)));
    TT=2*real(np.outer(T,exp(1j*k*x)));
    vorticityvorticity=2*real(np.outer(vorticity,exp(1j*k*x)))
    ax2.contourf(x,yy,TT,10); 
    ax2.quiver(x,yy,uu,vv,color='k');
    ax2.set_xlabel('x'); 
    z_str = f"{S.real:.3f} + {S.imag:.3f}j"
    plt.title("Structure of the eigenmode \n (Temperature component + velocity)\n" +" lambda = "+z_str,fontsize=10);

# Function to link the previous functin with the figure 1 
def on_pick(event):
    global EM
    ind = event.ind[0]  # Index du point cliqué
    EMM =EM[:,ind]
    z_str = f"{S[ind].real:.3f} + {S[ind].imag:.3f}j"
    print("Clic sur valeur propre "+str(ind)+" : Eigenvalue = "+z_str+ " ind " +str(ind))
    plotmode(event, EMM, S[ind], k)
    plt.show();plt.pause(0.01);
    



#%% Main program

#% parameters for the discretization
Nx = int(ES.myinput("==> Nx, number of points ? ",60)); # order of discretization
n = Nx+1 # number of grid points
DiscetizationType = ES.myinput("==> Discretization type ? ","Cheb")

#% building differentiation matriCes
y,dy,dyy,w = ES.dif1D(DiscetizationType,parameters={"a":0,"b":1},N=Nx); 


#%% Computing one spectrum for one set of parameters

print("\n Chapter 1 : compute one spectrum and plot the eigenvalues/eigenmodes \n")

# Selecting the parameters
k  = float(ES.myinput("==> Wavenumber k ?    ",3.117));    
Ra = float(ES.myinput("==> Rayleigh number ? ",1708));    
Pr = float(ES.myinput("==> Prandtl number ? ",10));       
# NB : According to the litterature, the neutral conditions are 
#      Ra = 1708 and k = 3.117 whatever the value of Pr

# compute spectrum 
print("   (Computing spectrum)")
S,EM = RB(k,Ra,Pr);

#We plot the spectrum in figure 1
fig, ax = plt.subplots(1)
h, = plt.plot(imag(S),real(S),'go',picker=5);
fig.canvas.mpl_connect('pick_event', on_pick) # to associate the previously defined function to click events
plt.plot([-1, 1],[0,0],'k:');
plt.ylabel('lambda_r');ylim([-200, 50]);
plt.xlabel('lambda_i');xlim([-1,1]);
plt.title('Spectrum  for Ra ='+str(Ra)+' ; k = '+f"{k:.4f}"+' ; Pr = '+str(Pr))
plt.show(); plt.pause(0.01)

# plot the leading eigenmode in figure 1 
plotmode(1,EM[:,0],S[0],k);
plt.figure(2);plt.show();plt.pause(0.1);

print("   => see spectrum in figure 1 and leading mode in figure 2 ")
print("        (click on figure 1 to plot other eigenmodes) ")


#%% Next : parametric study

"""
First we compute $\lambda(k)$ for several values of $Ra$ and plot the
results in figure 3.
"""

rep = ES.myinput("\n   Do you want to draw curves lambda(k) for several Ra ?","yes")
if rep=="yes":
    plt.figure();
    for Ra in linspace(1500,2500,5):
        print("Ra = "+str(Ra)+" :")
        ktab = linspace(1,5,41)
        smaxtab = []
        for k in ktab:
            s,U = RB(k,Ra,Pr);
            smaxtab.append(real(s[0]));
        plt.plot(ktab,smaxtab,label="Ra = "+str(Ra));
        plt.show;plt.pause(0.001);
    plt.plot(ktab,0*ktab,'k:')
    plt.legend();
    plt.xlabel('k');plt.ylabel('lambda')
    

  
#%% Marginal stability curve
  
"""
We now build the marginal stability curve $Ra_c(k)$ 
corresponding to the location in the $[Ra,k]$-plane where the 
leading eigenvalue is exactly zero. 
    
For this we do a loop over k and look for the location where
$\lambda_{max}$ is exactly zero, using function "root_scalar"
Results are plotted in figure 4.    
"""
    
rep = ES.myinput("\n   Do you want to Reconstuct the marginal stability curve in the (Re,k) plane ? ","yes")
if rep=="yes":
    Ra = 3000;
    Ratab=[];
    ktab = linspace(1.5,6,46);
    for k in ktab:
        sol = root_scalar(lambda a: np.real(leading_eigenvalue(k, a,Pr)),
                          x0=Ra, method='newton')
        Ratab.append(sol.root);
        print("k = " +str(k)+" ; Ra = "+str(sol.root))
        
        
    #%% Draw the figure
    plt.figure();
    plt.plot(ktab,Ratab);
    plt.xlabel('k');plt.ylabel('Ra');
    plt.ylim((1000,3000));
    plt.title('Neutral curve');
    plt.text(2, 1300, "Stable")
    plt.text(3.5, 2500, "Unstable")
    plt.show()

