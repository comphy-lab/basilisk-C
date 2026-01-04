# -*- coding: utf-8 -*-
"""
Created on Nov 28 2025

@author: D. Fabre (adapted from a previous Matlab program)


* How to use this program :
  - Launch the program (in spyder, after switching graphical mode to "automatic")
  - Specify the "problem type" you want to consider, among a list of predefined cases
  - Set the parameters of the problem (type enter to use default value)
  - Click on the figures to draw trajectories in phase space 
    (left-click = forward trajectory, right-click = backward trajectroy)

* Remarque importante :
  Pour utiliser ce programme de maniÃ¨re optimale il faut modifier le mode graphique de Python.
  Pour cela (avec Spyder) :
  - Ouvrir le menu "Preferences / Ipython Console / Graphics"
  - Dans le menu "graphics backend" sÃ©lectionner "automatic"'
    
"""
import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
plt.close("all")


def myinput(question,default):
    # Fonction pour selection de parametre avec valeur par defaut
    answer = input(f"{question} [default : {default}] : " )
    if not answer:
        return default 
    else:
        return answer

"""
# Selection of problem type
"""

print("#####################################################\n")
print("  Phase portrait generator for 2D dynamical systems\n")
print("#####################################################\n")
print("Select a dynamical system among the following list :")
print('    "Pendulum" \n    "Brusselator" \n    "VanDerPol" ')
print('    (or add you own case in the program) \n') 
problemtype = str(myinput("Your case : ","Pendulum"))


print(">>> Your choice : ",problemtype)
print(">>> Selection of parameters :")
Problem = {"type": problemtype} ;

"""
###############################################################################

# Big function defining the available problems

###############################################################################
"""


def F_PortraitdePhase(P,X):
  """
  this function defines the function F for the RHS of the dynamical system
  Usage of this function :
  At first call, Problem = F_PortraitdePhase(P,0) 
                 the first argument is expected as a dictionary , the second is unused 
                 and the output is a dictionary containing the parameters of the problem  
  At next calls, : F = F_PortraitdePhase(Problem,X)
                  the first argument is the dictionary, the second is vector X (dimension 2) 
                  and the output is the dynamical function F
  """

  x1 = X[0]; x2 = X[1];
  """    
  Selection according to the problem type      
  """    
  if (P["type"]=="Pendulum"):  
     """    
     ## Case "Pendulum" : basic damped pendulum
     """        
     if not "Init" in P.keys():
           # at first call : define parameters and ranges for figure
           P["mu"] = float(myinput("damping parameter mu = ",.3));     # damping parameter
           P["omega0"] = float(myinput("frequency omega0",1));  # frequency
           P["xmin"] = float(myinput("xmin",-3*np.pi)); 
           P["xmax"] = float(myinput("xmax",3*np.pi)); 
           P["ymin"] = float(myinput("ymin",-2*np.pi)); 
           P["ymax"] = float(myinput("ymax",2*np.pi));
           P["Init"] = "OK";
           return P
     else: 
           # Equation for the motion of a damped pendulum.
           # 
           #  m L^2 \ddot \theta = - \mu_f \dot \theta - m g \sin \theta
           #
           # Considering this problem as a dynamical system of order 2 for
           # $X = [x_1 ; x_2] = [\theta, \dot \theta]$ leads to :
           # 
           # d/dt [x_1 ; x_2]  = [ x_2 ; - \mu x_2 - \omega_0^2 \sin(x_1) ]
           # 
           # with $\omega_0^2 = \sqrt{\frac{g}{L}}$ ; $\mu = \frac{\mu_f}{M L^2}$.
           F = [x2 , -P["omega0"]**2*np.sin(x1)-P["mu"]*x2];
           return F
       
  elif (P["type"]=="VanDerPol"): # VanDerPol is well-known model of self-sustained oscillator
       """    
       ## Case "VanDerPol" : a classical model of nonlinear oscillator    
       """    
       if not "Init" in P.keys():
           # at first call : define parameters and ranges for figure
           P["omega0"] = float(myinput("frequency omega0",1)); # frequency
           P["r"] = float(myinput("damping rate r ",0.3));    # coef of linear term
           P["delta"] = float(myinput("nonlinear coeff delta",1));  # coef of nonlinear term
           print("Range for figure:")
           P["xmin"] = float(myinput("xmin",-3)) ; 
           P["xmax"] = float(myinput("xmax",3)); 
           P["ymin"] = float(myinput("ymin",-3)); 
           P["ymax"] = float(myinput("ymax",3));
           P["Init"] = "OK";
           return P
       else:
            # The VanDerPol oscillator is a well-known model of self-sustained
            # oscillator. It is defined as follows :
            #
            # $$ m \ddot x +\omega_0^2 x   = (r - \delta x^2) \dot x$$
            #
            F = [x2 , -P["omega0"]**2*x1+(P["r"]-P["delta"]*x1**2)*x2];
            return F
        
  elif (P["type"]=="Brusselator"):  
      """    
      ## Case "Brusselator", a model for funny chemical reactions
      """ 
      if not "Init" in P.keys():
          P["alpha"] = float(myinput("parameter alpha",1));
          P["beta"] = float(myinput("parameter beta",1));
          print("Range for figure:")
          P["xmin"] = float(myinput("xmin",0)) ; 
          P["xmax"] = float(myinput("xmax",6)); 
          P["ymin"] = float(myinput("ymin",0)); 
          P["ymax"] = float(myinput("ymax",6));
          P["Init"] ="OK";
          return P
      else:
        # 
        #  The 'Brusselator' (more precisely the homogeneous brusselator) is a simple model
        #  for a chemical reaction displaying unsteadiness. It is defined as follows :
        #
        # $$ \dot x_1 = -(\beta + 1) x_1 + x_1^2 x_2 + \alpha  ;$$
        # $$\dot x_2 = \beta x_1 - x_1^2 x_2 ;$$
        #
        F = [-(P["beta"]+1)*x1+x1**2*x2+P["alpha"],  P["beta"]*x1-x1**2*x2] ;
        return F

"""
###############################################################################
#
#  Main program
#
###############################################################################
"""


"""
## Trace de la figure et du champ de vecteur "flow"
"""
# Initialisation, to set the parameters of the problem
Problem = F_PortraitdePhase( {"type": problemtype} ,[0,0])

# Construction d'une grille de 21x21 points dans ces limites
xP = np.linspace(Problem["xmin"], Problem["xmax"], 21)
yP = np.linspace(Problem["ymin"], Problem["ymax"], 21)
[xG, yG] = np.meshgrid(xP,yP);
ux = np.zeros((21,21)); uy = np.zeros((21,21));

# Calcul du champ de vecteur sur cette grille pour visualisation sur figure
for i in range(len(xP)):
   for j in range(len(yP)):
       F = F_PortraitdePhase(Problem,[xP[i],yP[j]]);
       ux[j,i] = F[0];
       uy[j,i] = F[1];
      
# Affichage du champ de vecteurs
fig,ax = plt.subplots()
ax.quiver(xG, yG, ux, uy, color='b')
ax.set_xlim(Problem["xmin"], Problem["xmax"])
ax.set_ylim(Problem["ymin"], Problem["ymax"])

# Ajout d'un titre et d'une legende
plt.title('Phase portrait of the model problem "'+Problem["type"]+ '"') 
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()

# Definition de certains parametres
Tmax = 100; # temps maximal pour l'integration de la trajectoire
Tmaxplot = 5 ;# temps maximal pour l'affichage dans la figure 2
Nstep = 10000 ; # number of steps in integration
button = 1
compteur = 0

"""
## Fonction qui sera appelee lorsqu'on clique sur la figure
"""
def on_click(event):
    if event.button == 1:
        # Left-click: compute a forward (future) trajectory and draw with plain line
        xinit = np.array([event.xdata, event.ydata])
        plt.figure(1)
        t = np.linspace(0, Tmax, Nstep) 
        xtraj = np.zeros((2,Nstep))
        xtraj = odeint(lambda x, t: F_PortraitdePhase(Problem,x), xinit, t)
        plt.plot(xtraj[:,0], xtraj[:,1], 'r', xtraj[0,0], xtraj[0,1], 'ro')
        plt.pause(0.01);
        plt.show()
    if event.button == 3:
        # Right-click: compute a backwards (past) trajectory and draw with dotted line
        t = np.linspace(0, -Tmax, Nstep)
        xinit = np.array([event.xdata, event.ydata])
        plt.figure(1)
        xtrajB = np.zeros((2,Nstep)) 
        xtrajB = odeint(lambda x, t: F_PortraitdePhase(Problem,x), xinit, t)
        plt.plot(xtrajB[:,0], xtrajB[:,1], 'r:')
        plt.pause(0.01);
        plt.show()
        
## Association de la fonction on_click a la figure
fig.canvas.mpl_connect('button_press_event', on_click)

plt.savefig('PhasePortrait.png', dpi=80)
plt.show()



"""
###############################################################################

#  APPENDIX: the original function in Matlab with lots of other cases to reintroduce...
 
###############################################################################
"""

"""
function F = F_PortraitdePhase(X)

% Fonction pour tracer le portrait de phase d'un systeme 2-2 de la forme
% d X /d T = F(X).
%
global typeproblem xmin xmax ymin ymax;
x1 = X(1); x2 = X(2);



switch (typeproblem)

%{

### Case "Pendulum" :

Equation for the motion of a damped pendulum.

$$
  m L^2 \ddot \theta = - \mu_f \dot \theta - m g \sin \theta
$$

Considering this problem as a dynamical system of order 2 for
$X = [x_1 ; x_2] = [\theta, \dot \theta]$ leads to :

$$\frac{ d}{dt} [x_1 ; x_2]  = [ x_2 ; - \mu x_2 - \omega_0^2 \sin(x_1) ]$$

with $\omega_0^2 = \sqrt{\frac{g}{L}}$ ; $\mu = \frac{\mu_f}{M L^2}$.

**Exercice :**
Observe the structure of the phase portrait for various
values of the damping parameter mu. What do we observe for mu=0 ?

%}
    case('Pendulum')

        mu = 0;
        omega0 = 1; % omega_0 = sqrt(g/L)
        F = [x2;-omega0^2*sin(x1)-mu*x2];
        % adjust the figure axes
        xmin = -3*pi ; xmax = 3*pi; ymin = -2*pi; ymax = 2*pi;

%{

### Case "RotatingPendulum" :

Equation for the motion of a pendulum with imposed precession

$$
  m L^2 \ddot \theta = - \mu_f \dot \theta - m g L \sin \theta + m L^2 \Omega^2 \sin \theta \cos \theta
$$

Considering this problem as a dynamical system of order 2 for
$X = [x_1 ; x_2] = [\theta, \dot \theta]$ leads to :

$$\frac{ d}{dt} [x_1 ; x_2]  = [ x_2 ; - \mu x_2 - \omega_0^2 \sin x_1 (1 - R \cos x_1 )  ]$$

with $\omega_0^2 = \sqrt{\frac{g}{L}}$ ; $\mu = \frac{\mu_f}{M L^2}$ ; $R = \frac{\Omega^2}{\omega^2}$.

**Exercice :**
Observe the structure of the phase portrait for various values of the rotation parameter R.
Draw the corresponding bifurcation diagram.

%}
    case('RotatingPendulum')

        omega0 = 1; mu = .5;
        R = 1.1;
        F = [x2;-omega0^2*sin(x1)*(1-R*cos(x1))-mu*x2];
        % adjust the figure axes
        xmin = -pi/2 ; xmax = pi/2; ymin = -.5; ymax =.5;

        %{

### Case "InvertedPendulum" :

Equation for the motion of an inverted pendulum with a spring (offset angle alpha)

$$
  m L^2 \ddot \theta = - \mu_f \dot \theta + m g L \sin \theta - K (\theta-\alpha)
$$

Considering this problem as a dynamical system of order 2 for
$X = [x_1 ; x_2] = [\theta, \dot \theta]$ leads to :

$$\frac{ d}{dt} [x_1 ; x_2]  = [ x_2 ; - \mu x_2 + \omega_0^2 \sin x_1 - k ( x_1- \alpha)    ]$$

with $\omega_0^2 = \sqrt{\frac{g}{L}}$ ; $\mu = \frac{\mu_f}{m L^2}$ ; $k = \frac{K}{m L^2}$.

**Exercice :**
Observe the structure of the phase portrait for various values of the offset angle alpha.
Draw the corresponding bifurcation diagram.

%}
    case('InvertedPendulum')

        omega0 = 1; mu = .5;k = .9;
        alpha = -0.025;

        F = [x2;+omega0^2*sin(x1)-k*(x1-alpha)-mu*x2];
        % adjust the figure axes
        xmin = -pi/2 ; xmax = pi/2; ymin = -.5; ymax =.5;



%{
### Case 'VanDerPol' :

The VanDerPol oscillator is a well-known model of self-sustained
oscillator. It is defined as follows :

$$ m \ddot x +\omega_0^2 x   = (r - \delta x^2) \dot x$$


**Exercice :**
Observe the structure of the phase portrait for various choices of the parameters.
Reconstruct the corresponding bifurcation diagram in terms of the bifurcation parameter r.

%}
    case('VanDerPol');
        omega0 = 1;
        r = 0.5;
        delta = -1;
        F = [x2 ; -omega0^2*x1+(r-delta*x1^2)*x2];
        % adjust the figure axes
        xmin = -3 ; xmax = 3; ymin = -3; ymax = 3;



%{
### Case 'Brusselator' :

The 'Brusselator' (more precisely the homogeneous brusselator) is a simple model
for a chemical reaction displaying unsteadiness. It is defined as follows :

$$
\dot x_1 = -(\beta + 1) x_1 + x_1^2 x_2 + \alpha  ;
$$
$$
\dot x_2 = \beta x_1 - x_1^2 x_2 ;
$$



**Exercice :**
Observe the structure of the phase portrait for various choices of the parameters.
Construct the bifurcation diagram as function of the bifurcation parameter beta.

%}
    case('Brusselator')

        alpha = 1;
        beta = 2.2;
        F = [-(beta+1)*x1+x1^2*x2+alpha ; beta*x1-x1^2*x2] ;


        % adjust the figure axes
        xmin = 0 ; xmax = 6; ymin = 0; ymax = 6;


%{
### Case 'LotkaVolterra' :

The Lotka-Volterra is a simple system representing the competition of two
animal species (for instance hares and lynxes). It is defined as follows :

$$
\dot x_1 = x_1( \alpha - \beta x_2) ;
$$
$$
\dot x_2 = x_2( -\gamma + \delta x_1) ;
$$



**Exercice :**
Observe the structure of the phase portrait for various choices of the parameters.
Is this problem conservative or not ?

%}
    case('LotkaVolterra')
        % we take the same values as wikipedia...
        alpha = 2/3;
        beta = 4/3;
        gamma = 1;
        delta = 1;
        F = [(alpha-beta*x2)*x1 ; (delta*x1-gamma)*x2];
        % adjust the figure axes
        xmin = 0 ; xmax = 2.5; ymin = 0; ymax = 2.5;

        %{
### Case 'BuffaloWolf' :

This problem is a

$$
\dot x_1 = r x_1 - A x_1 (x_2+ E x_2^2) - B x_1^2;
$$
$$
\dot x_2 = -C x_2 + D x_1 (x_2+ E x_2^2) ) ;
$$



**Exercice :**
Observe the structure of the phase portrait for various choices of the parameters.
Reconstruct the corresponding bifurcation diagram as the parameter r is
varied

%}
    case('BuffaloWolf')
    r = .4;
    A = .3;
    B = 0.1;
    C = 1;
    D = .2;
    E = 0.15;
    F = [(r*x1 -A*x1*(x2+E*x2^2) - B*x1^2) ; (-C*x2 + D*x1*(x2+E*x2^2)) ];

    % adjust the figure axes
  xmin = -5 ; xmax = 20; ymin = -5; ymax = 10;

        %{
### Case 'Trefethen' :

This is a simple model which is asymptotically stable but may
  lead to sucritical transition with very small initial conditions thanks to
  transient growth.





$$
\dot x_1 = -(1/R) x_1  + x_2 -\sqrt{x_1^2+x_2^2}  x_2;
$$
$$
\dot x_2 = -(2/R) x_2 + \sqrt{x_1^2+x_2^2} x_1;
$$



**Exercice :**
observe the behaviour of the system for very small values of the initial
  condition (try with various values of $R$?)

%}
    case('Trefethen')
    R = 10;
    NL = 1;% select 0 or 1
    normX = sqrt(x1^2+x2^2)  ;
    F1 = -(1/R)*x1 + x2 - NL*normX*x2;
    F2 = -(2/R)*x2 + NL*normX*x1;
    F = [F1;F2];
  % adjust the figure axes
  xmin = -1 ; xmax = 1; ymin = -1; ymax = 1;

        %{
### Case 'SH_5.1' :

This case is a variant of Trefethen's original model,
and is taken from Schmid & Henningson (Exercice 5.1, p. 525)

%}
    case('SH_5.1')
    R = 10;
    NL = 0;% select 0 or 1
    normX = sqrt(x1^2+x2^2)  ;
    F1 = -(1/R)*x1 + NL*normX*x2;
    F2 = -(2/R)*x2 + x1 - NL*normX*x1;
    F = [F1;F2];
  % adjust the figure axes
  xmin = -1 ; xmax = 1; ymin = -1; ymax = 1;

  %{
### Case 'SaddleNode' :

We investigate a saddle-node bifurcation in dimension 2 by considering the
following model:

$$ m \ddot x + \mu \dot x  = r - x^2 $$

If the friction is dominant over the inertia this equation reduces to the classical
one-dimensional saddle-node bifurcation.

**Exercice :**
Observe the structure of the phase portrait for various
choices of the parameters. Reconstruct the corresponding bifurcation diagram.

%}
    case('SaddleNode')
          mass = .5; friction = 1; gravity = 1;
          r=  -.1;
          F = [x2 ; -friction*x2/mass + (r-x1^2)*gravity];
        % adjust the figure axes
        xmin = -3 ; xmax = 3; ymin = -1; ymax = 1;
  %{
### Case 'Pitchfork' :

We investigate a transcritical bifurcation in dimension 2 by considering the
following model:

$$ m \ddot x + \mu \dot x  = r x  - \delta x^3 $$

If the friction is dominant over the inertia this equation reduces to the classical
one-dimensional pitchfork bifurcation.

**Exercice :**
Observe the structure of the phase portrait for various choices of the parameters.
Reconstruct the corresponding bifurcation diagram.

%}
    case('Pitchfork')
          mass = .5; friction = 2;
          r= .5; delta = 1;
          F = [x2 ; -friction*x2/mass+(r*x1-delta*x1^3)/mass];
          % adjust the figure axes
          xmin = -1 ; xmax = 1; ymin = -1; ymax = 1;

%{
### Case 'TransCritical' :

We investigate a transcritical bifurcation in dimension 2 by considering the
following model:

$$ m \ddot x + \mu \dot x  = r x  - x^2 $$

If the friction is dominant over the inertia this equation reduces to the classical
one-dimensional transcritical bifurcation.

**Exercice :**
Observe the structure of the phase portrait for various choices of the parameters.
Reconstruct the corresponding bifurcation diagram.

%}

    case('TransCritical')
          mass = .5; friction = 1;
          r= .5;
          F = [x2 ; -friction*x2/mass+(r*x1-x1^2)/mass];
          % adjust the figure axes
          xmin = -1 ; xmax = 1; ymin = -.5; ymax = .5;

%{
### Case 'Exo3.1' :
%}

    case('Exo3.1')
          F = [x1*(1-x1) ; x2*(2-4*x1)];
          % adjust the figure axes
          xmin = -.5 ; xmax = 1.5; ymin = -.5; ymax = 1.5;


%{
### Case 'Exam2021' :
%}
    case('Exam2021')
          r = -.9; A = 100; B = 1;
          F = [x1*(r+1-(x1-1)^2)+A*x2 ; -B*x2];
          % adjust the figure axes
          xmin = -.5 ; xmax = 1.6; ymin = -.01; ymax = .02;

    case('Exam2025')
          r = .5; A = 0; 
          F = [x1*(r-x1^2)*(1-x1^2)+A*x2 ; -x2];
          % adjust the figure axes
          xmin = -1.5 ; xmax = 1.5; ymin = -1; ymax = 1;


    case('Custom')
        % add your own case here !

    otherwise
        error(['Error in PhasePortrait_NonLinear : unknown type ', typeproblem]);

end%swith

end%function
"""