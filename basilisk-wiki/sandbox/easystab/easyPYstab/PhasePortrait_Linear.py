# -*- coding: utf-8 -*-
"""
Created on Nov 28 2025

@author: D. Fabre (adapted from a previous Matlab program)


This program is designed to draw phase portraits of a linear homogeneous system
with the form 
$$ 
d/dt [X1;X2]  = [[A11,A12];[A21,A22]] [X1;X2] 
$$


* Remarque importante :
  Pour utiliser ce programme de manière optimale il faut modifier le mode graphique de Python.
  Pour cela (avec Spyder) :
  - Ouvrir le menu "Preferences / Ipython Console / Graphics"
  - Dans le menu "graphics backend" sélectionner "automatic"'
"""

import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
plt.close("all")

"""
#  Selection of parameters
"""


print("Trace le champ de vecteur et les trajectoire du flot défini par :")
print("     u = A11 x + A12 y")
print("     v = A21 x + A22 y\n")
# Valeur des paramètres a, b, c, d :
a = float(input("A11 = "))
b = float(input("A12 = "))
c = float(input("A21 = "))
d = float(input("A22 = "))
A11 = a;
A12 = b;
A21 = c;
A22 = d;

"""
# Trace du flot du systeme dynamique sous forme de champ de vecteur
"""


# Construction d'une grille de 21x21 points dans ces limites
xmin = -2; xmax = 2; ymin = -2; ymax = 2;
[xG, yG] = np.meshgrid(np.linspace(xmin, xmax, 21), np.linspace(ymin, ymax, 21))

# Calcul des composantes x et y du champ de vecteurs à partir des coefficients A11, A12, A21 et A22
ux = A11*xG + A12*yG
uy = A21*xG + A22*yG

# Affichage du champ de vecteurs
fig,ax = plt.subplots()
ax.quiver(xG, yG, ux, uy, color='b')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

# Ajout d'un titre et d'une légende
#ax.set_title('Illustration of a 2D linear flow.\nLeft-click to draw a trajectory; right-click to draw the deformation of a square')
ax.set_xlabel('x')
ax.set_ylabel('y')

# Affichage de la figure
plt.show()


# Définition de certains paramètres
Tmax = 10; # temps maximal pour l'intégration de la trajectoire
Tmaxplot = 5 ;# temps maximal pour l'affichage dans la figure 2
Tsquare = 0.2; # temps de déformation du carré
h = 1 ;# taille du carré
xA = [-h, -h]
xB = [-h, h]
xC = [h, h]
xD = [h, -h]
button = 1
compteur = 0


"""
# Fonction qui sera appelee lorsqu'on clique sur la figure
"""
def on_click(event):
    if event.button == 1:
        # Left-click: drawing of a few trajectories to illustrate the flow
        xinit = np.array([event.xdata, event.ydata])
        plt.figure(1)
        t = np.linspace(0, Tmax, 1000) #the future
        xtraj = odeint(lambda x, t: [A11*x[0]+A12*x[1], A21*x[0]+A22*x[1]], xinit, t)
        plt.plot(xtraj[:,0], xtraj[:,1], 'r', xtraj[0,0], xtraj[0,1], 'ro')
        plt.xlim((xmin,xmax));plt.ylim((ymin,ymax))
        plt.pause(0.01);
        plt.show()
    if event.button == 3:
        # Right-click: backward trajectories
        xinit = np.array([event.xdata, event.ydata])
        plt.figure(1)
        t = np.linspace(0, -Tmax, 1000) #the past
        xtraj = odeint(lambda x, t: [A11*x[0]+A12*x[1], A21*x[0]+A22*x[1]], xinit, t)
        plt.plot(xtraj[:,0], xtraj[:,1], 'r:')
        plt.xlim((xmin,xmax));plt.ylim((ymin,ymax))
        plt.pause(0.01);
        plt.show()
        
# Association de la fonction on_click à la figure
fig.canvas.mpl_connect('button_press_event', on_click)


plt.figure(1)
plt.title('Phase portrait of the linear system dX/dt = A X (Click to draw trajectories !) \n A = [[{},{}],[{},{}]]'.format(A11, A12, A21, A22)) 
#     ' Left-click to draw a trajectory ; right-click to stop'},'FontSize',14)


plt.savefig('Ecoulement2D_Linear.png', dpi=80)
plt.show()


A = np.array([[A11, A12],
              [A21, A22]])

print("Matrix A :")
print(A)

detA = A11*A22 - A12*A21
traceA = A11 + A22
print(f"Det(A) = {detA}")
print(f"Tr(A)  = {traceA}")
print()

# eigenvalues + eigenvectors
eigvals, eigvecs = np.linalg.eig(A)
l1, l2 = eigvals[0], eigvals[1]
e1 = eigvecs[:,0]
e2 = eigvecs[:,1]

print(f"Eigenvalues : lambda1 = {l1} ; lambda2 = {l2}")
print()
print("Eigenmodes :")
print("e1 =", e1)
print("e2 =", e2)
print()

"""
# Classification of equilibrium points according to eigenvalues/
"""

if (np.imag(l1)==0):
### Two real eigenvalues
    if(l1<0) and (l2<0):
      if (l1!=l2):
        print('Two real, negative, distinct eigenvalues : this is a stable node');
        plt.plot([-e1[0],e1[0]],[-e1[1],e1[1]],'m--');
        plt.plot([-e2[0],e2[0]],[-e2[1],e2[1]],'m--');
      elif  (l1==l2) and (abs(A12)+abs(A21)>0):
        print('Two real, negative, equal eigenvalues, nondiagonal matrix : this is an improper stable node');
        plt.plot([-e1[0],e1[0]],[-e1[1],e1[1]],'m--');
      elif (l1==l2) and (abs(A12)+abs(A21)==0):
        print('Two real, negative, equal eigenvalues, diagonal matrix : this is a stable star');
        plt.plot([-1,1],[0,0],'m--');plt.plot([0,0],[-1,1],'m--');
        plt.plot([-np.sqrt(.5),np.sqrt(.5)],[-np.sqrt(.5),np.sqrt(.5)],'m--');plt.plot([-np.sqrt(.5),np.sqrt(.5)],[np.sqrt(.5),-np.sqrt(.5)],'m--');
    elif(l1>0) and (l2>0): 
      if (l1!=l2):
        print('Two real, positive, distinct eigenvalues : this is an unstable node');
        plt.plot([-e1[0],e1[0]],[-e1[1],e1[1]],'r--');
        plt.plot([-e2[0],e2[0]],[-e2[1],e2[1]],'r--');
      elif (l1==l2) and (abs(A12)+abs(A21)>0):
        print('Two real, positive, equal eigenvalues, nondiagonal matrix : this is an improper unstable node');
        plt.plot([-e1[0],e1[0]],[-e1[1],e1[1]],'r--');
      elif (l1==l2) and (abs(A12)+abs(A21)==0):
        print('Two real, positive, equal eigenvalues, diagonal matrix : this is an unstable star');
        plt.plot([-1,1],[0,0],'r--');plt.plot([0,0],[-1,1],'r--');
        plt.plot([-np.sqrt(.5),np.sqrt(.5)],[-np.sqrt(.5),np.sqrt(.5)],'r--');plt.plot([-np.sqrt(.5),np.sqrt(.5)],[np.sqrt(.5),-np.sqrt(.5)],'r--');
    elif(l1*l2<0):
       print('Two real eigenvalues with opposite sign : this is a saddle');
       if(l1<0):
           plt.plot([-e1[0],e1[0]],[-e1[1],e1[1]],'m--');
           plt.plot([-e2[0],e2[0]],[-e2[1],e2[1]],'r--');
       else:
           plt.plot([-e1[0],e1[0]],[-e1[1],e1[1]],'r--');
           plt.plot([-e2[0],e2[0]],[-e2[1],e2[1]],'m--');
    elif(l2==0) and (l1>0):
        print('one null eigenvalue and one positive eigenvalue : this is an unstable degenerated node')
        plt.plot([-e1[0],e1[0]],[-e1[1],e1[1]],'g--');
        plt.plot([-e2[0],e2[0]],[-e2[1],e2[1]],'r--');
    elif(l2==0) and (l1<0):
        print('one null eigenvalue and one negative eigenvalue : this is a nonhyperbolic point of codimension 1')
        plt.plot([-e1[0],e1[0]],[-e1[1],e1[1]],'g--');
        plt.plot([-e2[0],e2[0]],[-e2[1],e2[1]],'m--');
    elif(l1==0) and (l2==0) and (abs(A12)+abs(A21)==0):
        print('two null eigenvalues, diagonal matrix : this is a nonhyperbolic point of codimension 2')
        plt.plot([-1,1],[0,0],'g--');plt.plot([0,0],[-1,1],'g--');
        plt.plot([-np.sqrt(.5),np.sqrt(.5)],[-np.sqrt(.5),np.sqrt(.5)],'g--');plt.plot([-np.sqrt(.5),np.sqrt(.5)],[np.sqrt(.5),-np.sqrt(.5)],'g--');
    elif(l1==0) and (l2==0) and (abs(A12)+abs(A21)>0):
        print('two null eigenvalues, non diagonal matrix : this is an algebraically unstable point')
        plt.plot([-e1[0],e1[0]],[-e1[1],e1[1]],'g--');
        plt.plot([e1[1],-e1[1]],[-e1[0],e1[0]],'o--');
else:
        ### Two complex eigenvalues
    er = np.real(e1);ei = np.imag(e1);
    if(A11==A22) and (A21==-A12): # this is a circular focus/center
        if(np.real(l1)>0):
            print('two complex conjugate eigenvalues with positive real part, all directions equivalent : this is an unstable circular focus')
            plt.plot([-1,1],[0,0],'r-.');plt.plot([0,0],[-1,1],'r-.');
            plt.plot([-np.sqrt(.5),np.sqrt(.5)],[-np.sqrt(.5),np.sqrt(.5)],'r-.');plt.plot([-np.sqrt(.5),np.sqrt(.5)],[np.sqrt(.5),-np.sqrt(.5)],'r-.');
        elif(np.real(l1)<0):
            print('two complex conjugate eigenvalues with negative real part, all directions equivalent : this is a stable circular focus')
            plt.plot([-1,1],[0,0],'m-.');plt.plot([0,0],[-1,1],'m.');
            plt.plot([-np.sqrt(.5),np.sqrt(.5)],[-np.sqrt(.5),np.sqrt(.5)],'m-.');plt.plot([-np.sqrt(.5),np.sqrt(.5)],[np.sqrt(.5),-np.sqrt(.5)],'m-.');
        elif(np.real(l1)==0):
            print('two complex conjugate eigenvalues with null real part :  all directions equivalent : this is a circular centre')
            plt.plot([-1,1],[0,0],'g-.');plt.plot([0,0],[-1,1],'g.-');
            plt.plot([-np.sqrt(.5),np.sqrt(.5)],[-np.sqrt(.5),np.sqrt(.5)],'g-.');plt.plot([-np.sqrt(.5),np.sqrt(.5)],[np.sqrt(.5),-np.sqrt(.5)],'g-.');
    else: # elliptical focus/center
       phi= -math.atan2(2*np.vdot(np.conj(er),ei),(np.vdot(np.conj(er),er)-np.vdot(np.conj(ei),ei)))/2;
       ec = er*np.cos(phi)-ei*np.sin(phi);
       es = -er*np.sin(phi)-ei*np.cos(phi);
       print('main directions for an elliptical focus/center:')
       ec
       es
       if(np.real(l1)>0):
           print('two complex conjugate eigenvalues with positive real part : this is an unstable elliptical focus')
           plt.plot([-ec[0],ec[0]],[-ec[1],ec[1]],'r-.');
           plt.plot([-es[0],es[0]],[-es[1],es[1]],'r-..');  
       elif(np.real(l1)<0):
           print('two complex conjugate eigenvalues with negative real part : this is an stable elliptical focus')
           plt.plot([-ec[0],ec[0]],[-ec[1],ec[1]],'m-.');
           plt.plot([-es[0],es[0]],[-es[1],es[1]],'m-..');  
       elif(np.real(l1)==0):
           print('two complex conjugate eigenvalues with null real part : this is an elliptical centre')
           plt.plot([-ec[0],ec[0]],[-ec[1],ec[1]],'g-.');
           plt.plot([-es[0],es[0]],[-es[1],es[1]],'g-..');  



# Exemple de tracé des directions propres (si utile)
#plt.figure(1)
#origin = np.array([[0, 0],[0, 0]])  # origine pour quiver
#vectors = np.stack((e1, e2), axis=1)
#plt.quiver(origin[0], origin[1], vectors[0], vectors[1],
#           angles='xy', scale_units='xy', scale=1, color=['r','b'])
#plt.axis('equal')
#plt.grid(True)
#plt.show()


