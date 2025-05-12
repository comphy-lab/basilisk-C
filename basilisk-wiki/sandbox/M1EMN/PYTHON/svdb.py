#
# Resolution des equations de Saint Venant instionnaires non visqueuses sur fond plat
# Solution de la rupture de barrage
# http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf
# http://basilisk.fr/sandbox/M1EMN/Exemples/svdb.c
#
#
# trucs de python
import numpy as np
import scipy.sparse as sp
import math
import matplotlib.pylab as plt
#parametres 
n    = 400
xmin = -8.
xmax = 8.
L=abs(xmax-xmin)
dx   = L/n
dt=dx*.25
# tableau des x, de 0 a n+1 soit n+2 points (n interieurs)
x=np.zeros(n+2)
for i in range(n+2):
  x[i]=(i-n/2)*dx

# hauteur initiale: retenue d'eau h=1 pour x<0; h=0 pour x>0
  u=np.zeros(n+2)
  h=np.zeros(n+2)
  for i in range(n+2):
      h[i]=1*(x[i]<0)
# definition du flux numerique pour la hauteur
def FR1(ug, ud, hg, hd):
  c=max(abs(ug)+np.sqrt(hg),abs(ud)+np.sqrt(ud))
  return (hg*ug+hd*ud)*0.5-c*(hd-hg)*0.5

# definition du flux numerique pour la quantite de mouvement
def FR2(ug, ud, hg, hd):
  c=max(abs(ug)+np.sqrt(hg),abs(ud)+np.sqrt(ud))
  return (ug*ug*hg + hg*hg/2. + ud*ud*hd + hd*hd/2.)*0.5 - c*(hd*ud-hg*ug)*0.5

# definition de l'avancee au temps t
def sol_num(t) :

  fp=np.zeros(n+2)
  fd=np.zeros(n+2)
  nt=int(t/dt)
  t=0

# boucle en temps
  for k in range(nt):
    t=t+dt
# pour chaque x, calcul des flux mis dasn un tableau 
# le cacul est pour h[1] a h[n] et pour u[1] a u[n]
# de 1 a n+1
    for i in range(1,n+2):
      fp[i]=FR1(u[i-1],u[i],h[i-1],h[i])
      fd[i]=FR2(u[i-1],u[i],h[i-1],h[i])
# avancee pour h, si h>0 avancee pour u      
    for i in range(1,n+1):
      hn=h[i]-dt*(fp[i+1]-fp[i])/dx;
      if(h[i]>0.):
        q=h[i]*u[i]-dt*(fd[i+1]-fd[i])/dx
        u[i]=q/hn
        h[i]=hn
      else:
        u[i]=0.
        h[i]=hn
#condition de neumann en 0 et en n         
    u[0]=u[1]
    h[0]=h[1]
    h[n+1]=h[n]
    u[n+1]=u[n]
    print(t)
  return h

# definition de la solution exacte
def sol_exact(x,t):
    h_exact=(((x-0)<=-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<=2*t))
    return h_exact

#test au temps 3 
h3n=np.zeros(n+2)
h3n=sol_num(3)
h3e=sol_exact(x,3)


# on compare la solution analytique et la solution numerique au temps 3
# verifier que plus n est grand plus la solution numerique est proche de la solution analytique (convergence)
# verifier que dt < dx pour que ce soit stable (condition de CFL)
plt.figure(figsize=(8,6))
plt.plot(x,h3n,'b',linestyle='--',label='num')
plt.plot(x,h3e,'r',label='exact')
plt.title('comparaison solution num et analytique probleme rupture barrage, t=3')
plt.xlabel('x')
plt.ylabel('h')
plt.grid()
plt.legend(loc='best')
plt.show()
