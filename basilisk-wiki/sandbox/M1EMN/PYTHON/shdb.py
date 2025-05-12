# 03/2020
#
# Resolution des equations de Savage Hutter
# aka Saint Venant instationnaires avec frottement de Coulomb sur fond plat
# Solution de la rupture de barrage
# http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf
# http://basilisk.fr/sandbox/M1EMN/Exemples/svdb.c
# http://basilisk.fr/sandbox/M1EMN/Exemples/savagestaron.c
# friction solide (http://basilisk.fr/sandbox/M1EMN/Exemples/granular_sandglass_muw.c)
#
#  voir Kerswell
#
# trucs de python
import numpy as np
import scipy.sparse as sp
import math
import matplotlib.pylab as plt
#parametres
n    = 100
xmin = -4.
xmax = 4.
L=abs(xmax-xmin)
dx   = L/n
dt=dx*.25
# coefficient de friction de Coulomb
mu=0.5
#tableau des x
x=np.zeros(n+2)
for i in range(n+2):
  x[i]=(i-n/2)*dx
h=np.zeros(n+2)
u=np.zeros(n+2)
# hauteur initiale: retenue de grains h=1 pour x<0; h=0 pour x>0
for i in range(n+2):
  h[i]=1*(x[i]<0)

# definition du flux numerique (Rusanov) pour la hauteur
def FR1(ug, ud, hg, hd):
  c=max(abs(ug)+np.sqrt(hg),abs(ud)+np.sqrt(ud))
  return (hg*ug+hd*ud)*0.5-c*(hd-hg)*0.5
# definition du flux numerique (Rusanov) pour la quantite de mouvement
def FR2(ug, ud, hg, hd):
  c=max(abs(ug)+np.sqrt(hg),abs(ud)+np.sqrt(ud))
  return (ug*ug*hg + hg*hg/2. + ud*ud*hd + hd*hd/2.)*0.5 - c*(hd*ud-hg*ug)*0.5


h0c=np.zeros(n+2)
h1c=np.zeros(n+2)
h2c=np.zeros(n+2)
h3c=np.zeros(n+2)
h4c=np.zeros(n+2)
h5c=np.zeros(n+2)



def sauve(t):
  global h0c,h1c,h2c,h3c,h4c,h5c
  if t>=0 and t < 0 + dt:
    h0c=h.copy()
  if t>=1 and t < 1 + dt:
    h1c=h.copy()
  if t>=2 and t < 2 + dt:
    h2c=h.copy()
  if t>=3 and t < 3 + dt:
    h3c=h.copy()
  if t>=4 and t < 4 + dt:
    h4c=h.copy()
  if t>=5 and t < 5 + dt:
    h5c=h.copy()




# definition de l'avancee au temps t
def sol_num(t) :
  global h
  global h0n,h1n,h2n,h3n,h4n,h5n
  fp=np.zeros(n+2)
  fd=np.zeros(n+2)
  nt=int(t/dt)
  t=0
  
# boucle en temps
  for k in range(nt):
    sauve(t)
    t=t+dt
# pour chaque x, calcul des flux mis dans un tableau
    for i in range(1,n+1):
      fp[i]=FR1(u[i-1],u[i],h[i-1],h[i])
      fd[i]=FR2(u[i-1],u[i],h[i-1],h[i])
# avancee pour h, si h>0 avancee pour u
    for i in range(1,n):
      hn=h[i]-dt*(fp[i+1]-fp[i])/dx;
      if(h[i]>0.):
        q=h[i]*u[i]-dt*(fd[i+1]-fd[i])/dx
        u[i]=q/hn
        h[i]=hn
      else:
        u[i]=0.
        h[i]=hn
# friction :
    for i in range(1,n+1):
      un = abs(u[i])
      if(un>0):
        u[i] = max(un -dt * mu,0)*u[i]/un
# le calcul est pour h[1] a h[n] et pour u[1] a u[n]
# les valeurs en 0 et n+1 sont pour les conditions aux limites
#condition de neumann en 0 et en n
    u[0]=u[1]
    h[0]=h[1]
    h[n+1]=h[n]
    u[n+1]=u[n]


#
#
# definition de la solution exacte de la rupture de barrage
# pour un fluide parfait
def sol_exact(x,t):
    h_exact=(((x-0)<=-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<=2*t))
    return h_exact
#
#test au temps 5
sol_num(5)

for i in range(0,n+2):
   print(i, x[i],h2c[i])


print("trace des tas aux temps 0 1 2 3 4 5 ")
print("on constate que pour t>3, ils sont superposes: arret du tas ")

plt.figure(figsize=(8,6))
plt.plot(x,h0c,'b',linestyle='-',label='t=0')
plt.plot(x,h1c,'b',linestyle='-',label='t=1')
plt.plot(x,h2c,'b',linestyle='-',label='t=2')
plt.plot(x,h3c,'b',linestyle='-',label='t=3')
plt.plot(x,h4c,'g',linestyle='-',label='t=4')
plt.plot(x,h,'r',linestyle='-',label='t=5')
plt.title('Solution num avalanche ')
plt.xlabel('x')
plt.ylabel('h')
plt.grid()
plt.legend(loc='best')
plt.show()
