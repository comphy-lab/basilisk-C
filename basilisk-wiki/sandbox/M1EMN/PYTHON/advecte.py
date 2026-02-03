# 02/2026
#
# Resolution de l'equation d'advection simple
# http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf
# http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf
# https://basilisk.fr/sandbox/M1EMN/Exemples/flood.c
# voir Whitham
#
# trucs de python
import numpy as np
import scipy.sparse as sp
import math
import matplotlib.pylab as plt
#parametres
n=500
L=16
dx=L/n
dt=dx*.25
# note: verifier que si CFL n'est pas 1, alors il faut augmenter le nombre de points

#tableau des x
x=np.zeros(n+2)
for i in range(n+2):
  x[i]=(i-n/2)*dx
# hauteur initiale: une exponentielle  
h=np.zeros(n+2)
for i in range(n+2):
  h[i]=np.exp(-x[i]*x[i])
# definition du flux numerique pour la hauteur
# cellule "gauche" et "droite"
def FR1(hg, hd):
  c=1.
  return (hg + hd)*0.5-c*(hd-hg)*0.5
#

# tableaux pour les plots
h0c=np.zeros(n+2)
h1c=np.zeros(n+2)
h2c=np.zeros(n+2)
h3c=np.zeros(n+2)
h4c=np.zeros(n+2)
h5c=np.zeros(n+2)

# sauvegardes
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
      fp[i]=FR1(h[i-1],h[i])
# avancee pour h, si h>0 avancee pour u
    for i in range(1,n):
      h[i]=h[i]-dt*(fp[i+1]-fp[i])/dx;
   
# les valeurs en 0 et n+1 sont pour les conditions aux limites
#condition de neumann en 0 et en n
    h[0]=h[1]
    h[n+1]=h[n]

#
#
#test au temps 5
sol_num(5)

# Solution exacte
x3e=np.zeros(n+2)
h3e=np.zeros(n+2)
for i in range(0,n+2):
  t=3
  h3e[i]=np.exp(-(x[i]-t)**2)


# plots et traces
# en rouge solution exacte
plt.figure(figsize=(8,6))
plt.plot(x,h0c,'b',linestyle='-',label='t=0')
plt.plot(x,h1c,'b',linestyle='-',label='t=1')
plt.plot(x,h2c,'b',linestyle='-',label='t=2')
plt.plot(x,h3c,'b',linestyle='-',label='t=3')
plt.plot(x,h4c,'b',linestyle='-',label='t=4')
plt.plot(x,h5c,'b',linestyle='-',label='t=5')
plt.plot(x,h3e,'r',label='exact t=3')
plt.title('Solution num advecte ')
plt.xlabel('x')
plt.ylabel('h')
plt.grid()
plt.legend(loc='upper left')
plt.show()
