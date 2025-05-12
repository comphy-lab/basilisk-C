"""
@author: Alexandre

Ce code permet d'avoir l'evolution du rim et le lambda max de la perturbation.

Code qui lit tous les fichiers crown* puis trie les pas de temps t/tau qui sont dans les noms des fichiers.
Les fichiers crown* contiennent les coordonnees (x,y,z=0) du bords du rim.
Je trace le rayon en fonction de theta. Puis je fais une PSD pour avoir l'evolution du lambda max.
"""

import numpy as np
from matplotlib import pyplot as plt
import glob
import sys

Oh = sys.argv[1] # Here type(Oh) = string
R0 = float(sys.argv[2])

txtfiles = []
kmax=[]
times=[]
time=[]
test=[]
# find files and add times value
for file in glob.glob("crown*"+str(Oh)):
    if ((file.split("crown")[1]).split(str(Oh))[0] == ""):
    	times.append(str(Oh))
    else :
    	times.append((file.split("crown")[1]).split(str(Oh))[0])
# check if int or float
for x in times:
	if (x.find(".") == -1):
		time.append(int(x))
	else :
		time.append(float(x))
time.sort()

theta=[]
Radius=[]
f2 = plt.figure(3)
ax2 = plt.gca()
ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
ticks_labels = [r'-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$']
plt.xticks(ticks, ticks_labels)
plt.xlabel("theta")
plt.ylabel("R/R(t=0)")

# add name files in list
for step in time:
	txtfiles.append("crown"+str(step)+str(Oh))
	data = np.loadtxt("crown"+str(step)+str(Oh))
	data = data.reshape((len(data),3))

	theta = np.arctan2(data[:,1],data[:,0]) #calcul of theta
	Radius = np.sqrt(data[:,1]**2 + data[:,0]**2)/R0

# trier theta et radius par ordre croissant selon theta
	ind = np.argsort(theta)
	theta = np.sort(theta)
	Radius = Radius[ind]

	if(step%2.==0.):
		plt.plot(theta,Radius,label="t/tau="+str(step))

# Debut calcul PSD
	Mradius = Radius - np.mean(Radius)
	Mradius = Mradius/R0
	n = len(theta) # length
	Radiushat = np.fft.fft(Mradius,n) #fft of the Radius
	PSD = Radiushat * np.conj(Radiushat)/n #Power spectral density
	delta = 0.0001

	k = (1/(delta*n)) * np.arange(n) #wavenumber
	indice = np.where(PSD==np.max(PSD)) #find position indice max PSD array
	kmax.append(k[indice][0])
	#print(k[indice][0]) #k_{max}
	L = np.arange(1,np.floor(n/2), dtype='int')

plt.legend()
f2.tight_layout()
plt.savefig("rim"+str(Oh)+".pdf")

# filtered the signal
"""
indices = PSD > np.max(PSD)/3. 
PDSclean = PSD * indices
Radiushatclean = indices * Radiushat
ffilt = np.fft.ifft(Radiushatclean)"""

print(kmax)

# plots
fig,axs = plt.subplots(2,1)
plt.sca(axs[0])
plt.plot(theta,Radius,label="Crown")
ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
ticks_labels = [r'-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$']
plt.xticks(ticks, ticks_labels)
"""plt.plot(theta,ffilt, label="Filtered")"""
plt.xlabel("theta")
plt.ylabel("R/R(t=0)")
plt.legend()
plt.title(r'k$_{max}$=%.2f' % (k[indice][0]))

plt.sca(axs[1])
plt.plot(k[L],PSD[L]/np.max(PSD))
plt.xlim(0,k[L[-1]])
plt.xlabel("k")
plt.ylabel("PSD")

fig.tight_layout()
#plt.show()

# lambda en fonction de t/tau
f1 = plt.figure(2)
ax = plt.gca()

max_value = np.max(time)
min_value = np.min(time)
number_of_steps = 5
l = np.arange(min_value, max_value+1, number_of_steps)

Lambda=[]
for x in kmax:
	Lambda.append(2*np.pi/x)
ax.set(xticks=l, xticklabels=l)
plt.plot(time,Lambda)
plt.xlabel("t/tau")
plt.ylabel(r"$\lambda_{max}$")
f1.tight_layout()
plt.savefig("lambda"+str(Oh)+".pdf",bbox_inches='tight')
#plt.show()


