# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
from scipy.optimize import curve_fit

Oh = sys.argv[1] # Here type(Oh) = string
H0= float(sys.argv[2])

times=[]
# find files and add times value
for file in glob.glob("interface*"+str(Oh)):
    if ((file.split("interface")[1]).split(str(Oh))[0] == ""):
    	times.append(str(Oh))
    else :
    	times.append((file.split("interface")[1]).split(str(Oh))[0])
# check if int or float
time=[]
for x in times:
	if (x.find(".") == -1):
		time.append(int(x))
	else :
		time.append(float(x))
time.sort()
kmax=[]
lambda_th=[]
f2 = plt.figure(2)
ax = plt.gca()
for step in time:
	data = np.loadtxt("interface"+str(step)+str(Oh))
	data = data.reshape((len(data),3))

	theta = np.round(np.arctan2(data[:,0],data[:,1]),2)

	Zmax=[]
	angle=[]
	for i in theta:
		if(angle.count(i)==0):
			indices = np.where(theta == i)
			Zmax.append(np.max(data[indices[:],2]))
			angle.append(i)


# calcul lambda_th
	lambda_th.append((np.max(Zmax)*4*np.pi)/np.sqrt(2))

# trier angle et Zmax par ordre croissant selon angle
	angle=np.array(angle)
	ind = np.argsort(angle)
	angle = np.sort(angle)
	Zmax=np.array(Zmax)
	Zmax = Zmax[ind]
	
# application filtre intervalle
	ZZ=[]
	ran=[]
	for j in range(0,len(angle)-2,2):
		ran.append(angle[j])
		ZZ.append(np.max(Zmax[j:j+2]))
	if(step%2.==0.):
		plt.plot(ran,ZZ,label="t/tau="+str(step))

# Debut calcul PSD
	Zmax=[]
	Zmax=np.array(ZZ)
	angle=[]
	angle=np.array(ran)

	MZmax = Zmax - np.mean(Zmax)
	MZmax = MZmax/H0
	n = len(angle) # length
	Zmaxhat = np.fft.fft(MZmax,n) #fft of the Radius
	PSD = Zmaxhat * np.conj(Zmaxhat)/n #Power spectral density
	delta = 0.0001

	k = (1/(delta*n)) * np.arange(n) #wavenumber
	indice = np.where(PSD==np.max(PSD)) #find position indice max PSD array
	kmax.append(k[indice][0])
	#print(k[indice][0]) #k_{max}
	L = np.arange(1,np.floor(n/2), dtype='int')
ticks = [0, np.pi/4.,np.pi/2.]
ticks_labels = [r'0', '$\pi/4$','$\pi/2$']
plt.xticks(ticks, ticks_labels)
plt.xlabel("Theta")
plt.ylabel("Z")
plt.legend()
f2.tight_layout()
plt.savefig("H(t)"+str(Oh)+".pdf",bbox_inches='tight')
#plt.show()



# filtered the signal

indices = PSD > np.max(PSD)/10.
PDSclean = PSD * indices
Zmaxhatclean = indices * Zmaxhat
ffilt = np.fft.ifft(Zmaxhatclean)

fig,axs = plt.subplots(2,1)
plt.sca(axs[0])

ticks = [0, np.pi/4.,np.pi/2.]
ticks_labels = [r'0', '$\pi/4$','$\pi/2$']
plt.xticks(ticks, ticks_labels)

plt.plot(angle,ffilt+np.mean(Zmax), label="Filtered")
plt.plot(angle,Zmax, label="Thickness")
plt.xlabel("Angle")
plt.ylabel(r"Z_{max}")

plt.sca(axs[1])
plt.plot(k[L],PSD[L]/np.max(PSD))
plt.xlim(0,k[L[-1]])
plt.xlabel("k")
plt.ylabel("PSD")

fig.tight_layout()
plt.savefig("TF_H"+str(Oh)+".pdf",bbox_inches='tight')
#plt.show()

# lambda en fonction de t/tau
f1 = plt.figure(1)
ax = plt.gca()

max_value = np.max(time)
min_value = np.min(time)
number_of_steps = 5
l = np.arange(min_value, max_value+1, number_of_steps)

# Calcul de lambda_max grace PSD_max
Lambda=[]
for x in kmax:
	Lambda.append(2*np.pi/x)

# fit function
def func(x,a):
	return a*x

# intervalle fit
"""
b1 = time.index(1)
b2 = time.index(2)
popt, pcov = curve_fit(func,time[b1:b2],Lambda[b1:b2])
print(popt)#coeff a"""

# plot lambda_max PSD
ax.set(xticks=l, xticklabels=l)
plt.plot(time,Lambda, label=r"$\lambda_{max}$")
#plt.plot(time, lambda_th, label=r"λ$_{RP}$")
plt.xlabel("t/tau")
plt.ylabel(r"$\lambda_{max}$")
f1.tight_layout()
plt.legend()
plt.savefig("lambda_H"+str(Oh)+".pdf",bbox_inches='tight')
#plt.show()
#print(lambda_th)



