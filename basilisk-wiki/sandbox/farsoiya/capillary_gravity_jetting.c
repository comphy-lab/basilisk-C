
/**
# Nonlinear Capillary-gravity wave
We investigate the nonlinear evolution of capillary-gravity wave
Ref : 

![Wave geometry](./Jetting_capillary_gravity/3dwave.png)

*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "tension.h" 
#include "navier-stokes/perfs.h" 
#include "maxruntime.h"

int resolution =256;
#define epsilon 0.5
#define l_5 16.47063
#define mean 2.5

//BOUNDARY CONDITIONS FOR FACE-CENTERED VELOCITY FIELD//
//NO PENETRATION B.C//
uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

int main(int argc, char * argv[]){
  //resolution = atoi(argv[1]);
	//PARAMETERS//
	G.x = -981.;
	rho1 = 1.0, rho2 = 0.0010;
	f.sigma = 72.0;
	mu1 = 0.00, mu2 = 0.00;
	TOLERANCE = 1e-6;
	
	size(0.26*l_5);
	init_grid(resolution);
	
	
		
	run();
}

//INITIAL COLOR FUNCTION FIELD//
event init (t = 0){
	double domain = 0.26*l_5;
	double a0 = epsilon*domain/l_5;
	//~ printf("a0 = %g\n",a0);
	fraction (f, -((x - mean) - a0*( j0((l_5/domain)*y) )));
}

event vof (i++, first);

double pos_t_0=0;
event interface (t = 0; t <= 0.1; t+=0.0001){
	
	scalar pos[];
	position (f, pos, {1,0});
	//~ double max = statsf(pos).max;
	//~ double max = interpolate (pos, x, y);
	double y_min=100., pos_r_0=0;
	
	foreach(){
		if(pos[]!=nodata){
			if (y_min>y){
				y_min = y;
				pos_r_0 = pos[]-mean;
			}
		}
		
	}
	if (i==0){
		pos_t_0 = pos_r_0;
	}

	char name[80];
	sprintf (name, "wave-%d", N);
	static FILE * fp = fopen (name, "w");
	fprintf (fp, "%g %g\n", t, pos_r_0/pos_t_0);
	fflush (fp);	

}
/**
# Results

 ~~~pythonplot Amplitude vs time
	import numpy as np
	import matplotlib
	import matplotlib.pyplot as plt
	matplotlib.use('Agg')
	#from matplotlib import rc
        # Analytical data can be reproduced using this script on local machine where scipy package is available
	#import scipy.special as sp 

	# SURFACE TENSION
	T = 72.

	# GRAVITY
	g = 981.

	# DENSITY
	rho = 1.

	
	df = np.loadtxt("../J1_roots.csv")
	
	eps = 0.5
	q = 5

	l_q = df[q - 1]
	domain = 0.26*l_q

	# PRIMARY WAVENUMBER
	k = l_q/domain

	# INVERSE BOND NUMBER
	sigma = T*(k**2)/(rho*g)
	
	# NUMERICAL INTEGRATION TRAPEZOIDAL RULE
	
	#~ def integrateTrap(r,alpha_j,toggle):
		#~ if (toggle == 'Type1'):
			#~ func = np.array([num*(sp.jv(0,num)**2)*sp.jv(0,alpha_j*num) for num in r])
		#~ elif (toggle == 'Type2'):
			#~ func = np.array([num*(sp.jv(1,num)**2)*sp.jv(0,alpha_j*num) for num in r])
		#~ elif (toggle == 'Type3'):
			#~ func = np.array([num*sp.jv(0,num)*sp.jv(1,num)*sp.jv(1,alpha_j*num) for num in r])
		#~ elif (toggle == 'Type4'):
			#~ func = np.array([num*(sp.jv(0,num)**4) for num in r])
		#~ elif (toggle == 'Type5'):
			#~ func = np.array([num*(sp.jv(0,num)**2)*(sp.jv(1,num)**2) for num in r])
		#~ elif (toggle == 'Type6'):
			#~ func = np.array([num*sp.jv(0,num)*(sp.jv(1,num)**2)*sp.jv(2,num) for num in r])
		#~ else:
			#~ print('ERROR IN integrateTrap ~ SELECT CORRECT TYPE !!!')
			#~ sys.exit()
    
		#~ delta = r[1] - r[0]
		#~ return 0.5*delta*(2*np.sum(func) - func[0] - func[len(func) - 1])
	
	#~ # $\textrm{Defining}$ $\eta = \epsilon \eta_1(r,t) + \epsilon^2 \eta_2(r,t)$
	
	#~ def etaNonlinear(time,omega,Omega_2,zeta,g=g,k=k,eps=eps,q=q):
		#~ tau = np.sqrt(g*k)*( 1 + (eps**2)*Omega_2 )*time
		#~ amp = 0.
		#~ for j in range(len(omega)):
			#~ amp += (eps**2)*0.5*(zeta[j][0]*np.cos(omega[j]*tau) + zeta[j][1]*np.cos(2*omega[q - 1]*tau) + zeta[j][2]) 
		#~ amp += eps*np.cos(omega[q - 1]*tau)
		#~ #... normalized ...#
		#~ return amp/eps

	#~ def etaLinear(time,omega,g=g,k=k,q=q):
		#~ #... normalized ...#
		#~ return np.cos(omega[q - 1]*np.sqrt(g*k)*time)
	
	#~ # $\textrm{Computing}$ $\left\{ \zeta_1^{(j)}, \zeta_2^{(j)}, \zeta_3^{(j)} \right\}, \left\{ \xi_1^{(j)}, \xi_2^{(j)} \right\}, \Omega_2$
	
	#~ N = len(df)
	#~ r = np.linspace(0,l_q,5000)

	#~ alpha = np.array([val/l_q for val in df])
	#~ omega = np.sqrt(alpha*(1 + (alpha**2)*sigma))
	#~ I_1 = np.array(list(map(lambda val:integrateTrap(r,val,'Type1'),alpha)))
	#~ I_2 = np.array(list(map(lambda val:integrateTrap(r,val,'Type2'),alpha)))
	#~ I_3 = np.array(list(map(lambda val:integrateTrap(r,val,'Type3'),alpha)))
	#~ I_4 = integrateTrap(r,-10000000,'Type4')
	#~ I_5 = integrateTrap(r,-10000000,'Type5')
	#~ I_6 = integrateTrap(r,-10000000,'Type6')

	#~ zeta = np.zeros((N,3))
	#~ xi = np.zeros((N,2))
	#~ Omega_2 = 0.
	#~ for i in range(N):
		#~ zeta[i][0] = ( (alpha[i] - 2) + (alpha[i]**3 - alpha[i]**2 - 1)*sigma )*I_1[i] + ( 2 + (1 + alpha[i]**2)*sigma )*I_2[i]
		#~ zeta[i][0] = -( 4*(1 + sigma)/(1 + (alpha[i]**2)*sigma) )*zeta[i][0]
		#~ zeta[i][0] = zeta[i][0]/( omega[i]**2 - 4*(omega[q - 1]**2) )
		#~ zeta[i][0] = zeta[i][0]/( (l_q**2)*(sp.jv(0,df[i])**2) )

		#~ zeta[i][1] = ( (3*alpha[i] - 4) + (3*(alpha[i]**3) - 4*(alpha[i]**2))*sigma )*I_1[i] + ( (alpha[i] + 4) + (alpha[i]**3 + 4*(alpha[i]**2))*sigma )*I_2[i]
		#~ zeta[i][1] = ( (1 + sigma)/(1 + (alpha[i]**2)*sigma) )*zeta[i][1]
		#~ zeta[i][1] = zeta[i][1]/( omega[i]**2 - 4*(omega[q - 1]**2) )
		#~ zeta[i][1] = zeta[i][1]/( (l_q**2)*(sp.jv(0,df[i])**2) )

		#~ zeta[i][2] = I_1[i] - I_2[i]
		#~ zeta[i][2] = ( (1 + sigma)/(1 + (alpha[i]**2)*sigma) )*zeta[i][2]
		#~ zeta[i][2] = zeta[i][2]/( (l_q**2)*(sp.jv(0,df[i])**2) )

		#~ xi[i][0] = ( (alpha[i] - 2) + (alpha[i]**3 - alpha[i]**2 - 1)*sigma )*I_1[i] + ( 2 + (1 + alpha[i]**2)*sigma )*I_2[i]
		#~ xi[i][0] = ( 2*omega[q - 1]/omega[i] )*xi[i][0]
		#~ xi[i][0] = xi[i][0]/( omega[i]**2 - 4*(omega[q - 1]**2) )
		#~ xi[i][0] = xi[i][0]/( (l_q**2)*(sp.jv(0,df[i])**2) )

		#~ xi[i][1] = ( 2 + (3 - alpha[i]**2)*sigma )*I_1[i] + ( 2 + (1 + alpha[i]**2)*sigma )*I_2[i]
		#~ xi[i][1] = -xi[i][1]/( omega[i]**2 - 4*(omega[q - 1]**2) )
		#~ xi[i][1] = xi[i][1]/( (l_q**2)*(sp.jv(0,df[i])**2) )

		#~ #... NONLINEAR FREQUENCY CORRECTION ...#
		#~ Omega_2 += ( ((2 + sigma)*alpha[i] - alpha[i]**2)*xi[i][1] + sigma*zeta[i][1] )*I_1[i] 
		#~ Omega_2 += -( 2*(1 + sigma)*alpha[i]*xi[i][1] + 0.5*(alpha[i]**3)*sigma*zeta[i][1] - (alpha[i]**3)*sigma*zeta[i][2] )*I_3[i]
	    
	#~ Omega_2 += - 0.5*(1 + sigma)*I_4 + ( (1 - 4*sigma + sigma**2)/(4*(1 + sigma)) )*I_5 + 0.75*( (1 + 3*sigma + sigma**2)/(1 + sigma) )*I_6
	#~ Omega_2 = 0.5*Omega_2/( (l_q**2)*(sp.jv(0,df[q - 1])**2) )
	#print('Omega_2 = ',Omega_2)

	#~ t = np.arange(0,0.1,0.001)
	#signal_WNLS = etaNonlinear(t,omega,Omega_2,zeta)
	#signal_linear = etaLinear(t,omega)

	
	#rc('font', family='serif')
	#rc('text',usetex=True)
	filename = 'wave-256'
	ts,a = np.loadtxt(filename,delimiter=' ',unpack=True)
	plt.plot(ts[::10],a[::10],'ko');
	#plt.plot(t,etaNonlinear(t,omega,Omega_2,zeta))
	#plt.plot(t,etaLinear(t,omega))
	t,etaN,etaL = np.loadtxt('../analytical.csv',delimiter=' ',unpack=True)
	plt.plot(t,etaN)
	plt.plot(t,etaL)
	plt.legend(['Simulation','upto O($\epsilon^2$)','Linear'])
	
	plt.xlabel('$\hat{t}$',fontsize=7)
	plt.ylabel('$\eta(0,t)/\epsilon$',fontsize=7,rotation=0)
	plt.savefig('plot.png')
	
   ~~~


*/
