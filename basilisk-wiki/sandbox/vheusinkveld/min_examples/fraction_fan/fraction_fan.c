/** We investigate the volume error caused by fraction.h of a sliced sphere (a 'fan').*/
#include "grid/octree.h"
#include "utils.h"
#include "fractions.h"

/** Defining a starting and an ending level, and the global scalar field fan[], since this field is iteratively used to refine where the fan is. */
int startlevel = 1;
int endlevel = 7;
int ilevel;
scalar fan[];

/** For the volume of a sliced sphere one can visit [wikipedia](https://en.wikipedia.org/wiki/Spherical_cap) */

int main() {
	init_grid(2<<startlevel);
	for(ilevel = startlevel; ilevel<=endlevel;ilevel++){
  
	double VolEst = 0.;	// Estimated volume from fractions.h
	double VolCalc = 0.; 	// Calculated volume 

  	double L0 = 1.;
  	X0 = Y0 = Z0 = 0.;
 
  	double R = L0/4.;
  	double w = R/5.;
	
	VolCalc = 4./3.*M_PI*pow(R, 3.) - 2*M_PI*pow(R-w/2., 2.)/3.*(3*R - (R-w/2.));
  
  	double xf = L0/2.;
  	double yf = L0/2.;
  	double zf = L0/2.;

  	double rTheta = 0.;  
  	double rPhi = 0.;
	
    	coord fn = {sin(rTheta)*cos(-rPhi + M_PI/2.), sin(rTheta)*sin(-rPhi + M_PI/2.), cos(rTheta)};

	scalar sph[], planeup[], planedown[];
	fan.prolongation = fraction_refine;
	refine(fan[]>0.000001 && level < ilevel); // Refined where needed, 0.0001 is arbitrarily chosen since otherwise due to truncation errors it will refine the whole domain.

    	fraction(sph, -sq((x - xf)) - sq((y - yf)) - sq((z - zf)) + sq(R));
    	fraction(planeup, fn.x*(x-xf) + fn.y*(y-yf) + fn.z*(z-zf) + w/2.);
    	fraction(planedown, -fn.x*(x-xf) - fn.y*(y-yf) - fn.z*(z-zf) + w/2.);
    	foreach () {
      		fan[] = sph[]*planeup[]*planedown[];
	}
	int nn = 0.;
	int nc = 0.;
	foreach (reduction(+:VolEst) reduction(+:nn) reduction(+:nc)) {
		VolEst += dv()*fan[];
		nc++;
		if(fan[]>0.){
			nn++;
		}
	}
	// printing restults
 	printf("VolCalc=%g\tVolEst=%g\terr=%.2g%%\tnn=%d\tnc=%d\tlvl=%d\tred=%g \n", VolCalc, VolEst, 100.*(VolCalc-VolEst)/VolCalc, nn, nc, ilevel, (double)((2<<(ilevel*3))/nc));
	
/** Also we would like to see how the error varies with fan radius to delta ratio */ 
	if(ilevel == endlevel){

  	double L0 = 1.;
	double startR = L0/4.;
	double endR = L0/100.;
	
	for(R=startR; R>= endR; R-=L0/100.){

	double VolEst = 0.;	// Estimated volume from fractions.h
	double VolCalc = 0.; 	// Calculated volume 

  	X0 = Y0 = Z0 = 0.;
 
  	double w = R/5.;
	
	VolCalc = 4./3.*M_PI*pow(R, 3.) - 2*M_PI*pow(R-w/2., 2.)/3.*(3*R - (R-w/2.));
  
  	double xf = L0/2.;
  	double yf = L0/2.;
  	double zf = L0/2.;

  	double rTheta = 0.;  
  	double rPhi = 0.;
	
    	coord fn = {sin(rTheta)*cos(-rPhi + M_PI/2.), sin(rTheta)*sin(-rPhi + M_PI/2.), cos(rTheta)};

	scalar sph[], planeup[], planedown[];
	fan.prolongation = fraction_refine;

    	fraction(sph, -sq((x - xf)) - sq((y - yf)) - sq((z - zf)) + sq(R));
    	fraction(planeup, fn.x*(x-xf) + fn.y*(y-yf) + fn.z*(z-zf) + w/2.);
    	fraction(planedown, -fn.x*(x-xf) - fn.y*(y-yf) - fn.z*(z-zf) + w/2.);
    	foreach () {
      		fan[] = sph[]*planeup[]*planedown[];
	}
	int nn = 0.;
	int nc = 0.;
	foreach (reduction(+:VolEst) reduction(+:nn) reduction(+:nc)) {
		VolEst += dv()*fan[];
		nc++;
		if(fan[]>0.){
			nn++;
		}
	}
	// printing restults
 	printf("VolCalc=%g\tVolEst=%g\terr=%.2g%%\tnn=%d\tlvl=%d\tratio=%g\n", VolCalc, VolEst, 100.*(VolCalc-VolEst)/VolCalc, nn, ilevel, (double)(R/(L0/(2<<ilevel))));
	}
	}
	}	
}

/**## output resolution changes
VolCalc=0.00978475      VolEst=0.0046875        err=52% nn=3    nc=64   lvl=1   red=0

VolCalc=0.00978475      VolEst=0.0046875        err=52% nn=3    nc=64   lvl=2   red=2

VolCalc=0.00978475      VolEst=0.00288628       err=71% nn=12   nc=85   lvl=3   red=12

VolCalc=0.00978475      VolEst=0.00879942       err=10% nn=100  nc=2031 lvl=4   red=4

VolCalc=0.00978475      VolEst=0.00969119       err=0.96%       nn=448  nc=11558        lvl=5   red=5

VolCalc=0.00978475      VolEst=0.00975704       err=0.28%       nn=3424 nc=84120        lvl=6   red=6

VolCalc=0.00978475      VolEst=0.00977837       err=0.065%      nn=26592        nc=384063       lvl=7   red=10

## output change in radius
VolCalc=0.00978475      VolEst=0.00977837       err=0.065%      nn=26592        lvl=7   ratio=64

VolCalc=0.00865692      VolEst=0.00865155       err=0.062%      nn=24640        lvl=7   ratio=61.44

VolCalc=0.00761927      VolEst=0.00761443       err=0.064%      nn=17048        lvl=7   ratio=58.88

VolCalc=0.00666803      VolEst=0.00666276       err=0.079%      nn=15552        lvl=7   ratio=56.32

VolCalc=0.00579946      VolEst=0.00579374       err=0.099%      nn=14128        lvl=7   ratio=53.76

VolCalc=0.00500979      VolEst=0.00500403       err=0.12%       nn=12968        lvl=7   ratio=51.2
  
VolCalc=0.00429527      VolEst=0.00428985       err=0.13%       nn=11760        lvl=7   ratio=48.64
  
VolCalc=0.00365214      VolEst=0.00364722       err=0.13%       nn=10504        lvl=7   ratio=46.08
  
VolCalc=0.00307664      VolEst=0.00307226       err=0.14%       nn=9408 lvl=7   ratio=43.52
  
VolCalc=0.00256501      VolEst=0.0025616        err=0.13%       nn=8360 lvl=7   ratio=40.96

VolCalc=0.00211351      VolEst=0.00211035       err=0.15%       nn=4928 lvl=7   ratio=38.4
  
VolCalc=0.00171836      VolEst=0.00171487       err=0.2%        nn=4304 lvl=7   ratio=35.84
  
VolCalc=0.00137581      VolEst=0.00137222       err=0.26%       nn=3744 lvl=7   ratio=33.28
  
VolCalc=0.00108212      VolEst=0.00107871       err=0.32%       nn=3248 lvl=7   ratio=30.72

VolCalc=0.000833504     VolEst=0.000830407      err=0.37%       nn=2720 lvl=7   ratio=28.16
  
VolCalc=0.000626224     VolEst=0.000623392      err=0.45%       nn=2240 lvl=7   ratio=25.6
  
VolCalc=0.000456517     VolEst=0.000454223      err=0.5%        nn=1872 lvl=7   ratio=23.04
  
VolCalc=0.000320627     VolEst=0.000319048      err=0.49%       nn=1488 lvl=7   ratio=20.48
  
VolCalc=0.000214795     VolEst=0.000213283      err=0.7%        nn=568  lvl=7   ratio=17.92
  
VolCalc=0.000135264     VolEst=0.000133753      err=1.1%        nn=432  lvl=7   ratio=15.36
  
VolCalc=7.8278e-05      VolEst=7.69493e-05      err=1.7%        nn=312  lvl=7   ratio=12.8
  
VolCalc=4.00783e-05     VolEst=3.89817e-05      err=2.7%        nn=224  lvl=7   ratio=10.24
  
VolCalc=1.69081e-05     VolEst=1.59719e-05      err=5.5%        nn=120  lvl=7   ratio=7.68

VolCalc=5.00979e-06     VolEst=4.38076e-06      err=13% nn=64   lvl=7   ratio=5.12

## Discussion
It can be seen that the error quickly decreases with increasing resolution. The results shown here is for an grid alligned fan, these numbers are subjec to change when it fan is placed under an angle.   
This is very rough data atm. 
*/
