/**
# Basic advection scheme
# Corrected advection scheme

![Basic advection scheme](advec_circle4/basic.mp4)
![Corrected advection scheme](advec_circle4/corrected.mp4)


 */


//#define Corrected 1

#include "advection.h"
#if Corrected
#include "headers/vof4.h"
#else
#include "vof.h"
#endif
#include "headers/norm-4p.h"
#include "view.h"


scalar f2[],f3[],f1[],f4[];
int nb_phases=4;
scalar * interfaces = {f1,f2,f3,f4}, * tracers = NULL;

#define theta0 0
#define theta1 60
#define theta2 60
#define dtheta 10

#define R 0.3
#define cercle(x,y,R) (sq(R) - sq(x) - sq(y))
#define uppery(x,y) y-x
#define lowery(x,y) x-y
//#define fluid1(x,y) max(x*tan(45*pi/180)-y,-cercle(x,y,R)) //intersection(cercle(x,y-0.001,R),lowery(x,y))
//#define fluid2(x,y) min(x*tan(45*pi/180)-y,cercle(x,y,R)) //intersection(cercle(x,y-0.001,R),lowery(x,y))



#define fluid1(x,y) max(max(x*tan((theta0+dtheta)*pi/180),x*tan((theta1+dtheta)*pi/180)) - y,-cercle(x,y,R))
#define fluid2(x,y) min(min(x*tan((theta0+dtheta)*pi/180),-x*tan((theta2-dtheta)*pi/180)) - y,cercle(x,y,R))
#define fluid3(x,y) max(min(x*tan((theta0+dtheta)*pi/180),-x*tan((theta2-dtheta)*pi/180)) - y,-cercle(x,y,R))

int level=5;
double dtuser=0.1;
double tmax=6.;
#define grid_size level

scalar normU[];

//#define LISSAGE
//#define FILTERED
#define R_VOFLIMIT 1e-3
#define maxlevel 7
#define minlevel level

stats inits1;
stats inits2;
stats inits3;
stats inits4;

int main(int argc, char * argv[]){
	origin(-0.5,-0.5);
	size(1);
	init_grid(pow(2,level));
	run();
}

scalar champ[],floc_err[],f1e[],f2e[],f3e[],f4e[];
char path[80],tempo[160],comm1[160];

event init(i=0) {

	scalar f5[],f6[];

	fraction(f5,fluid1(x,y));
	fraction(f2,fluid2(x,y));
	fraction(f6,fluid3(x,y));

	foreach(){
		f1[]=1-f5[];
		f3[]=f6[]-f2[];
		f4[]=clamp(1.-f1[]-f2[]-f3[],0.,1.);
		f1e[] = f1[];
		f2e[] = f2[];
		f3e[] = f3[];
		f4e[] = f4[];
   	 	//rhov[]=rho(f2[],f1[]);
	}
	boundary ((scalar *){f1,f2,f3,f4});
/*
	FILE * fp4 = fopen("volume_initial","w");
	double v1=0.;
	foreach(){
		v1+=dv();
		fprintf(fp4,"%f\t,%f\t,%f\t,%f\t,%f\t,%f\t,%f\t,%f\n",x,y,(f1[]+f2[]+f3[]+f4[])-1,v1,f1[],f2[],f3[],f4[]);

	}
	fclose(fp4);

	FILE * fp1 = fopen("fraction_error","w");
    //fprintf(fp1, "Temps\t, x\t , y\t, fraction_error\n,f\n",grid_size );
	fclose(fp1);

	FILE * fp3 = fopen("norm2","w");
	fprintf(fp3, "%d\n",level);
	fclose(fp3);
	inits1=statsf(f1);
	inits2=statsf(f2);
	inits3=statsf(f3);
	inits4=statsf(f4);
*/
}


event velocity (i++){
	vertex scalar psi[];
	foreach_vertex(){
		psi[] = - 1.5*sin(2.*pi*t/tmax)*sin((x-0.5)*pi)*sin((y-0.5)*pi)/pi;
	}
	trash ({u});
	struct { double x, y; } f1 = {-1.,1.};
	foreach_face(){
		u.x[] = f1.x*(psi[0,1] - psi[])/Delta;
	}	
}

event movie(t=0. ;t+=dtuser ; t<=tmax) {
	clear();
	view (fov=20 ,tx = -0, ty =-0);
	draw_vof ("f1",filled=1, fc = {0.7,0.,0.}, lw = 2);
	draw_vof ("f3",filled=1, fc = {0.7,0.7,0.} ,lw = 2);
	draw_vof ("f2",filled=1, fc = {0.,0.7,0.7},lw = 2);
	draw_vof ("f4",filled=1, fc = {0.7,0.,0.7},lw = 2);
	cells();
	#if Corrected
	save("corrected.mp4");
	#else
	save("basic.mp4");
	#endif

/*
	if (t==tmax){
		FILE * fp5 = fopen("volume_final","w");
		double v1=0.;

		foreach(){
			v1+=dv();
			fprintf(fp5,"%f\t,%f\t,%f\t,%f\t,%f\t,%f\t,%f\t,%f\n",x,y,(f1[]+f2[]+f3[]+f4[])-1,v1,f1[],f2[],f3[],f4[]);

		}
		fclose(fp5);
	}
*/
}
