
/**
# **Wave**
## Set up for the code 
* Download [bash.sh](http://basilisk.fr/sandbox/aaykin/wave_code/bash.sh) and wave.c
* Run:. */

/**
## Code

### Variables and Libraries used */

/**
# **Wave** 
## Set up for the code 
* Download [bash.sh](http://basilisk.fr/sandbox/aaykin/wave_code/bash.sh) and wave.c
* Run:. */

/**
## Code

### Libraries used */

#define P 0
#define multiple_run_LEVEL 0
#define multiple_run_time 1
#define multiple_run_nl 0
#define multiple_run_KH 0
#define multiple_run_CAM 0
#define mutliple_run_RD 0
#define multiple_run_image 0
#define multiple_run_TOL 0
#define data_recover 1
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/check_eta.h"
#if P
#include "layered/perfs.h"
#endif



/**
### Initialisation of variables, main part of the code and exit condition
*/

//Physical variables
double A,om,H,AMP; //A: Constant link to the Altiude, om: Wave's pulsation, H: Water's depth; A: Wave's Alitude
double k,CAM=0.001,RD=1,KH; //DZ: Vertical space step


//Numerical variables
int POINT; //POINT: Number of point on the grid according x
double DZ,DX; //DZ: Vertical space step, DX: Horizontal space step
int t_end=10,LEVEL=10; //t_end: Final time of the simulation, LEVEL: Linked to the number such as POINT=2^LEVEL; Defines the number of points 

//Intern variables, no physical or numerical meanings
double tb;
double z1=0.0,z2=0.0,err_avg_time=0.0,s1=0,s2=0,powerb=1,nlb=0.01,power=2,tc,zs,tb2,Lz,Dtb,told;
int c=1,c2,c3,c4=0,nm=0;

event init (i = 0){
  foreach(serial){
    double z = zb[];
    foreach_layer(){
      h[] = (H - ((A*k)/om)*sin(k*x)*sinh(k*H))/nl;
      //h[]= z0(x,0)/nl;
      //printf("%g\n",z);
      z += h[]/(2.0);
      u.x[] = -(A*k)*(sin(k*x)*cosh(k*z));
      w[]   = (A*k)*(cos(k*x))*sinh(k*z);
      z += h[]/(2.0);
      //w[]   = (A*k)*(cos(k*x))*sinh(k*z);
    }
  }
}

int main(){
  
  periodic(right);
  CFL_H=0.5;
  //linearised = true;
  G=9.81; 
					  

					  
  //Creation of saving files
  FILE * fp = fopen ("error.txt", "w");
  fprintf(fp,"%s %s %s %s %s %s %s %s %s","x_position","Time","POINT_number","LEVEL","Amplitude","Layer_number","Error","Z0","A");
  fprintf(fp,"\n\n");
  fflush(fp);
  FILE * fp2 = fopen ("surface_x.txt", "w");
  fprintf(fp2,"%s %s %s %s %s %s %s %s %s %s %s","Time","x_position","h","ux","w","POINT_number","LEVEL","H","Alitude","Layer number","Surface");
  fprintf(fp2,"\n\n");
  fflush(fp2);
  FILE * fp8 = fopen ("surface_y.txt", "w");
  fprintf(fp8,"%s %s %s %s %s %s %s %s %s %s %s","Time","x_position","h","ux","w","POINT_number","LEVEL","H","Alitude","Layer number","Surface");
  fprintf(fp8,"\n\n");
  fflush(fp8);
  FILE * fp9 = fopen ("surface.txt", "w");
  fprintf(fp9,"%s %s %s %s %s %s %s %s","Time","x_position","POINT_number","LEVEL","H","Alitude","Layer number","Surface");
  fprintf(fp9,"\n\n");
  fflush(fp9);
  FILE * fp3 = fopen ("phase.txt", "w");
  fprintf(fp3,"%s %s %s %s %s %s %s %s %s %s","x_position","time","Number_point","LEVEL","CAM","KH","H","Layer_number","Surface","Tolerance");
  fprintf(fp3,"\n\n");
  fflush(fp3);
  FILE * fp4 = fopen ("error_avg.txt", "w");
  fprintf(fp4,"%s %s %s %s %s %s %s %s %s","Time","POINT_number","LEVEL","Altitude","H","Layer_number","Camber","Surface","Average_Error");
  fprintf(fp4,"\n\n");
  fflush(fp4);
  FILE * fp5 = fopen ("variable.txt", "w");
  fprintf(fp5,"%s %s %s %s %s %s %s %s %s %s %s","RD","CAM","H","AMP","POINT","L0","nl","DZ","DT","KH","TOL");
  FILE * fp6 = fopen ("eta_time.txt", "w");
  fprintf(fp6,"%s %s %s\n\n","time","eta","x");
  fflush(fp6);
  FILE * fp7 = fopen ("phase_time.txt","w");
#if multiple_run_TOL
  for (TOLERANCE = pow(10,-4); TOLERANCE >= pow(10,-10) ; TOLERANCE *= 0.1)
#else
    TOLERANCE = 1e-8;
#endif
  {

#if multiple_run_nl
    for (nl = 900; nl <= 2000;  nl += 100)
#else
      nl = 40;
#endif
      {
#if multiple_run_CAM
	for (CAM = pow(10,-4)*5 ; CAM <=  pow(10,-2); CAM *= 1.2)
#endif
	  {
#if multiple_run_LEVEL
	    for (LEVEL = 1; LEVEL <= 10; LEVEL += 1)
#endif
	      {
#if multiple_run_KH
		for (KH = 0.1 ; KH <=  100; KH += (KH>=3)?1:0.1)
#endif
		  {
#if multiple_run_RD
		    for (RD = 0.1 ; RD <=  10; RD *= 10 )
#endif
		      {
			POINT = pow(2,LEVEL);
			init_grid(POINT);
			//AMP = 1;
			//L0 = AMP/CAM;
			k = (2*pi)/(10.0);
			L0 = (2*pi)/k; // Domain set as one wave length
			H = RD*L0;
			//H = KH/k;
			AMP = (CAM*L0);
			//om = sqrt( G*k*tanh(k*H) * (1 + ( (9 - 10*pow((tanh(k*H)),2) + 9*pow((tanh(k*H)),4))/(8*pow((tanh(k*H)),4)))*pow((k*AMP),2) ) );
			om = sqrt(G*k*tanh(k*H));
			A = (AMP*om)/(k*sinh(k*H));
			DX = L0/POINT;
			DZ = H/nl;						
			fprintf(fp5,"\n\n");
			fprintf(fp5,"%g %g %g %g %i %g %i %g %.10g %g",RD,CAM,H,A,POINT,L0,nl,DZ,dt,KH,TOLERANCE);
			fflush(fp5);

			c2=0;
			c3=0;
			c4=0;
			run();
		      }
		  }
	      }
	  }
      }
  }
}
/**
### Code's events and functions
*/

/**
* function z0
 */

double z0 (double q,double s){
  return (H - ((A*k)/om)*sin(-s*om+k*q)*sinh(k*H));
}



/**
* output's events
 */
	
#if data_recover

event error (t = t_end){	 
  static FILE * fp = fopen ("error.txt", "a");
  static FILE * fp2 = fopen ("error_avg.txt", "a");
  double err_avg=0.0;

  foreach(serial){ 
    err_avg += (fabs(eta[]-z0(x,t))/AMP);
    fprintf (fp,"%g %g %g %i %i %g %i %g %g %g\n",x,L0,t,POINT,LEVEL,A,nl,fabs(eta[]-z0(x,t)),z0(x,t),((A*k)/om)*sinh(k*H));
    fflush (fp);  
  }

  fprintf (fp2,"%g %i %i %g %i %g %g\n",t,POINT,LEVEL,A,nl,CAM,((err_avg)/POINT));
  fflush(fp2);
}



#if multiple_run_time
event surface (i++)
{
  static FILE * fp = fopen ("surface.txt", "a");

  foreach(serial)
  {
    //foreach_layer()
      {
	fprintf (fp,"%g %g %i %i %g %g %i %.10g\n",t,x,POINT,LEVEL,H,A,nl,eta[]);
	fflush (fp);

      }
  
  }
}

event surface_x (t=t_end){
  double sd = 0;
  static FILE * fp = fopen ("surface_x.txt", "a");

  foreach(serial){
    sd = 0;
    foreach_layer(){
      sd += h[];
      if((x >= 0.499) & (x<= 0.501)){
	fprintf (fp,"%g %g %g %g %g %i %i %g %g %i %.10g\n",t,x,sd,u.x[],w[],POINT,LEVEL,H,A,nl,eta[]);
	fflush (fp);

	}
  
    }
  }
}
event surface_y (t=t_end)
{
  double sd = 0;
  double sdb = 0;
  double xb = 0;
  double c = 0;
  double d = 0;
  static FILE * fp2 = fopen ("surface_y.txt","a");
  foreach(serial)
  {
    sd = zb[];
    sdb = zb[];

    c = 0;
    d = 0;
      
    foreach_layer()
      {
	sd += 0.025;
	sdb += h[]/2;
         
	if((sdb >= 0.45) &(sdb <= 0.55) &  (c == 0))// &*/ (t==t_end))
	  {//& (t>=0.5)){
	    fprintf (fp2,"%g %g %g %g %g %i %i %g %g %i %.10g\n",t,x,sdb,u.x[],w[],POINT,LEVEL,H,A,nl,eta[]);
	    fflush (fp2);
	    	c += 1.0;
	  }
	sdb += h[]/2;
	d += 1;


      }
    //xb = x;
  }
}
#endif

#endif 



#if multiple_run_time
event phase(i++){
  int c=0,o=0; 
      static FILE * fp = fopen ("phase.txt", "a");
      static FILE * fp2 = fopen ("eta_time.txt", "a");
      foreach(serial){ 
	c += 1;
	if(c == 550){
	  if(i == 0){
	    z1 = eta[];
	  }
	}
      }
      c = 0;
      foreach(serial){ 
	c += 1;
	if((x>L0*0.49) & (x < 0.54*L0) & (o==0)){
	  o = 1;
	  fprintf(fp2,"%g %g %g\n",t,eta[],x);
	  if(eta[]-H > 0) s2=1;
	  if(eta[]-H < 0) s2=0;
	  if((s1!=s2) & (s1==0) ){
	    s1=1;
   
	    if (c2==1+nm*2){
	      Lz = zs + fabs(eta[]-H);
	      Dtb = t-told;
	      tc = t - (1 - (zs/Lz))*Dtb - tb;
	      //tc = (t - fabs( (eta[]-H) / (fabs(eta[]-H)+zs) ) * dt - tb);
	    }
	
	    if (c2==1+nm*2) fprintf(fp,"%g %.10g %.10g %i %i %g %g %g %g %i %.10g %g %g\n",x,tc,(t-tb2),POINT,LEVEL,CAM,k*H,k,H,nl,eta[]-H,tb,TOLERANCE);  
	    fflush(fp);

	    c2+=1;
	
	    /* c3+=1; */
	    /* if (c3==2){ */
	    /*   c2=0; */
	    /*   c3=0; */
	    /* } */

	    tb2 = t;
	    /* tb = (t - fabs( (eta[]-H) / (fabs(eta[]-H)+zs) ) * dt); */

	    Lz = zs + fabs(eta[]-H);
	    Dtb = t-told;
	    tb = t - (1 - (zs/Lz))*Dtb;


	  }
	  if((s1!=s2) & (s1==1)){
	    s1=0;
	  }
	  zs=fabs(eta[]-H);
	  told=t;
	}
    
      }  
      /*if(nlb!=k){
	fflush(fp);
	nlb += 0.01;
	}
	if(powerb!=power){
	nlb = 1;
	}  */
    
}
#endif

event phase_time(i++)
{
  static FILE*fp = fopen("phase_time.txt","w");
  double fb=0;
  foreach(serial)
    {
      if (fb <(eta[]-H))
	{
	  fprintf(fp,"%g %g %g %g",t,x,y,eta[]-H);
	}
      fb = eta[]-H;
    }
}

      

#if multiple_run_KH 
event end (t = 1.05*(2.*pi/sqrt(G*k*tanh(k*H)))) {}
#else
event end (t = t_end) {}
#endif



/**
 * Animation of the surface analytically and numerically according to time
 */

#if multiple_run_image
event plot (i++) {
  static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
  if(t == 0) fprintf (fp,"set term gif animate size 1920,1080;set output 'wave.gif';set size ratio 0.5\n");
  fprintf (fp,"\nset grid\n");
  fprintf (fp,"set title 'Wave at t= %.2lf '\n"
	   "G = 9.81 \n"
	   "k = %.2lf \n"
	   "L0 = (2*pi)/k \n"
	   "H = %.2lf \n"
	   "om = sqrt(G*k*tanh(k*H)) \n"
	   //"om = sqrt( G*k*tanh(k*H) * (1 + ( (9 - 10*(tanh(k*H))**2 + 9*(tanh(k*H))**4)/(8*(tanh(k*H)**4)))*(k*%.2lf)**2) ) \n"
	   "A = %.10g \n"	
	   "set xrange[0:L0] \n"
	   "set yrange[H-(A*k/om)*sinh(k*H)*1.1:H+(A*k/om)*sinh(k*H)*1.1] \n"
	   "set samples 10000 \n"
	   "set xlabel 'x' \n"
	   "set ylabel 'z' \n"
	   "z(x) = H - ((A*k)/om)*sin(-t*om+k*x)*sinh(k*H); t= %.2lf ;\n"
	   "p[%g:%g]  '-' u 1:2 t'Numerical result' w lp ,"
	   "z(x) t 'Analytical result'\n",
	   t,k,H,A,t,X0,X0+L0);
  foreach(serial) fprintf (fp,"%g %.15g %g\n",x,eta[],t);
  fprintf (fp,"e\n\n");
}
#endif

