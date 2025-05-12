/**
# **Wave**
## Set up for the code 
* Download [bash.sh](http://basilisk.fr/sandbox/aaykin/wave_code/bash.sh) and wave_3D.c
* Run:. */

/**
## Code

### Libraries used */


#define P 1
#define multiple_run_LEVEL 0
#define multiple_run_time 0
#define multiple_run_nl 0
#define multiple_run_KH 0
#define multiple_run_CAM 0
#define multiple_run_image 0
#define multiple_run_TOL 0
#define multiple_waves 1
#define data_recover 1
#define output_image 1
#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/check_eta.h"
#if P
#include "layered/perfs.h"
#endif

/*
### Initialisation of variables, main part of the code and exit condition
*/

//Physical variables
double A,A1,A2,om,om1,om2,H,AMP,AMP1,AMP2;  //A: Constant link to the Altiude, om: Wave's pulsation, H: Water's depth; AMP: Wave's Alitude
double k,k1,k2,CAM=0.009,RD=1,RL=3.0,KH; //k: Wave number, CAM: Steepness of the wave, RD: Relative depth, RL: Relative lenght

//Numerical variables
int POINT,NB;                                  //POINT: Number of point on the grid according x
double DZ,DX;                               //DZ: Vertical space step, DX: Horizontal space step
int t_end=200,LEVEL=8;                      //t_end: Final time of the simulation, LEVEL: Linked to the number such as POINT=2^LEVEL; Defines the number of points 

//Intern variables, no physical or numerical meanings
double tb;
double z1=0.0,z2=0.0,err_avg_time=0.0,s1=0,s2=0,powerb=1,nlb=0.01,power=2,tc,zs,tb2,Lz,Dtb,told;
int c2,c3,c4=0,nm=0;
scalar maxa[];


//Initialisation of the code
event init (i = 0){
  foreach(serial){
    maxa[]=0;
    double z = zb[];
    foreach_layer(){
#if multiple_waves
      h[] = (H - (((A2*k2)/om2)*sin(k2*x)*sinh(k2*H)) -  (((A1*k1)/om1)*sin(k1*y)*sinh(k1*H)))/nl;
      z += h[]/(2.0);
      u.x[] = -((A2*k2)*sin(k2*x)*cosh(k2*z));
      u.y[] =  -((A1*k1)*sin(k1*y)*sinh(k1*z));
      w[]   = (A1*k1)*(cos(k1*y))*sinh(k1*z) + (A2*k2)*(cos(k2*x))*sinh(k2*z);  
      z += h[]/(2.0);
#else
      h[] = (H - (((A1*k1)/om1)*sin(k1*y)*sinh(k1*H)))/nl;
      z += h[]/(2.0);
      u.x[] = -(A1*k1)*(sin(k1*x)*cosh(k1*z));
      w[]   = (A1*k1)*(cos(k1*x))*sinh(k1*z);  
      z += h[]/(2.0);
#endif
    }
  }
}


//Main part of the code
int main()
{

  //Creation of saving files
  FILE * fp = fopen ("error.txt", "w");
  fprintf(fp,"%s %s %s %s %s %s %s %s %s","x_position","Time","POINT_number","LEVEL","Amplitude","Layer_number","Error","Z0","A");
  fprintf(fp,"\n\n");
  fflush(fp);
  FILE * fp3 = fopen ("phase.txt", "w");
  FILE * fp2 = fopen ("surface_2.txt", "w");
  fprintf(fp3,"%s %s %s %s %s %s %s %s %s","x_position","y_position","time","Number_point","LEVEL","KH","H","Layer_number","Surface");
  fprintf(fp3,"\n\n");
  fflush(fp3);
  FILE * fp4 = fopen ("error_avg.txt", "w");
  fprintf(fp4,"%s %s %s %s %s %s %s %s %s","Time","POINT_number","LEVEL","Altitude","H","Layer_number","Camber","Surface","Average_Error");
  fprintf(fp4,"\n\n");
  fflush(fp4);
  FILE * fp5 = fopen ("variable.txt", "w");
  fprintf(fp5,"%s %s %s %s %s %s %s %s %s %s","RD","CAM","H","AMP","POINT","L0","nl","DZ","DT","KH");
  FILE * fp6 = fopen ("third_wave_formation_time.txt","w");
  FILE * fp7 = fopen ("fourier_t0.txt", "w");
  //Boundary conditions
  periodic(right);
  periodic(top);

  //CFL and tolerance
#if multiple_run_TOL
  for (TOLERANCE=pow(10,-1) ; TOLERANCE >=pow(10,-10); TOLERANCE *= 0.1)
#else
    TOLERANCE = 1e-7;
#endif
  {
    CFL_H=0.5;

    linearised = true;

    //Multiple run loops
#if multiple_run_nl
    for (nl = 1; nl <= 5; nl += 1)
#else
      nl=3;
#endif
    {
#if multiple_run_CAM
      for (CAM = pow(10,-3) ; CAM <=  pow(10,-2); CAM *= 10)
#endif
	{
#if multiple_run_LEVEL
	  for (LEVEL = 1; LEVEL <= 14; LEVEL += 1)
#endif
	    {
#if multiple_run_KH
	      for (KH = 0.1 ; KH <=  100; KH += (KH>=3)?1:0.1)
#endif
		{

		  //Initialisation of the grid
		  POINT = pow(2,LEVEL);
		  init_grid(POINT);
		  Y0    = 0;
		  X0    = 0;
		  NB     = POINT*POINT;

		  //Settings of the physical and numerical parameters
		  G     = 9.81;
		  k2    =  2*pi;
		  k1    = (k2*RL);
		  L0    = (2*pi)/k2;
#if multiple_run_KH
		  RD    = KH/(2*pi);
#endif
		  H     = RD*L0;
		  om    = sqrt(G*k*tanh(k*H));
		  om1   = sqrt(G*k1*tanh(k1*H));
		  om2   = sqrt(G*k2*tanh(k2*H));
		  AMP   = CAM*L0;
		  AMP1  = 0.0057*L0;
		  AMP2  = CAM*L0;
		  // AMP2 = 0;
		  A     = (AMP*om)/(k*sinh(k*H));
		  A1    = (AMP1*om1)/(k1*sinh(k1*H));
		  A2    = (AMP2*om2)/(k2*sinh(k2*H));
		  DX    = L0/POINT;
		  DZ    = H/nl;

		  //Printing of the parameters
		  fprintf(fp5,"\n\n");
		  fprintf(fp5,"%g %g %g %g %i %g %i %g %.10g %g",RD,CAM,H,A,POINT,L0,nl,DZ,dt,KH);
		  fflush(fp5);
		  printf("KH = %g",KH);

		  //Intern variables of the code
		  c2    = 0;
		  c3    = 0;
		  c4    = 0;

		  //Run of the code
		  run();
		}
	    }
	}
    }
  }
}


event maximum (i++)
{
  foreach()
    if (fabs(eta[]) > maxa[])
      maxa[] = fabs(eta[]);
}


//Save of the eta fields
#if output_image
/*event output (i++)
  {
  char filename[64];  
  FILE *fp;
  sprintf(filename, "end_%.2lf", t);
  fp = fopen( filename, "w");
  FILE * fp2 = fopen ("section2", "w");
  FILE * fp5 = fopen ("section5", "w");
  FILE * fp7 = fopen ("section7", "w");
  output_field ({eta,zb,maxa},fp, linear = true);
  for (double y = -10.; y <= 10.; y += 0.02)
  {
  fprintf (fp2, "%g %g\n", y, interpolate (maxa, 3., y));
  fprintf (fp5, "%g %g\n", y, interpolate (maxa, 9., y));
  }
  for (double x = -10.; x <= 10.; x += 0.02)
  fprintf (fp7, "%g %g\n", x, interpolate (maxa, x, 0.));
  }

*/

//Save of the movie, top view
event movie(i++)
{
output_ppm (eta,file = "eta.mp4",n = 2000,min=-2*AMP+H,max=2*AMP+H);
}

#endif
 /**
 ### Code's events and functions
 */

 /**
  * function z0
  */

double z0 (double q,double s)
{
  return (H - ((A*k)/om)*sin(-s*om+k*q)*sinh(k*H));
}


/**
 * output's events */
#if data_recover

//Calcul of Fourier transfrom at t0
event fourier (t = 0; t <t_end; t+=1.0)
{
  static FILE * fp = fopen ("fourier_t0.txt", "a");
  static FILE * fp2 = fopen ("third_wave_formation_time.txt","a");
  double vx,SR,SRB,SI,tb,px,py,pxb,pyb;
  int cy,cx;
  cx = 0;
  cy = 0;
  SRB = 0;
  pxb = 0;
  pyb = 0;
  tb = 0;
  for (px = 0; px <= POINT; px++)
    {
      for (py = 0; py <= POINT; py++)
	{
	  cx = 0;
	  vx = 0;
	  SR = 0;
	  SI = 0;
	  foreach(serial)
	  {
	    if ((vx != x))
	      {
		cy = 0;
		cx += 1;
		//printf(fp,"ok");
	      }
	    SR += eta[]*cos(-2*pi*(((px*cx*1.0)/POINT) + ((py*cy*1.0)/POINT)));
	    SI += eta[]*sin(-2*pi*(((px*cx*1.0)/POINT) + ((py*cy*1.0)/POINT)));
	    vx = x;
	    cy += 1;
	    //printf("%g %g %g %g\n",t,x,y,cos(-2*pi*(((px*cx*1.0)/POINT) + ((py*cy*1.0)/POINT))));
	  }
	  if ((pxb != px)); /*& ((t==0) || (t==t_end))/)
	    {
	      fprintf(fp,"\n");
	      fflush(fp);
	      //fprintf(fp2,"\n");
	      //fflush(fp2);
	    }
	  /*	  if ((pyb != py) & ((t==0) || (t==t_end)))
	    {
	      fprintf(fp2,"\n");
	      fflush(fp2);	     
	      }*/
	  pxb = px;
	  pyb = py;
	  if (((px != 0) || (py !=0)))
	    {
	      fprintf(fp,"%g %i %i %g %g %g\n",t,cx,cy,px,py,sqrt(pow(SR,2)+pow(SI,2)));
	      fflush(fp);
	      //printf("\n");
	      if (((sqrt(pow(SR,2)+pow(SI,2))) > 1) && ((sqrt(pow(SR,2)+pow(SI,2))) < 1000))
		{
		  if (tb != t)
		    {
		      //fprintf(fp2,"\n");
		      //fflush(fp2);
		    }
		
		  tb = t;
		  //printf("%g",tb);
		  fprintf(fp2,"\n\n");
		  fflush(fp2);
		  fprintf(fp2,"%g %g %g %g\n",t,px,py,0.0);
		  fprintf(fp2,"%g %g %g %g\n",t,px,py,sqrt(pow(SR,2)+pow(SI,2)));
		  fflush(fp2);
		}
	      /*else
		{
		fprintf(fp2,"%g %g %g %g\n",t,px,py,0.0);
		fflush(fp2);
		}*/
	      SRB =sqrt(pow(SR,2)+pow(SI,2));
		
	    }

	}
    }

}

#endif
/*
//Calcul of the error
event error (t = t_end)
{	 
  static FILE * fp = fopen ("error.txt", "a");
  static FILE * fp2 = fopen ("error_avg.txt", "a");
  double err_avg=0.0;

  foreach(serial)
  { 
    err_avg += (fabs(eta[]-z0(x,t))/AMP);
    fprintf(fp,"%g %g %g %i %i %g %i %g %g %g\n",x,L0,t,POINT,LEVEL,A,nl,fabs(eta[]-z0(x,t)),z0(x,t),((A*k)/om)*sinh(k*H));
    fflush(fp);  
  }
  fprintf(fp2,"%g %i %i %g %i %g %g\n",t,POINT,LEVEL,A,nl,CAM,((err_avg)/POINT));
  fflush(fp2);
}
*/

#if multiple_run_time
//Calculation of the period
event phase(i++)
{
#if multiple_waves
  if (t >= t_end*0.9)
#endif
    {
      int c = 0; 
      static FILE * fp = fopen ("phase.txt", "a");
      c = 0;
      foreach(serial)
      { 
	c += 1;
	if (c == 33410)
	  {
	    //if (c==2210){
	    if(eta[] - H > 0) s2 = 1;
	    if(eta[] - H < 0) s2 = 0;
	    if((s1!=s2) & (s1==0))
	      {
		s1=1;
		if (c2==1+nm*2)
		  {
		    Lz = zs + fabs(eta[]-H);
		    Dtb = t-told;
		    tc = t - (1 - (zs/Lz))*Dtb - tb;
		  }
		if (c2==1+nm*2)
		  {
#if multiple_waves
		    fprintf(fp,"%g %g %g %.10g %i %i %g %g %.6g %g %i %.10g %g\n",x,y,tc,(t-tb2),POINT,LEVEL,CAM,(2*k1-k2)*H,2*k1-k2,H,nl,eta[]-H,tb);  
		    fflush(fp);
#else
		    fprintf(fp,"%g %g %g %.10g %i %i %g %g %.6g %g %i %.10g %g\n",x,y,tc,(t-tb2),POINT,LEVEL,CAM,k1*H,k1,H,nl,eta[]-H,tb);  
		    fflush(fp);
#endif
		  }
		c2 += 1;
		tb2 = t;
		Lz = zs + fabs(eta[]-H);
		Dtb = t-told;
		tb = t - (1 - (zs/Lz))*Dtb;
	      }
	    if((s1!=s2) & (s1==1))
	      {
		s1 = 0;
	      }
	    zs  =fabs(eta[]-H);
	    told = t;
	  }
    
      }  
    }
}
#endif

//End condition
event end (t = t_end){}


