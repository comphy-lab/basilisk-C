/**
# Rupture de barrage bien balancée / Dam Break Well Balanced
 
 An example in simple C (not a `Basilisk` example).
 We solve here the dambreak problem  with a topography down stream
 and a wall
 with HLL flux and Well Balanced scheme.
 
 At the left a impermeable wall, at the right a constant slope.
 
 */

/* | wall                                    /
/* |––––––––––—––—––—                      / constant slope
/* |  initial        |                   /
/* |  water          |                 /
/* |  at rest        |               /
/*  ------------------------------
/*   x=-L/2                         x=0           x=L
*/
/**
see details of the code
http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf

 usage :  cc svdbWB.c -lm  ;./a.out | gnuplot
 results are in file "sol.OUT" plotted at the end
 
The same in python is there:
https://basilisk.fr/sandbox/M1EMN/PYTHON/svdbWB.py

 
*/
#include <stdio.h>
#include <stdlib.h>
//#include "math.h"
#include <math.h>
#include <string.h>
// pyl 05/12 
//     03/26
//
// M1EMN lecture Lagree
// "A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows"
// Emmanuel Audusse, Franccois Bouchut, Marie-Odile Bristeau, Rupert Klein and Benoıt Perthame
// see as well These O. Delestre
    double*x=NULL,*h=NULL,*u=NULL,*Q=NULL,*U=NULL,*Z=NULL,*dZ=NULL,*Zp=NULL;
    double t,dt,tmax,dx,Cf,h0,xbosse,xtas,h0,alpha,e;
    int cas,nx,icl;
// the heart of WB technology, reconstruction hydrostatique
 void  reconsetat(double hl,double hr,double dz,double *hil,double *hir)
{
      *hil=fmax(0.,hl-fmax(0.0,dz));
      *hir=fmax(0.,hr-fmax(0.0,-dz));
}
// flux h  
double FHLL1(double ug,double ud,double hg,double hd)
 { double c1,c2,f1;
 //HLL  flux for h
     c1=fmin(ug -sqrt(hg),ud-sqrt(hd));
     c2=fmax(ug +sqrt(hg),ud+sqrt(hd));
      if (c1>=0.){
      f1=hg*ug;}
      else{ 
      	if ((c1<0.)&&(0.<c2)) {
      f1=(c2*hg*ug-c1*hd*ud)/(c2-c1)+c1*c2*(hd-hg)/(c2-c1);
      	} else {
      f1=hd*ud;} }    
   return f1;
 }
// flux q=hu 
double FHLL2(double ug,double ud,double hg,double hd)
 { double c1,c2,f2;
 //HLL flux for hu
     c1=fmin(ug -sqrt(hg),ud-sqrt(hd));
     c2=fmax(ug +sqrt(hg),ud+sqrt(hd));
      if (c1>=0.){
      f2=hg*ug*ug+hg*hg/2;}
      else{ 
      	if ((c1<0.)&&(0.<c2)) {
      f2=(c2*(hg*ug*ug+hg*hg/2)  -c1*(hd*ud*ud+hd*hd/2))/(c2-c1) +c1*c2*(hd*ud-hg*ug)/(c2-c1);
      	} else {
      f2=hd*ud*ud+hd*hd/2.;} }    
   return f2;
 }

/** 
------------------------------------------------------------------------*/
void ecrit()
{    
	extern double  dx,t; 
  extern int nx;
	int i;
	double y1,y2;
	extern double *x,*h,*Z,*Q,e;
 
	/* Saving the fields */ 
	 // sortie GNUPLOT		
	
	printf("t=%lf\n",t);
	printf("set xlabel \"x\" ; set title \"t=%lf  masse=%lf \"\n",t,e);
    y1=-.25,y2=fmin(2,2);
 	 
	printf("p[%lf:%lf][%lf:%lf] '-' u 1:2 t'Q'  w l,'-'  t'h' w d,'-'  t'Z' w l linec -1\n ",
			   -nx*dx/2,nx*dx/2,y1,y2); 
		for (i=0; i<=nx;i++)
		{
			printf("%lf %lf \n",x[i],Q[i]);}
		printf("e \n");
		for (i=0; i<=nx;i++)
		{
			printf("%lf %lf \n",x[i],h[i]+Z[i]);}
		printf("e \n");
		for (i=0; i<=nx;i++)
		{
			printf("%lf %lf \n",x[i],Z[i]);}
		printf("e \n");	
} 
     
/** 
------------------------------------------------------------------------*/
    
int main (int argc, const char *argv[]) {
    int  i,it=0;         
  double*fp=NULL,*fd=NULL,*hg=NULL,*hd=NULL,*hid=NULL,*hig=NULL,*ug=NULL,*ud=NULL,*un=NULL,*hn=NULL;
    double q,f1,f2;
    
 // init 
    dt=0.01;
    tmax=40; //64
    dx=0.025;
    nx=2000;
    Cf=0.01*0;
    h0=1;
    alpha=0.1;
    xbosse =0;
    xtas=-2;

 // affectation 
      x= (double*)calloc(nx+2,sizeof(double));
      h= (double*)calloc(nx+2,sizeof(double));
      Q= (double*)calloc(nx+2,sizeof(double));
      u= (double*)calloc(nx+2,sizeof(double));
      Z= (double*)calloc(nx+2,sizeof(double));
      dZ=(double*)calloc(nx+2,sizeof(double));
      Zp=(double*)calloc(nx+2,sizeof(double));
      fp=(double*)calloc(nx+2,sizeof(double));
      fd=(double*)calloc(nx+2,sizeof(double));
      hg=(double*)calloc(nx+2,sizeof(double));
      hd=(double*)calloc(nx+2,sizeof(double));
      hig=(double*)calloc(nx+2,sizeof(double));
      hid=(double*)calloc(nx+2,sizeof(double));
      ug=(double*)calloc(nx+2,sizeof(double));
      ud=(double*)calloc(nx+2,sizeof(double));
      un=(double*)calloc(nx+2,sizeof(double));
      hn=(double*)calloc(nx+2,sizeof(double));
      
      FILE *g = fopen("sol.OUT", "w");
      fclose(g);
   // initialisation       
    for(i=0;i<=nx+1;i++)
    { x[i]=(i-nx/2)*dx;
	  Z[i]=alpha*(x[i]>xbosse)*(x[i]-xbosse);
      h[i]=h0*(x[i]<-10);
      u[i]=0;
      Q[i]=0;
	     }        
    // slope, topography
 for(i=0;i<nx+1;i++)
    {  dZ[i]= Z[i+1]-Z[i];	}
      dZ[nx]=0;

    while(t<=tmax){
      t=t+dt;
      it=it+1;
      
// variables pour la reconstruction
       for(i=0;i<=nx+1;i++)
    { 
      hd[i]=h[i];
      ud[i]=u[i];
      hg[i]=h[i];
      ug[i]=u[i];}
/**
-----  la reconstruction "hydrostatique" 
 */ 
 for(i=1;i<=nx+1;i++)
    { 
       reconsetat(hd[i-1],hg[i],dZ[i-1],&f1,&f2);
       hid[i-1]=f1;
       hig[i]=f2;
      }
    
	
// construction de flux intermediaires a l'interface : i-1/2  fp:premier  fd:deuxieme

   for(i=1;i<=nx+1;i++)
    {    
    	  fp[i]=FHLL1(ud[i-1],ug[i],hid[i-1],hig[i]);
    	  fd[i]=FHLL2(ud[i-1],ug[i],hid[i-1],hig[i]); 
    }
  
    
// boucle principale
//   U=(h,hu)
//  	Unouveau = Uancien  + dt (Flux gauche-Flux droit)
		
   for(i=1;i<nx+1;i++)
    {	
      hn[i]=h[i]- dt*(fp[i+1]-fp[i])/dx;   //conservation de la masse
		
  if(h[i]>0.){                             //conservation qunatité de mouvement
      q=h[i]*u[i]-dt*(fd[i+1]-fd[i])/dx 
       -dt*( hig[i]*hig[i]/2 - hg[i]*hg[i]/2
            + hd[i]*hd[i]/2- hid[i]*hid[i]/2)/dx;
    
     // frottement     
      if((q!=0)&&(Cf>0)) q = q/(1+dt*Cf/h[i]/h[i]);
	  
      un[i]=q/hn[i];}
      else{
      un[i]=0.;}
    }
		
	hn[0]=hn[1];
	un[0]=-un[1];  // wall
        //un[0]=un[1];  // free
		
	hn[nx+1]=hn[nx];
	un[nx+1]=un[nx];
	//u[nx+1]=u[nx]-sqrt(h[nx]) + sqrt(h[nx+1]);
	//; 
	e=0; // controle conservation masse
  for(i=0;i<=nx+1;i++)
    { 
      e=e + hn[i]*dx;
      h[i]=hn[i];
      u[i]=un[i];
      Q[i]=hn[i]*un[i];
    }

      if(it%200==0){
        ecrit();
        FILE *g = fopen("sol.OUT", "a");
        for (i=0; i<=nx;i++) {
          fprintf(g,"%lf %lf %lf %lf \n",x[i],h[i],Q[i],Z[i]);} 
         fprintf(g,"\n");
       fclose(g);}
      }
    return 0;
}
/**

~~~gnuplot
reset; 
set xlabel "x";  
 p[-25:25][0:2] 'sol.OUT'  u 1:($2+$4) t'eta'w l,'' u 1:4 t' zb'w l
~~~ 


*/

 
