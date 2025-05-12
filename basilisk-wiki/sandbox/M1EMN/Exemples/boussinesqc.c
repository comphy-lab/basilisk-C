/**
## Test Boussinesq Linearized

 
 This example is not Basilisk but standard C for M1 course. We use explicit resolution
 with centered finite derivatives to solve the Boussinesq linearized wave problem.
 
 
 
# Problem

The shallow water equations are non linear bit they are not dispersives. But water waves in moderate depths are dispersives.
To fix this problem, a simple way is to consider Boussinesq correction.
 
# Equations
 

We look at the solution of linearized Boussinesq equations as then come from analysis (see M2 course)
$$\frac{\partial u}{\partial t}  = -\frac{\partial \eta }{\partial x}$$
$$\frac{\partial \eta }{\partial t} = - \frac{\partial u}{\partial x} + 
  \frac{1}{3} \frac{\partial^3 \eta}{\partial x^2 \partial t }  $$
 note that a better point of view, which presserves mass conservation, is to write instead
 $$\frac{\partial \eta}{\partial t}  = -\frac{\partial u }{\partial x}$$
$$\frac{\partial u }{\partial t} = - \frac{\partial \eta}{\partial x} + 
  \frac{1}{3} \frac{\partial^3 u}{\partial x^2 \partial t }  $$
 as the problem is linearized, it is the same problem.
 The equations are solved using a factorisation of the $\frac{\partial^2 }{\partial x^2 }$ term: 
 $$\frac{\partial u}{\partial t}  = -\frac{\partial \eta }{\partial t}$$
 $$(1 - 
   \frac{1}{3} \frac{\partial^2 }{\partial x^2 } )
 \frac{\partial \eta }{\partial t} = - \frac{\partial u}{\partial x}  
 $$
 In version of $(1 - 
   \frac{1}{3} \frac{\partial^2 }{\partial x^2 } )$ will be implicit whereas  the scheme is explicit.
 Again, this is done in simple C for  comparison.  

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#define pi 3.151592
/** 
starting form here, it can be in a 
`invSOR.h`

--------------------------------------------------------------------------------
  
routine alloc and free de OOura pour créeer matrice n x n */ 
double **alloc_2d_double(int n1, int n2)
{
    double **dd, *d;
    int j; 
    dd = (double **) malloc(sizeof(double *) * n1); 
    d = (double *) malloc(sizeof(double) * n1 * n2); 
    dd[0] = d;
    for (j = 1; j < n1; j++) {
        dd[j] = dd[j - 1] + n2;
    }
    return dd;
}
void free_2d_double(double **dd)
{
    free(dd[0]);
    free(dd);
}
/** 
   solve    $A.X = B$
   by  "JOR/SOR  Jacobi and Gauss-Seidel  over-relaxation Method" 
   Alfio Quarteroni Riccardo Sacco Fausto Saleri,    page 136
*/ 
void invSOR(int sor,double **A,double *X,double *B,int n)
{ 
  int iter=0,i,j;
  double s,omega=.5;
  double  r0,err;
  static double *Xbiel=NULL,*R=NULL;  
  Xbiel=(double*)calloc(n,sizeof(double));
  R=(double*)calloc(n,sizeof(double)); 
  s=0;
  for(i=0;i<n;i++){
      for(j=0;j<n;j++)s+=A[i][j]*X[j];
      R[i]=B[i]-s;
   } 
   r0=0;
   for(i=0;i<n;i++)r0+=(R[i]*R[i]);
   r0=sqrt(r0);
   err=r0;
   do{  
    iter++;
    for(i=0;i<n;i++)Xbiel[i]=X[i];
    for(i=0;i<n;i++){
      s=0;  
      if(sor==1){  
      if(i>0){
      for(j=0  ;j<=i-1;j++){s+=A[i][j]*X[j]; }
      }
      if(i<n-1){
        for(j=i+1;j<n  ;j++){s+=A[i][j]*Xbiel[j];}
      }
      } 
      else{     
      for(j=0 ;j<n;j++){if(i!=j) {s+=A[i][j]*X[j]; }}
      }
      X[i]=  omega*(B[i]-s)/A[i][i]+(1-omega)*Xbiel[i] ;
     }
    for(i=0;i<n;i++){
          s=0;
      for(j=0 ;j<n;j++){  s+=A[i][j]*X[j];  }
      R[i]=B[i]-s;      
    } 
   err=0;
   for(i=0;i<n;i++)err+=(R[i]*R[i]);
    err=sqrt(err)/r0;   
      }
    while((err > 1e-12)&&(iter<20000));
  if((err > 1e-12)||(iter>20000-1)) fprintf(stderr," iter=%d err=%lf\n", iter,err*1e+12);
  free(Xbiel);
  free(R);
} 
/**
fin de la  routine d'inversion SOR

--------------------------------------------------------------------------------

end of `invSOR.h`

##Main

Main program, value of paerameters, initial wave
*/    
int main (int argc, const char *argv[]) {
    int  i,nx,it=0;//,sor=1;
  //  char file;    
    double*x=NULL,*u=NULL,*F=NULL,*uo=NULL,*eta=NULL,*etao=NULL,*detadt=NULL;
    double dt,dx,x0,L,s,t,sigma,tmax; 
    double **a=NULL;
    FILE *g;
 
     nx = 375;
     L = 125;
     x0=-L/3.; 
     dx = L/nx;
     dt = fmin(dx,.25) ;// OK  0.1; 
     t = 0;
     tmax=60;
     sigma=1./3; //   1./3 for real Boussinesq
      
      x=     (double*)calloc(nx+1,sizeof(double));
      u=     (double*)calloc(nx+1,sizeof(double));
      eta=   (double*)calloc(nx+1,sizeof(double));
      uo=    (double*)calloc(nx+1,sizeof(double));
      etao=  (double*)calloc(nx+1,sizeof(double));
      detadt=(double*)calloc(nx+1,sizeof(double));
      F=     (double*)calloc(nx+1,sizeof(double)); 
     fprintf(stderr,"solution de  du/dt = - dh/dx; \n");
     fprintf(stderr,"solution de  dh/dt = - du/dx + sigma d^3 h/dx^2dt \n");
     fprintf(stderr,"dx=%lf,  sigma*dt//dx^3=%lf  dt/dx=%lf\n",dx,sigma*dt/dx*dx*dx,dt/dx);
// init 
     for(i=0;i<=nx;i++)
     { u[i]=0;
       x[i]= x0 +i*dx;
       uo[i]=0;
       u[i]=uo[i];
       etao[i]=exp(-pi*(x[i])*(x[i]));;
      // etao[i]=.01*(1-tanh(x[i]));;         
       u[i]=etao[i];
       uo[i]=u[i];
       eta[i]=etao[i];
     } 
     g = fopen("sol.OUT", "w");
     fclose(g);
/**
Begin time loop while time less than `tmax` */
  do{  
    t=t+dt; 
    it++;
/** Explicit resolution 
with centered finite derivatives. 
 Compute $\dfrac{\partial \eta}{\partial x} \sim \dfrac{  \eta(x+\Delta x)-\eta(x-\Delta x)}{2  \Delta x }$ so that 
 $$u (x, t +\Delta t ) = u (x, t ) - \Delta t \dfrac{  \eta(x+\Delta x)-\eta(x-\Delta x)}{2  \Delta x }  $$
*/  
    for(i=1;i<nx-1;i++) u[i] = uo[i]   - dt*(eta[i+1]-eta[i-1])/2/dx ;                            
            u[0]=uo[1];
            u[nx-1]=0;
            u[nx-2]=0;
// swap               
        for(i=0;i<nx;i++) uo[i]=u[i];        
// prepare    (1 - (1/3) d^2/dx^2) dh/dt = - du/dx       
      a = alloc_2d_double(nx,nx);  
      s=sigma/(dx*dx);
      for(i=2;i<nx-2;i++){a[i][i-2]=0;a[i][i-1]=-s;a[i][i]=1+2*s;a[i][i+1]=-s;a[i][i+2]=0;}
      a[0][0]=1;  
      a[1][1]=1; 
      a[nx-2][nx-2]=1;
      a[nx-1][nx-1]=1;
      detadt[0]=0;
      detadt[1]=0;
      detadt[nx-2]=0.;
      detadt[nx-1]=0.;
/** compute implicitely the increase in time of $\eta$
$$F= -\dfrac{\partial u}{\partial x} \sim -\dfrac{  u(x+\Delta x)-u(x-\Delta x)}{2  \Delta x }  $$
*/ 
    for(i=1;i<nx-1;i++) F[i] = - (uo[i+1] -uo[i-1] )/dx/2;
    F[0]=0;
    F[nx-1]=0;
/** by inversion of 
$(1 - 
   \frac{1}{3} \frac{\partial^2 }{\partial x^2 } )
 f  = F$
 */
#ifdef bla
      invSOR(sor,a,detadt,F,nx);
    //  for(i=1;i<nx;i++)detadt[i]=F[i];
#else
      double eps=.0,detadtn=0.,omega=.25,s=1./3;
      do
      { eps=0;
           for(i=1;i<nx;i++)
          {
              detadtn = (dx*dx*F[i] + s*(detadt[i-1] + detadt[i+1]))/(2*s + dx*dx);
              eps = eps +  (detadt[i]-detadtn)*(detadt[i]-detadtn);
              detadt[i] =   omega * detadtn + (1-omega) * detadt[i] ;
          }
          // fprintf (stderr, "t=%g eps=%g \n",t, sqrt(eps));
      }while(sqrt(eps)>1.e-11);
   
#endif
      free_2d_double(a);
/** 
we update 
$\eta (x, t+\Delta t) =  \eta (x, t) + (\Delta t) f$
*/
    for(i=1;i<nx;i++) eta[i] = etao[i]   +  detadt[i]*dt ;
/** 
swap an simple BC */
    
      eta[0] = eta[1];  
      eta[nx-1]=0;  
      eta[nx-2]=0;
      eta[nx-3]=0;
      for(i=1;i<nx;i++) etao[i]=eta[i];
/**      affichage  */      
    if((it%4)==0){ 
    g = fopen("sol.OUT", "a");
  for (i=0; i<=nx;i++)
  {      fprintf(g,"%lf %lf %lf \n",x[i],u[i],t);}
  fprintf(g,"\n\n");
  fclose(g);
  
  fprintf(stderr,"  t=%lf\n",t);
  
  printf(" p[%lf:%lf][%lf:%lf]  '-' u 1:2 t'u'  w l linec 3 , '-' u 1:2 t'h' w l\n ",
         x[1],x[nx],-1.,1.); 
    for (i=0; i<=nx;i++)
    {
    printf("%lf %lf \n",x[i],u[i]);}
    printf("e \n");
     for (i=0; i<=nx;i++)
    {
    printf("%lf %lf \n",x[i],eta[i]);}
    printf("e \n");
    }
  
  }while(t<tmax);

  // liberation finale
   free(x);
   free(u);
   free(uo);
   free(eta);
   free(etao);
   free(F);
   free(detadt);  
}
/**
## Run
to compile with standard C

~~~bash
 cc boussinesqc.c -lm -o boussinesqc; ./boussinesqc | gnuplot
~~~ 

note that using `qcc` it works like with `cc`.


##Results
One analytical solution is the self similar solution
$$\eta = \eta_0((t/2)^{-1/3}) Airy(2^{1/3} (x-t)t^{-1/3},t) $$
 the Airy function is approximated by its two asymptotic descirptions:
 $$sin(2/3.*((-x)^{3./2}+\pi/4))/\sqrt{(\pi*\sqrt{(-x)})}\text{ and }exp(-2/3 (x)^{(3./2)})/\sqrt{(4*\pi*\sqrt x)}$$
 and is in gnuplot as well `airy`
~~~gnuplot
set xlabel "x"
plot [-20:]for [IDX=0:100:10] 'sol.OUT' i IDX u  1:2 not w lines 
~~~   
    
which gives $U(x,t)$ plotted here for t=0 1 2 3  and $-3<x<9$ 
~~~gnuplot
 dp3=2**(1./3)
 Airy(x) = (x<-.5? sin(2/3.*((-x)**(3./2)+pi/4))/sqrt(pi*sqrt(-x)): (x >.5 ? exp(-2/3.*(x)**(3./2))/sqrt(4*pi*sqrt(x)): NaN))
 set xlabel "(x-t)/t^{1/3}" 
 plot [-10:5]for [IDX=1:100:10] 'sol.OUT' i IDX u (($1-$3)/(($3/2)**(1./3))):($2*(($3/2)**(1./3))) not w lines, Airy(x),airy(x)
~~~
  

## Links
 
* [airy_watertrainfront.c]() non linear case with the multilayer code
 
 
 
## biblio
  
* Magnar Bjørkavåg, Henrik Kalisch "Wave breaking in Boussinesq models for undular bores" 
 Physics Letters A 375 (2011) 1570–1578
 
* [PYL](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEhoule.pdf) M1 course on waves

* [PYL](http://www.lmm.jussieu.fr/~lagree/COURS/M2MHP/kdv.pdf) Kdv Equation M2 course

*  Alfio Quarteroni, Riccardo Sacco, Fausto Saleri,   "JOR/SOR  Jacobi and Gauss-Seidel  over-relaxation Method" 
   Numerical Mathematics,  Springer 2000
 

   
 Version 14 Jan 2018, initial: Pont de Mauvert, janvier 2014  
 */
  
