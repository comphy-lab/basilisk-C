/**
# SV Multilayer hydrolic Jump

 this is not Basilisk, example in simple C
 
*/
#include <stdio.h>
#include <stdlib.h>
//#include "math.h"  // uncomment 
#include <math.h>
#include <string.h>
// Emmanuel Audusse, Marie-Odile Bristeau, Benoıt Perthame, and Jacques Sainte-Marie
// A MULTILAYER SAINT-VENANT SYSTEM WITH MASS EXCHANGES FOR SHALLOW WATER FLOWS. DERIVATION AND NUMERICAL VALIDATION
//
// From NR
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

double*x=NULL,*h=NULL,*u=NULL,*Q=NULL,*Z=NULL,*Fp=NULL,*hn=NULL,*un=NULL,*c0=NULL;
double**h_alpha=NULL,**u_alpha=NULL,**Q_alpha=NULL,*l_alpha=NULL;
double**z_alphap12=NULL,**u_alphap12=NULL,**G_alphap12=NULL,**w_alphap12=NULL;

double t,dt,tmax,dx,nu;
double hi,hf,U0;

int nx,N;
/**
  Definition of the velocity, and of the descrete flux, simple Rusanov
*/
double C(double ug,double ud,double hg,double hd,double g)
{ double c;
    //Rusanov
    c=fmax(fabs(ug)+sqrt(g*hg),fabs(ud)+sqrt(g*hd));
    //c=fmax(  c,1*dx/2/dt);
    return c;
}
double FR1(double ug,double ud,double hg,double hd,double g,double c)
 { double cR;
 //Rusanov 
   cR=fmax(fabs(ug)+sqrt(g*hg),fabs(ud)+sqrt(g*hd));
  //c=1*dx/2/dt;
   return (hg*ug+hd*ud)*0.5-c*(hd-hg)*0.5;
 }
double FR2(double ug,double ud,double hg,double hd,double g,double c)
 { double cR;
 //Rusanov 
   cR=fmax(fabs(ug)+sqrt(g*hg),fabs(ud)+sqrt(g*hd));
 //  c=1*dx/2/dt;
   return (ug*ug*hg + g*hg*hg/2. + ud*ud*hd + g*hd*hd/2.)*0.5 - c*(hd*ud-hg*ug)*0.5;
 }
  
/*     -------------------------------------------------    */        
int main (int argc, const char *argv[]) {
    int  it=0;
    double**fp=NULL,*Fp=NULL,**fd=NULL,**un_alpha=NULL,**hn_alpha=NULL;
/**
## Initialisations
*/
// parametres --------------------------------------------------------------  
    dt=0.001;
    tmax=5;
    dx=0.005;
    nx=2*50*2-(.14/dx);
  N=25;//128*2;
    t=0;
    nu=1;
    hi=.2;
    hf=.02;
    U0=5.5;
    
    fprintf(stderr,"                      ---------------------     \n");
    fprintf(stderr,"                     |       -->   -->          \n");
    fprintf(stderr,"  --------------------                          \n");
    fprintf(stderr,"     --->     --->           ->    ->           \n");
    fprintf(stderr,"  ----------------------------------------------\n");
/* ------------------------------------------------------------------------*/      
    x= (double*)calloc(nx+1,sizeof(double));
    h= (double*)calloc(nx+1,sizeof(double));
    hn=(double*)calloc(nx+1,sizeof(double));
    Q= (double*)calloc(nx+1,sizeof(double));
    Z= (double*)calloc(nx+1,sizeof(double));
    u= (double*)calloc(nx+1,sizeof(double));
    un=(double*)calloc(nx+1,sizeof(double));
    Fp=(double*)calloc(nx+1,sizeof(double));
    c0=(double*)calloc(nx+1,sizeof(double));
    
    l_alpha= (double*)calloc(N+1,sizeof(double));
/**
We consider that the flow domain is divided in the vertical direction into N layers 
of thickness $h_\alpha$ with N + 1   interfaces $z_{\alpha+1/2}$
*/
    h_alpha=alloc_2d_double(nx+1,N+1);
    u_alpha=alloc_2d_double(nx+1,N+1);
    Q_alpha=alloc_2d_double(nx+1,N+1);
    hn_alpha=alloc_2d_double(nx+1,N+1);
    un_alpha=alloc_2d_double(nx+1,N+1);
    fp=alloc_2d_double(nx+1,N+1);
    fd=alloc_2d_double(nx+1,N+1);
    u_alphap12=alloc_2d_double(nx+1,N+1);
    w_alphap12=alloc_2d_double(nx+1,N+1);
    z_alphap12=alloc_2d_double(nx+1,N+1);
    G_alphap12=alloc_2d_double(nx+1,N+1);
// initialisation cond init ----------------------------      
  for(int i=0;i<=nx;i++)
    {x[i]=0.14+(i)*dx;
     Z[i]=0;
     h[i]=hi+0*fmax(1.8*x[i],1);
     u[i]=U0;//1.1/hi;
     Q[i]=u[i]*h[i];
    }
/** 
  $h =\Sigma_{\alpha=1}^N h_\alpha$
  and each layer depth $h_\alpha$ is then deduced from the total water height by the relation $h_\alpha=l_\alpha H$
 (2.20)
*/       
     for(int alpha=1;alpha<=N;alpha++)
     {
         l_alpha[alpha]= (1.00)/N;
         for(int i=0;i<=nx;i++)
                  h_alpha[i][alpha]=h[i]*l_alpha[alpha];
     }
/**
  We consider the average velocities $u_\alpha$, $\alpha = 1,...,N $
*/     
   for(int alpha=0;alpha<=N;alpha++)
     {
       for(int i=0;i<=nx;i++)
       {
           u_alpha[i][alpha]=u[i]*2.*alpha/N;
           w_alphap12[i][alpha]=0;
       }
     }
/**
the value of the velocity at the interface $z_{\alpha+1/2}$ is $u_{\alpha+1/2} = u(x,z_{\alpha+1/2},t)$ eq. (2.10)
*/     
   for(int alpha=0;alpha<=N;alpha++)
     { 
       for(int i=0;i<=nx;i++)
           u_alphap12[i][alpha]=0;
     }
/**
  The definition of :
  $G_{\alpha+1/2}$=$\partial_t z_{\alpha+1/2} + u_{\alpha+1/2} \partial_t z_{\alpha+1/2}$-$w(z,z_{\alpha+1/2},t)$ (2.13)
  The relation (2.13) gives the mass flux leaving/entering the layer $\alpha$ through the interface $z_{\alpha+1/2}$
*/     
   for(int alpha=0;alpha<=N;alpha++)
     { 
       for(int i=0;i<=nx;i++)
           G_alphap12[i][alpha]=0;
     }
/**
## Time advance
*/
    while(t<tmax){   // boucle en temps
        t+=dt;
        it++;
/**
 Estimating a global wave velocity
*/
        for(int i=1;i<=nx;i++)
            c0[i]=C(u[i-1],u[i],h[i-1],h[i],1.);
        
/**
  flux corresponding to mass conservation accross the full layer $F_{p}= \Sigma_\alpha h_\alpha u_\alpha = Q $
*/
        for(int i=1;i<=nx;i++)
            Fp[i]=FR1(u[i-1],u[i],h[i-1],h[i],1.,c0[i]);
        for(int i=1;i<nx;i++)
            hn[i]=h[i]-dt*(Fp[i+1]-Fp[i])/dx;   //conservation de la masse
/**
  Flux of mass $f_{p\alpha}$ in each layer (2.11): $f_{p\alpha}= h_\alpha u_\alpha $
*/
        for(int alpha=1;alpha<=N;alpha++)
            for(int i=1;i<=nx;i++)
                fp[i][alpha]=FR1(u_alpha[i-1][alpha],u_alpha[i][alpha],h_alpha[i-1][alpha],h_alpha[i][alpha],1/l_alpha[alpha],c0[i]);
/**
 Flux of  of momentum
    $f_{d\alpha}= h_\alpha u_\alpha^2 + \frac{1}{l_\alpha} \frac{g h_\alpha^2}{2} $
       eq.(2.21) has a $l_\alpha$ in flux $\partial_x(h_\alpha u_\alpha^2 + \frac{1}{l_\alpha} \frac{g h_\alpha^2}{2})$
*/
        for(int alpha=1;alpha<=N;alpha++)
            for(int i=1;i<=nx;i++)
                fd[i][alpha]=FR2(u_alpha[i-1][alpha],u_alpha[i][alpha],h_alpha[i-1][alpha],h_alpha[i][alpha],1/l_alpha[alpha],c0[i]);
/**
have to compute $G$, Using (2.18), (2.19), the expression of $G_{\alpha+1/2}$  given by (2.16) can also be written as (2.22)
 $$G_{\alpha+1/2}= \Sigma_{j=1}^\alpha(\partial_x(F_{pj}) -l_j \partial_x(f_{pj}))$$
*/
    for(int i=0;i<nx;i++)
        {
            for(int alpha=1; alpha<N; alpha++)
            {
                double dxhjuj=0;
                double slj=0;
                for(int j=1;j<=alpha;j++)
                {
                    slj=slj+l_alpha[j];
                    dxhjuj =  dxhjuj + (fp[i+1][j]-fp[i][j])/dx;
                }
                G_alphap12[i][alpha]= -slj*(Fp[i+1]-Fp[i])/dx + dxhjuj ;
            }
//          fprintf(stderr,"%d dh/dt=%lf  G %lf %lf %lf\n",i,-(Fp[i+1]-Fp[i])/dx,G_alphap12[i][0],G_alphap12[i][1],G_alphap12[i][2]);
        }
/**
$G_{1/2}=0$ $G_{N+1/2}=0$  (2.15)
the equations just express that there is no loss/supply of mass through the bottom and the free surface.
*/
        for(int i=1;i<nx;i++)
            G_alphap12[i][0]=0;
        for(int i=1;i<nx;i++)
            G_alphap12[i][N]=0;
/**
The velocities $u_{\alpha+1/2}$ are obtained using an upwinnding (2.23),
 if $G_{\alpha+1/2}>=0$ then $u_{\alpha+1/2}= u_{\alpha}$,
 if $G_{\alpha+1/2}<0$ then  $u_{\alpha+1/2}= u_{\alpha+1}$
*/
        for(int alpha=1;alpha<N;alpha++)
        {
            for(int i=1;i<nx;i++)
            {
                if( G_alphap12[i][alpha]>=0) {
                    u_alphap12[i][alpha] =  u_alpha[i][alpha];
                }else{
                    u_alphap12[i][alpha] =  u_alpha[i][alpha+1];
                }
               // u_alphap12[i][alpha]=(u_alpha[i][alpha]+u_alpha[i][alpha])/2;
            }
        }
/** 
 Final inviscid update
  $q_\alpha^{n+1} $
  $=q_\alpha^{n} - $
 $\Delta t(\partial_x (f_{d \alpha}) + u_{\alpha+1/2} G_{\alpha+1/2} - u_{\alpha-1/2} G_{\alpha-1/2})$
*/
      for(int alpha=1;alpha<=N;alpha++)
        {
            double q;
        for(int i=1;i<nx;i++)
            {
            hn_alpha[i][alpha]=l_alpha[alpha]*hn[i];
            if(hn_alpha[i][alpha]>0.){                             //conservation qunatité de mouvement
                q=h_alpha[i][alpha]*u_alpha[i][alpha]
                      -dt*(fd[i+1][alpha]-fd[i][alpha])/dx
                      +dt*(u_alphap12[i][alpha]*G_alphap12[i][alpha]-u_alphap12[i][alpha-1]*G_alphap12[i][alpha-1]);
                un_alpha[i][alpha]=q/hn_alpha[i][alpha];}
            else{
                un_alpha[i][alpha]=0.;}
            }
        }
/** 
 Simple Neumann BC
*/ 
     /*
        FILE *g = fopen("F.IN", "r");
        fscanf(g,"hi=%lf          \n",&hi);
        fscanf(g,"hf=%lf          \n",&hf);
        fscanf(g,"U0=%lf          \n",&U0);
        fscanf(g,"nu=%lf          \n",&nu);
        fclose(g);
       */  
        for(int alpha=0;alpha<=N;alpha++)
        {
            // flux nul en entree sortie
            hn_alpha[0][alpha]= l_alpha[alpha]*hi;//hn_alpha[1][alpha];
            un_alpha[0][alpha]= U0;//un_alpha[1][alpha];
            hn_alpha[nx][alpha]=l_alpha[alpha]*hf;//
          //  hn_alpha[nx][alpha]=hn_alpha[nx-1][alpha];
            un_alpha[nx][alpha]=un_alpha[nx-1][alpha];
        }
        hn[0]=hi;//hn[1];
        hn[nx]=hf;//hn[nx-1];
/**
## Viscous step
*/
#define vis
#ifdef vis  
        
        double *uN=(double*)calloc(N+2,sizeof(double));
        double *diags = (double*)calloc(N+1,sizeof(double));
        double *diagp = (double*)calloc(N+1,sizeof(double));
        double *diagi = (double*)calloc(N+1,sizeof(double));
        double *rhs   = (double*)calloc(N+1,sizeof(double));
        static double dz,NU;

        for(int i=1;i<nx;i++)
        {
            dz = hn[i]/N;
            if(hn[i]>0.){
                NU = dt*nu/dz/dz;
        /**
          Consturtion of the tridiagonal matrix corresponding to the "heat equation"
         $$(u_j^n - u_j^o)/\Delta t =   \frac{\nu }{\Delta z^2}  (u_{j-1}^n-2 u_j^n +u_{j+1}^n ) $$
         written as:
         
        $$  - NU u_{j-1}^n + (1+ 2 NU)  u_{j}^n - NU u_{j+1}^n = u_j^o$$
        
       
        
         with $NU = \frac{\nu \Delta t}{\Delta z^2}$
         */
        for(int j=2; j<=N; ++j){
            diags[j] =  -NU;
            diagp[j] = 1+2*NU;
            diagi[j] =  -NU;
            rhs[j] =    un_alpha[i][j];
        }
        /**
         values at the wall
        $u_1^n-u_1^o = \frac{\nu \Delta t}{\Delta z^2}  (u_{-1}^n-2 u_1^n +u_2^n$ 
         with $u_1^n=-u_{-1}^n$ so that the velocity is 0 at the wall
        */
            diagp[1]= 1+3*NU;
            diags[1]= -NU;
            rhs[1]= un_alpha[i][1];
        /**
         value at the surface,
                 $u_N^n-u_N^o = \frac{\nu \Delta t}{\Delta z^2}  (u_{N-1}^n-2 u_N^n +u_{N+1}^n$ with
         here   $u_N^n=-u_{N+1}^n$ so that the derivatice is 0 at the top
        */
            diagi[N]= -NU;
            diagp[N] = 1 + NU;
            rhs[N] = un_alpha[i][N];
            /*     inversion de matrice tridiag                               */
            /*     ai ui-1 + bi ui +ci ui+1 = rhs i                           */
            /*     diagi[j] u[j-1] +  diagp[j] u[j] + diags[j] u[j+1]= rhs[j] */
            /*    	u[j] = a[j] * u[j + 1] + b[j];                            */
            /*     Peyret Taylor p 20                                         */
            /*     l indice varie de 1 a N                                    */
            /* Using the tridiagonal matrix algorithm (TDMA), Thomas Method   */
            /* http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm) */
            /* Forward elimination phase  */
            for (int j = 2; j <= N; j++) {
                double m = diagi[j]/diagp[j-1];
                diagp[j] -= m*diags[j-1];
                rhs[j] -= m*rhs[j-1];
            }
            /* Backward substitution phase  */
            uN[N] = rhs[N]/diagp[N];
            for (int j = N - 1; j >= 1; j--)
                uN[j] = (rhs[j] - diags[j]*uN[j+1])/diagp[j];
             for (int j = 1; j <= N; ++j)
                 un_alpha[i][j] = uN[j];
        }
        }
        free(uN);
        free(diags);free(diagp);free(diagi);free(rhs);
#endif
       //final swap
        for(int i=0;i<=nx;i++)
        {
            h[i]=hn[i];
            Q[i]=0;
            for(int alpha=1;alpha<=N;alpha++)
            {
                h_alpha[i][alpha]=hn_alpha[i][alpha];
                u_alpha[i][alpha]=un_alpha[i][alpha];
                Q[i]=Q[i]+hn_alpha[i][alpha]*un_alpha[i][alpha];
            }
            if(hn[i]>0.){  un[i]=Q[i]/hn[i];} else {un[i]=0;}
            u[i]=un[i];
            /**
             $z_{alpha+1/2} = Z + \Sigma_{j=1}^\alpha h_j$ (2.7)
             */
            z_alphap12[i][0] = Z[i];
            for(int alpha=1;alpha<=N;alpha++){
                for(int j=1;j<=alpha;j++)
                {
                    z_alphap12[i][j]=z_alphap12[i][j-1]+h_alpha[i][j];
                }
            }
        }  
/**
## Results
*/
  if(it%100==0){
      /* Saving the fields */
    //  if(i==nx/2){
          FILE *gq= fopen("profil0.OUT", "w");
          for (int j = 1; j <= N; ++j) {
              double s;
              s=0;
              fprintf(gq," %lf %lf %lf %lf\n",j*dz,un_alpha[nx/2][j],w_alphap12[nx/2][j],t);
          }
          fclose(gq);
     // }
      
      
      FILE *g= fopen("solxhQt.OUT", "w");
      for (int i=0; i<=nx;i++)
      {      fprintf(g,"%lf %lf %lf %lf %lf \n",x[i],h[i],Q[i],t,u_alpha[i][1]/h_alpha[i][1] );}
      fprintf(g,"\n");
      fclose(g);
 
        printf("set key left;t=%lf\n",t);
        printf("set arrow nohead  from 0,0 to 0.2,0 lw 1 \n");
        printf("set arrow nohead  from 0,0 to 0.2,0.2*1.8138 lw 1 \n");
        printf("set xlabel \"x\" ; set title \"t=%lf \"\n",t);
        double y1=-.1,y2=2,scu=.05;;
        printf("p[%lf:%lf][%lf:%lf]  '-' u 1:2 t'Q'  w l linec 2,'-'  t'h' w l linec 3,'-'  t'h_1' w l,'-'  t'Z' w l linec -1,'-' not w l linec 1,'-' not w l linec 1,'-' w l not linec 1,'-' w l t 'tau', 1\n ",
			   0*x[0],x[nx],y1,y2);
		for (int i=0; i<=nx;i++)
		{
			printf("%lf %lf \n",x[i],Q[i]);}
		printf("e \n");
		for (int i=0; i<=nx;i++)
		{
			printf("%lf %lf \n",x[i],h[i]+Z[i]);}
		printf("e \n");
		for (int i=0; i<=nx;i++)
		{
			printf("%lf %lf \n",x[i],h_alpha[i][1]+Z[i]);}
		printf("e \n");
		for (int i=0; i<=nx;i++)
		{
			printf("%lf %lf \n",x[i],Z[i]);}
		printf("e \n");
        
        
        for (int j=1; j<=N;j++)
		{
			printf("%lf %lf \n",x[nx/4]+scu*u_alpha[nx/4][j],z_alphap12[nx/4][j]);}
        printf("%lf %lf \n",x[nx/4] ,z_alphap12[nx/4][N]);
        printf("%lf %lf \n",x[nx/4],0.);
		printf("e \n");
        for (int j=1; j<=N;j++)
		{
			printf("%lf %lf \n",x[nx/2]+scu*u_alpha[nx/2][j],z_alphap12[nx/2][j]);}
        printf("%lf %lf \n",x[nx/2] ,z_alphap12[nx/2][N]);
        printf("%lf %lf \n",x[nx/2],0.);
		printf("e \n");
        
        for (int j=1; j<=N;j++)
		{
			printf("%lf %lf \n",x[3*nx/4]+scu*u_alpha[3*nx/4][j],z_alphap12[3*nx/4][j]);}
        printf("%lf %lf \n",x[3*nx/4] ,z_alphap12[3*nx/4][N]);
        printf("%lf %lf \n",x[3*nx/4],0.);
		printf("e \n");
      for (int i=0; i<=nx;i++)
      {
          printf("%lf %lf \n",x[i],u_alpha[i][1]/h_alpha[i][1]);}
      printf("e \n");
  }
    //getchar();
      }

     free(x);
     free(fp);
     free(fd); 
     free(un);
     free(hn);
     free_2d_double(h_alpha);
    return 0;
}
/**

##Run

Programme en C simple fichier (en cours on a regardé un fichier plus compliqué). Cas de rupture de barrage flux Rusanov (essayer HLC) compilation
on le compile et on crée l’exécutable db

~~~bash
 cc -O3 -ffast-math -std=c99 -lm svdbvismult_hydrojump.c -o svdbvismult_hydrojump
~~~

pour lancer le programme:

~~~bash
 ./svdbvismult_hydrojump | gnuplot
~~~ 

 

##Results


 profil de vitesse avec z u v  et t est créé, pour le tracer, on lance gnuplot:


~~~bash
 set ylabel "z"; set xlabel "u";
 p 'profil0.OUT' u 1:2 not w l linec 3
~~~ 
   

~~~gnuplot
p[0:][:2]'solxhQt.OUT' t 'h' w l,''u 1:5 t'tau'w l,0 not
~~~
   
 Comparing with Higuera
 
~~~gnuplot
set output 'f.png'
 reset
 X0=169
 X1=604
 Y0=222.24
 Y1=528
 unset tics
 p[0:][0:605]'../Img/higueragraph.png' binary filetype=png   with rgbimage not,\
 'solxhQt.OUT'u (X0+$1*(X1-X0)):($2/2*(Y1-Y0)+Y0) t 'h' w lp linec 3,\
 '' u (X0+$1*(X1-X0)):($5/15*(Y1-Y0)+Y0) t 'tau' w l
  
~~~

## Links
  
  * [http://basilisk.fr/src/test/higuera.c]() with Basilisk
  
## Bibliography

 * Emmanuel Audusse, Marie-Odile Bristeau, Benoıt Perthame, and Jacques Sainte-Marie
   A multilayer saint-venant system with mass exchanges for shallow water flows. derivation and numerical validation
   ESAIM: M2AN 45 (2011) 169–200 DOI: 10.1051/m2an/2010036
 
 * Higuera, F. 1994 The hydraulic jump in a viscous laminar flow. J. Fluid Mech. 274, 69–92.
 
 * PYL
 ["exam"](http://www.lmm.jussieu.fr/~lagree/COURS/M2MHP/exam2013.pdf)
 
  
* Francesco De Vita, Pierre-Yves Lagrée, Sergio Chibbaro, Stéphane Popinet
[Beyond Shallow Water: appraisal of a numerical approach to hydraulic jumps based upon the Boundary Layer Theory](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/jump19.pdf) to appear
 
 
PYL version 1, Manchester, Besançon fev 2016, Twente 2018
*/


 
