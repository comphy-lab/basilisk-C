/**
# Rupture de barrage Multi Couches

 exemple en C simple  Sans Viscosité
 
*/
#include <stdio.h>
#include <stdlib.h>
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

double t,dt,tmax,dx;
int nx,N;
/**
  Definition of the velocity, and of the discrete flux
*/
double C(double ug,double ud,double hg,double hd,double g)
{ double c;
    //Rusanov
    c=fmax(fabs(ug)+sqrt(g*hg),fabs(ud)+sqrt(g*hd));
    //c=fmax(  c,1*dx/2/dt);
    return c;
}
double FR1(double ug,double ud,double hg,double hd,double g,double c)
 { //double cR;
 //Rusanov 
   //cR=fmax(fabs(ug)+sqrt(g*hg),fabs(ud)+sqrt(g*hd));
  //c=1*dx/2/dt;
   return (hg*ug+hd*ud)*0.5-c*(hd-hg)*0.5;
 }
double FR2(double ug,double ud,double hg,double hd,double g,double c)
 { //double cR;
 //Rusanov 
  // cR=fmax(fabs(ug)+sqrt(g*hg),fabs(ud)+sqrt(g*hd));
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
// parameters --------------------------------------------------------------  
    dt=0.05/8;
    tmax=1;
    dx=0.1/8;
    nx=2*50*8;
    N=16;
    t=0;
    fprintf(stderr,"  ---------------------                         \n");
    fprintf(stderr,"                   <-- |                        \n");
    fprintf(stderr,"  ---------------------------                   \n");
    fprintf(stderr,"                             |--->              \n");
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
    {  x[i]=(i-nx/2)*dx;
   //dambreak
     Z[i]=0;
     h[i]=1*( x[i]<0)+.00;
     u[i]=0.00;
     Q[i]=u[i]*h[i];
    }
/** 
  $H =\Sigma_{\alpha=1}^N h_\alpha$  (2.6) 
  and each layer depth $h_\alpha$ is then deduced from the total water height by the relation $h_\alpha=l_\alpha H$
 (2.19)
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
           u_alpha[i][alpha]=0;
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
  flux corresponding to mass conservation accross the full layer $F_{p }= \Sigma_{\alpha=1}^N h_\alpha u_\alpha = Q$
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
 Final update
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
        for(int alpha=0;alpha<=N;alpha++)
        {
            // flux nul en entree sortie
            hn_alpha[0][alpha]=hn_alpha[1][alpha];
            un_alpha[0][alpha]=un_alpha[1][alpha];
            hn_alpha[nx][alpha]=hn_alpha[nx-1][alpha];
            un_alpha[nx][alpha]=un_alpha[nx-1][alpha];
        }
        hn[0]=hn[1];
        hn[nx]=hn[nx-1];
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
         $z_{\alpha+1/2} = z_b + \Sigma_{j=1}^\alpha h_j$ (2.7)
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
  if(it%10==0){
      /* Saving the fields */
      FILE *g= fopen("solxhQt.OUT", "w");
      for (int i=0; i<=nx;i++)
      {      fprintf(g,"%lf %lf %lf %lf %lf %lf \n",x[i],h[i],Q[i],t,h_alpha[i][1],h_alpha[i][2]);}
      fprintf(g,"\n");
      fclose(g);
  }
        printf("t=%lf\n",t);
        printf("set xlabel \"x\" ; set title \"t=%lf \"\n",t);
        double y1=-.1,y2=1.25;
        printf("set label \"*\" at 0,.296296296296296 \n");
        printf("h(x)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t)) \n");
        printf("p[%lf:%lf][%lf:%lf] h(x) t'h exact','-' u 1:2 t'Q'  w l,'-'  t'h' w lp,'-'  t'h_1' w lp,'-'  t'Z' w l linec -1\n ",
			   -nx*dx/2,nx*dx/2,y1,y2);
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
cc -O3 -ffast-math -std=c99 -lm svdbmult.c  -o svdbmult
~~~

pour lancer le programme:

~~~bash
 ./svdbmult | gnuplot
~~~ 

##Results



un fichier appelé solxhQt.OUT avec x h Q et est créé, pour le tracer, on lance gnuplot:

 In this case there is no viscosity, so no coupling, each layer behaves as it is alone, we recover the exact Ritter solution of dam break.
 

 
~~~gnuplot     
  set ylabel "h"; set xlabel "x";  
  set label "+" at 0,.296296296296296 
  h(x,t)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t))
  p  'solxhQt.OUT' u 1:2 w p, h(x,1)
~~~ 

## Links

  * [damb.c]() the same with only one layer 

 
##Bibliography

 * Emmanuel Audusse, Marie-Odile Bristeau, Benoıt Perthame, and Jacques Sainte-Marie
   A multilayer saint-venant system with mass exchanges for shallow water flows. derivation and numerical validation
   ESAIM: M2AN 45 (2011) 169–200 DOI: 10.1051/m2an/2010036
 
 
 
 
PYL version 1, Manchester, fev 2016
*/






 



 
