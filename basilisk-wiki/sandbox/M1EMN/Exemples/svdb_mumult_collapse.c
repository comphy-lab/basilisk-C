/**
# Rupture de barrage Multi Couches pour avancée de Front en $\mu(I)$
 
Exemple en C simple

les numéros des équations renvoient à Audusse et al. 2011

## Code  
 */
#include <stdio.h>
#include <stdlib.h>
//#include "math.h"
#include <math.h>
#include <string.h>

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

double*x=NULL,*h=NULL,*u=NULL,*Q=NULL,*Z=NULL,*dZ=NULL,*Fp=NULL,*hn=NULL,*un=NULL,*c0=NULL;
double *Q2=NULL,*Gamma=NULL,*Y=NULL;
double *hd,*hg,*ud,*ug,*hid,*hig;
double**h_alpha=NULL,**u_alpha=NULL,**Q_alpha=NULL,*l_alpha=NULL;
double**z_alphap12=NULL,**u_alphap12=NULL,**G_alphap12=NULL,**w_alphap12=NULL;

double t,dt,tmax,dx,nu;
double dry=1e-5,dry2=0.;
int nx,N;

double mus=.35,I0=.4,Deltamu=.21,d_g=0.043,tauy;
double tantheta=(-.45);
/**
 Definition of the velocity, and of the descrete flux
 */
void  reconsetat(double hl,double hr,double dz,double *hil,double *hir)
{
    *hil=fmax(0.,hl-fmax(0.0,dz));
    *hir=fmax(0.,hr-fmax(0.0,-dz));
}

double C(double ug,double ud,double hg,double hd,double g)
{ double c;
    //Rusanov
    c=fmax(fabs(ug)+sqrt(g*hg),fabs(ud)+sqrt(g*hd));
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
double U_Bagnold(double z, double h0,double tantheta)
{
    double u,Itheta;
    Itheta=(fabs(tantheta)-mus)/(mus+Deltamu-fabs(tantheta));
    u = (2./3)*Itheta*sqrt(d_g *pow(h0/d_g,3.))*(1-pow(1-(z/h0),3./2));
    return u;
}
/*     -------------------------------------------------    */
int main (int argc, const char *argv[]) {
    int  it=0;
    double**fp=NULL,*Fp=NULL,**fd=NULL,**un_alpha=NULL,**hn_alpha=NULL;
    double **hid_a,**hig_a,**ud_a,**ug_a,**hd_a,**hg_a;
    /**
## Initialisations
     
     parametres */
    dt=0.0025;
    tmax=350;      tmax=35;
    dx=0.02;
    nx=1000*5;     nx=2500;
    N=64*2;        N=32;
    t=0;
    nu=1;
    fprintf(stderr,"  ---------------------                       \n");
    fprintf(stderr,"                        \\  -->                \n");
    fprintf(stderr,"                         \\ -->                \n");
    fprintf(stderr,"                           |--->              \n");
    fprintf(stderr,"  --------------------------------------------\n");
    fprintf(stderr,"%lf %lf\n", U_Bagnold(1,1,tantheta), 3./5*U_Bagnold(1,1,tantheta));
 /** affectation */
    /* ------------------------------------------------------------------------*/
    x= (double*)calloc(nx+1,sizeof(double));
    h= (double*)calloc(nx+1,sizeof(double));
    hn=(double*)calloc(nx+1,sizeof(double));
    Q= (double*)calloc(nx+1,sizeof(double));
    Q2=(double*)calloc(nx+1,sizeof(double));
    Z= (double*)calloc(nx+1,sizeof(double));
    dZ=(double*)calloc(nx+1,sizeof(double));
    u= (double*)calloc(nx+1,sizeof(double));
    un=(double*)calloc(nx+1,sizeof(double));
    Fp=(double*)calloc(nx+1,sizeof(double));
    c0=(double*)calloc(nx+1,sizeof(double));
    hd=(double*)calloc(nx+1,sizeof(double));
    hg=(double*)calloc(nx+1,sizeof(double));
    hid=(double*)calloc(nx+1,sizeof(double));
    hig=(double*)calloc(nx+1,sizeof(double));
    ud=(double*)calloc(nx+1,sizeof(double));
    ug=(double*)calloc(nx+1,sizeof(double));
    Gamma=(double*)calloc(nx+1,sizeof(double));
    Y=(double*)calloc(nx+1,sizeof(double));
    
    
    l_alpha= (double*)calloc(N+1,sizeof(double));
    /**
     We consider that the flow domain is divided in the vertical direction into N layers
     of thickness $h_\alpha$ with $N + 1$   interfaces $z_{\alpha+1/2}$
     */
    h_alpha=alloc_2d_double(nx+1,N+1);
    u_alpha=alloc_2d_double(nx+1,N+1);
    Q_alpha=alloc_2d_double(nx+1,N+1);
    hn_alpha=alloc_2d_double(nx+1,N+1);
    un_alpha=alloc_2d_double(nx+1,N+1);
    hd_a=alloc_2d_double(nx+1,N+1);
    hg_a=alloc_2d_double(nx+1,N+1);
    hid_a=alloc_2d_double(nx+1,N+1);
    hig_a=alloc_2d_double(nx+1,N+1);
    ud_a=alloc_2d_double(nx+1,N+1);
    ug_a=alloc_2d_double(nx+1,N+1);
    fp=alloc_2d_double(nx+1,N+1);
    fd=alloc_2d_double(nx+1,N+1);
    u_alphap12=alloc_2d_double(nx+1,N+1);
    w_alphap12=alloc_2d_double(nx+1,N+1);
    z_alphap12=alloc_2d_double(nx+1,N+1);
    G_alphap12=alloc_2d_double(nx+1,N+1);
    /** initialisation cond init ----------------------------   */
    for(int i=0;i<=nx;i++)
    {
        x[i]=i*dx;
        Z[i]=tantheta*x[i];
        h[i]=(x[i]<2) +dry;
        u[i]=(x[i]<2)*0;
        Q[i]=u[i]*h[i];
    }
    for(int i=0;i<nx;i++)
    {  dZ[i]= Z[i+1]-Z[i];    }
    dZ[nx]=0;
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
     We consider the average velocities $u_\alpha$, for  $\alpha = 1,...,N$
     */
    for(int alpha=0;alpha<=N;alpha++)
    {
        for(int i=0;i<=nx;i++)
        {
            u_alpha[i][alpha]=u[i] ;
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
     The definition of $G$ is equ (2.13):
     $$G_{\alpha+1/2}=\partial_t z_{\alpha+1/2} + u_{\alpha+1/2} \partial_t z_{\alpha+1/2}-w(z,z_{\alpha+1/2},t)$$
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
         $z_{alpha+1/2} = Z + \Sigma_{j=1}^\alpha h_j$ (2.7)
         */
        for(int i=0;i<=nx;i++)
        {
            z_alphap12[i][0] = Z[i];
            for(int alpha=1;alpha<=N;alpha++){
                for(int j=1;j<=alpha;j++)
                    z_alphap12[i][j]=z_alphap12[i][j-1]+h_alpha[i][j];
            }
        }
        
        // variables pour la reconstruction
        for(int i=0;i<=nx;i++)
        {for(int alpha=1;alpha<=N;alpha++){
            hd_a[i][alpha]=h_alpha[i][alpha];
            ud_a[i][alpha]=u_alpha[i][alpha];
            hg_a[i][alpha]=h_alpha[i][alpha];
            ug_a[i][alpha]=u_alpha[i][alpha];}
            hd[i]=h[i];
            ud[i]=u[i];
            hg[i]=h[i];
            ug[i]=u[i];
            
        }
        /*
         c-----la reconstruction "hydrostatique"
         c-----construction d une section a gauche a l interface A_G(i-1/2) notee Ad(i-1)
         c-----construction d une section a droite a l interface A_D(i-1/2) notee Ag(i)
         call reconsetat(Ar(i-1),Al(i),delRaA0(i-1),Ad(i-1),Ag(i))
         */
        for(int i=1;i<=nx;i++)
        {
            double f1,f2;
            reconsetat(hd[i-1],hg[i],dZ[i-1],&f1,&f2);
            hid[i-1]=f1;
            hig[i]=f2;
            for(int alpha=1;alpha<=N;alpha++){
                //reconsetat(hd_a[i-1][alpha],hg_a[i][alpha],dZ[i-1],&f1,&f2);
                hid_a[i-1][alpha]=(1./N)*hid[i-1];
                hig_a[i][alpha]  =(1./N)*hig[i];
            }
        }
        /**
         Estimating a global wave velocity
         */
        for(int i=1;i<=nx;i++){
            c0[i]=C(ud[i-1],ug[i],hid[i-1],hig[i],1.);
        }
        /**
         flux corresponding to mass conservation accross the full layer $F_{p }= \Sigma_\alpha h_\alpha u_\alpha = Q$
         */
        for(int i=1;i<=nx;i++)
            Fp[i]=FR1(ud[i-1],ug[i],hid[i-1],hig[i],1.,c0[i]);
        for(int i=1;i<nx;i++)
            hn[i]=h[i]-dt*(Fp[i+1]-Fp[i])/dx;   //conservation de la masse
        /**
         Flux of mass $f_{p\alpha}$ in each layer (2.11): $f_{p\alpha}= h_\alpha u_\alpha$
         */
        for(int alpha=1;alpha<=N;alpha++)
            for(int i=1;i<=nx;i++)
                fp[i][alpha]=FR1(ud_a[i-1][alpha],ug_a[i][alpha],hid_a[i-1][alpha],hig_a[i][alpha],1/l_alpha[alpha],c0[i]);
        /**
         Flux of  of momentum
         $f_{d\alpha}= h_\alpha u_\alpha^2 + \frac{1}{l_\alpha} \frac{g h_\alpha^2}{2}$
         eq.(2.21) has a $l_\alpha$ in flux $\partial_x(h_\alpha u_\alpha^2 + \frac{1}{l_\alpha} \frac{g h_\alpha^2}{2})$
         */
        for(int alpha=1;alpha<=N;alpha++)
            for(int i=1;i<=nx;i++)
                fd[i][alpha]=FR2(ud_a[i-1][alpha],ug_a[i][alpha],hid_a[i-1][alpha],hig_a[i][alpha],1/l_alpha[alpha],c0[i]);
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
            }
        }
        /**
         Final inviscid update
         $q_\alpha^{n+1}$
         $=q_\alpha^{n} -$
         $\Delta t(\partial_x (f_{d \alpha}) + u_{\alpha+1/2} G_{\alpha+1/2} - u_{\alpha-1/2} G_{\alpha-1/2})$
         */
        for(int alpha=1;alpha<=N;alpha++)
        {
            double q;
            for(int i=1;i<nx;i++)
            {
                hn_alpha[i][alpha]=l_alpha[alpha]*hn[i];
                if(hn_alpha[i][alpha]>dry2){                             //conservation qunatité de mouvement
                    q=h_alpha[i][alpha]*u_alpha[i][alpha]
                    -dt*(fd[i+1][alpha]-fd[i][alpha])/dx
                    -dt*1./l_alpha[alpha]*( hig_a[i][alpha]*hig_a[i][alpha]/2 - hg_a[i][alpha]*hg_a[i][alpha]/2
                                           + hd_a[i][alpha]*hd_a[i][alpha]/2- hid_a[i][alpha]*hid_a[i][alpha]/2)/dx
                    +dt*(u_alphap12[i][alpha]*G_alphap12[i][alpha]-u_alphap12[i][alpha-1]*G_alphap12[i][alpha-1]);
                    un_alpha[i][alpha]=q/hn_alpha[i][alpha];}
                else{
                    un_alpha[i][alpha]=0.;}
            }
        }
        /**
         Simple Neumann BC
         */
        double hi= 1;//hi=1.2;//,hf=.01,U0=0;
        // hf= .02;
        /*
         FILE *g = fopen("F.IN", "r");
         fscanf(g,"hi=%lf          \n",&hi);
         fscanf(g,"hf=%lf          \n",&hf);
         fscanf(g,"U0=%lf          \n",&U0);
         fscanf(g,"nu=%lf          \n",&nu);
         fclose(g);
         */
        double zz=0,qq=0;
        for(int alpha=1;alpha<=N;alpha++)
        {
            // flux nul en entree sortie
            hn_alpha[0][alpha]= l_alpha[alpha]*hi;//
            zz= zz+hn_alpha[0][alpha];
            //hn_alpha[0][alpha]= hn_alpha[1][alpha];
            un_alpha[0][alpha]=0;
            //un_alpha[0][alpha]=  un_alpha[1][alpha];
            // un_alpha[0][alpha]=  3./5*U_Bagnold(1,1,tantheta)   ;//U_Bagnold((alpha-.5)/N*hi, hi, tantheta);
            qq=qq+un_alpha[0][alpha]*hn_alpha[0][alpha];
            hn_alpha[nx][alpha]= hn_alpha[nx-1][alpha];
            un_alpha[nx][alpha]= un_alpha[nx-1][alpha];
        }
        hn[0]=hn[1];
        hn[0]=hi;
        un[0]=qq/hn[0];
        hn[nx]=hn[nx-1];
        un[nx]=un[nx-1];
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
        double *dz    = (double*)calloc(N+1,sizeof(double));
        double *mu    = (double*)calloc(N+1,sizeof(double));
        double *mu12  = (double*)calloc(N+1,sizeof(double));
        double *dudz  = (double*)calloc(N+1,sizeof(double));
        double *a     = (double*)calloc(N+1,sizeof(double));
        // static double NU;
        /*
         compute shear and viscosity
         */
        
        for(int i=1;i<nx;i++)
        {
            //hn[i]=1;
            //    for (int j = 0; j <= N; ++j){
            //        u_alpha[i][j]=1;
            //        un_alpha[i][j]=1;}
        }
        
        for(int i=1;i<nx;i++)
        {
            if(h[i]>dry2){
                
                for (int j = 0; j <= N; ++j){
                    dz[j] = h[i]/N;
                }
                
                for (int j = 1; j < N; ++j)
                    dudz[j] =(u_alpha[i][j+1]-u_alpha[i][j])/dz[j];
                dudz[N] = 0;
                
                double signe,In,muI,pz;
                for (int j = 1; j < N; ++j){
                    //pz = hn[i]-((z_alphap12[i][j]+z_alphap12[i][j+1])/2.-Z[i]) ; //
                    //pz= hn[i]*((1.-(j-0.5)/N)/2. + (1.-(j+0.5)/N)/2.);
                    signe=1;
                    pz= h[i]*(1.-(j-0.5)/N);
                    pz= h[i]*(1.-(j*1.)/N);
                    pz= h[i]-(z_alphap12[i][j]-Z[i]);
                    In = 1./32*fabs(dudz[j])/sqrt(pz+.0000001);
                    muI = (mus + Deltamu*In/(I0+In));
                    tauy=muI*pz;
                    mu12[j] = tauy*signe/fabs(dudz[j]+.000001);              
                    if(mu12[j]>1000) Y[i]=j*1./N; // Y[] is in fact use less
                    mu12[j]= fmin(mu12[j], 1000.) ;
                }
                mu12[0] =mu12[1];
                mu12[N] =0;
                
                for (int j = 0; j <= N; ++j){
                    mu[j]=(mu12[j]+mu12[j])/2;}
                //mu[N]=mu12[N];
                //mu[N-1]=0*mu12[N];
                
                /**
                 
                 $h_j (u_j^n - u_j^o)/\Delta t = \tau_{j+1/2}-\tau_{j-1/2}$
                 
                 $\tau_{j+1/2}=$
                 
                 Consturtion of the tridiagonal matrix corresponding to the "heat equation"
                 
                 $(u_j^n - u_j^o)/\Delta t =   \frac{\nu }{\Delta z^2}  (u_{j-1}^n-2 u_j^n +u_{j+1}^n)$
                 so that we have to solve
                 $u_j^n +  \Nu  (-u_{j-1}^n+2 u_j^n -u_{j+1}^n = u_j^o$
                 
                 with $\Nu = \frac{\nu \Delta t}{\Delta z^2}$
                 */
                
                for(int j=1; j<N; ++j)
                    a[j]= dt*( 2* mu[j] + 0*mu[j+1])/(dz[j]*(dz[j] + dz[j+1]));
                //a[N] = 0*dt*mu[N]/(dz[N]*dz[N]);
                //  a[N-1] = dt*(  mu12[N-1]*2 +  0*mu[N])/(dz[N-1]*(dz[N-1] + dz[N]));
                
                double a1 = dt*mu[1]/(dz[1]*dz[1]);
                
                for(int j=2; j<N; ++j){
                    diags[j] =  -a[j];
                    diagp[j] = 1 + a[j] + a[j-1] ;
                    diagi[j] =  -a[j-1];
                    rhs[j] =    un_alpha[i][j];
                }
                /**
                 values at the wall
                 $u_1^n-u_1^o = \frac{\nu \Delta t}{\Delta z^2}  (u_{-1}^n-2 u_1^n +u_2^n)$
                 with $u_1^n=-u_{-1}^n$ so that the velocity is 0 at the wall
                 */
                diagp[1]= 1 + a[1] + 2*a1;
                diags[1]= -a[1];
                rhs[1]= un_alpha[i][1];
                /**
                 value at the surface,
                 $u_N^n-u_N^o = \frac{\nu \Delta t}{\Delta z^2}  (u_{N-1}^n-2 u_N^n +u_{N+1}^n)$ with
                 here   $u_N^n=-u_{N+1}^n$ so that the derivatice is 0 at the top
                 */
                
                diagp[N] = 1 + a[N-1];
                diagi[N]= -a[N-1];
                rhs[N] = un_alpha[i][N];
                
                
                /*     inversion de matrice tridiag                               */
                /*     ai ui-1 + bi ui +ci ui+1 = rhs i                           */
                /*     diagi[j] u[j-1] +  diagp[j] u[j] + diags[j] u[j+1]= rhs[j] */
                /*        u[j] = a[j] * u[j + 1] + b[j];                            */
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
                
            }else{
                for (int j = 1; j <= N; ++j)
                    un_alpha[i][j]=0;}
            
            if((it%100==0)&&(i==nx/4)){
                FILE *g1= fopen("profil1.OUT", "w");
                for (int j = 1; j <= N; ++j) {
                    // double s;
                    // s=0;
                    fprintf(g1," %lf %lf  %lf %lf %lf %lf\n",
                            z_alphap12[i][j]-Z[i],un_alpha[i][j],
                            hn[i]*(1-(j-.0)/N), dudz[j-1],
                            (mu[j-1])/2*dudz[j-1],mu12[j]);
                }
                fclose(g1);
                fprintf(stderr, " %d \n",i);
                // getchar();
            }
        }
        free(uN);
        free(a);
        free(diags);free(diagp);free(diagi);free(rhs);
#endif
        
        //final swap
        for(int i=0;i<=nx;i++)
        {
            h[i]=hn[i];
            Q[i]=0;
            Q2[i]=0;
            Gamma[i]=0;
            for(int alpha=1;alpha<=N;alpha++)
            {
                h_alpha[i][alpha]=hn_alpha[i][alpha];
                u_alpha[i][alpha]=un_alpha[i][alpha];
                Q[i]=Q[i]+hn_alpha[i][alpha]*un_alpha[i][alpha];
                Q2[i]=Q2[i]+hn_alpha[i][alpha]*un_alpha[i][alpha]*un_alpha[i][alpha];
            }
            if(hn[i]>dry2){
                un[i]=Q[i]/hn[i];
                Gamma[i]=(Q[i]>0? hn[i]*Q2[i]/Q[i]/Q[i] : 0.0);
            } else {un[i]=0;
                Gamma[i]=0;}
            for(int alpha=1;alpha<=N;alpha++)
                u_alpha[i][alpha]=un_alpha[i][alpha];
            u[i]=un[i];
        }
        
        
        /**
## Results
         */
        if(it%100==0){
            /* Saving the fields */
            //  if(i==nx/2){
            
            FILE *gq= fopen("profil0.OUT", "w");
            for (int j = 1; j <= N; ++j) {
                fprintf(gq," %lf %lf %lf %lf %lf %lf \n",
                        (j-.5)*hn[nx/2]/N,un_alpha[nx/2][j],w_alphap12[nx/2][j],t,z_alphap12[nx/2][j]-Z[nx/2],
                        hn[nx/2]-(z_alphap12[nx/2][j]-Z[nx/2]));
            }
            fclose(gq);
            
            FILE *gn= fopen("profilN.OUT", "w");
            for (int i=0; i<=nx/5;i++){
                for (int j = 1; j <= N; ++j) {
                    fprintf(gn," %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
                            (j-.5)*hn[5*i]/N,un_alpha[5*i][j],w_alphap12[5*i][j],t,z_alphap12[5*i][j]-Z[5*i],
                            hn[5*i]-(z_alphap12[5*i][j]-Z[5*i]),x[5*i],un_alpha[5*i][N],hn[5*i]);
                }
                fprintf(gn,"\n");
            }
            fclose(gn);
            
            // }
            
            
            FILE *g= fopen("solxhQt.OUT", "w");
            for (int i=0; i<=nx;i++)
            {      fprintf(g,"%lf %lf %lf %lf %lf %lf  \n",x[i],h[i],Q[i],t,Gamma[i],Y[i]);}
            fprintf(g,"\n");
            fclose(g);
#ifdef gnuX
            printf("set key left;t=%lf\n",t);
            printf("set arrow nohead  from 0,0 to 0.2,0 lw 1 \n");
            printf("set arrow nohead  from 0,0 to 0.2,0.2*1.8138 lw 1 \n");
            printf("set xlabel \"x\" ; set title \"t=%lf \"\n",t);
            double y1=-1,y2=2,scu=5./100*x[nx],flat=1;;
            printf("p[%lf:%lf][%lf:%lf]  '-' u 1:2 t'Q'  w l linec 2,'-'  t'h' w l linec 3,'-'  t'h_1' w l,'-'  t'Z' w l linec -1,'-' not w l linec 1,'-' not w l linec 1,'-' w l not linec 1,'-' w l t 'Gamma', 1,5./4\n ",
                   0*x[0],x[nx],y1,y2);
            for (int i=0; i<=nx;i++)
            {
                printf("%lf %lf \n",x[i],Q[i]);}
            printf("e \n");
            for (int i=0; i<=nx;i++)
            {
                printf("%lf %lf \n",x[i],h[i]+(1.0-flat)*Z[i]);}
            printf("e \n");
            for (int i=0; i<=nx;i++)
            {
                printf("%lf %lf \n",x[i],h_alpha[i][1]+(1.0-flat)*Z[i]);}
            printf("e \n");
            for (int i=0; i<=nx;i++)
            {
                printf("%lf %lf \n",x[i],Z[i]);}
            printf("e \n");
            
            // vitesses
            for (int j=1; j<=N;j++)
            {
                printf("%lf %lf \n",x[nx/4]+scu*u_alpha[nx/4][j],z_alphap12[nx/4][j]-(flat)*Z[nx/4]);}
            printf("%lf %lf \n",x[nx/4] ,z_alphap12[nx/4][N]-(flat)*Z[nx/4]);
            printf("%lf %lf \n",x[nx/4],0.);
            printf("e \n");
            for (int j=1; j<=N;j++)
            {
                printf("%lf %lf \n",x[nx/2]+scu*u_alpha[nx/2][j],z_alphap12[nx/2][j]-(flat)*Z[nx/2]);}
            printf("%lf %lf \n",x[nx/2] ,z_alphap12[nx/2][N]-(flat)*Z[nx/2]);
            printf("%lf %lf \n",x[nx/2],0. );
            printf("e \n");
            
            for (int j=1; j<=N;j++)
            {
                printf("%lf %lf \n",x[3*nx/4]+scu*u_alpha[3*nx/4][j],z_alphap12[3*nx/4][j]-(flat)*Z[3*nx/4]);}
            printf("%lf %lf \n",x[3*nx/4] ,z_alphap12[3*nx/4][N]-(flat)*Z[3*nx/4]);
            printf("%lf %lf \n",x[3*nx/4],0.-0*(flat)*Z[3*nx/4]);
            printf("e \n");
            for (int i=0; i<=nx;i++)
            {
                printf("%lf %lf \n",x[i],Gamma[i]);}
            printf("e \n");
#endif
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
 
## Run
 
 Programme en C, flux Rusanov (essayer HLC) compilation
 on le compile et on crée l’exécutable db
 
~~~bash
 cc -O3 -ffast-math -std=c99 -lm svdb_mumult_collapse.c -DgnuX=1 -o svdb_mumult_collapse
~~~
 
 pour lancer le programme:
 
~~~bash
 ./svdb_mumult_collapse | gnuplot
~~~
 
## Results
 
 
 profil de vitesse avec z u v  et t est créé, pour le tracer, on lance gnuplot:
 
 
~~~gnuplot
 set ylabel "z"; set xlabel "u";
 p 'profil0.OUT' u 1:2 not w l linec 3
~~~
 
 un fichier appelé solxhQt.OUT avec x h Q et t est créé, pour le tracer, on lance gnuplot:
 
 Plot of $Q$ and $h$
~~~gnuplot
 set ylabel "x";
 p 'solxhQt.OUT' u 1:2 t'h'  w l,''u 1:3 t 'Q' linec 3
~~~

 Plot of $Q$ and $h$, with $\Gamma$ compared with 5./4 the Bagnold value, + $Q/h$ the mean velocity.
~~~gnuplot
 set ylabel "x";
 p[:][0:2.5] 'solxhQt.OUT' u 1:2 t'h'  w l,''u 1:3 t 'Q' linec 3,''u 1:5 t 'Gamma' linec 3,5./4,'' u 1:($3/$2) t 'u'
~~~
 
velocity profiles $u/U_{bagnold}$ near the nose
 
~~~gnuplot
 set ylabel "z";
 p 'profilN.OUT'u ($1/$9):(($9<.2&& $9>.05?$2:NaN)/$8) w lp
~~~
 
 
 ~~~gnuplot
  set ylabel "x";
  p[:][0:2.5] 'solxhQt.OUT' u 1:2 t'h'  w l,''u 1:3 t 'Q' linec 3,''u 1:5 t 'Gamma' linec 3,5./4,'' u 1:($3/$2) t 'u','' u 1:6  t'Y'
 ~~~
 
## Bibliography
 
* Emmanuel Audusse, Marie-Odile Bristeau, Benoıt Perthame, and Jacques Sainte-Marie
 A multilayer saint-venant system with mass exchanges for shallow water flows. derivation and numerical validation
 ESAIM: M2AN 45 (2011) 169–200 DOI: 10.1051/m2an/2010036
 
* NR [http://www.nrbook.com/a/bookcpdf.php]()

* Pierre-Yves Lagrée et al, "Granular front for flow down a rough incline: about the value of the shape factor in depths averaged models" [http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/epjconf162335.pdf]()




 
 
 
 PYL version 1, Manchester, fev 2016 v2 07/2019
*/









