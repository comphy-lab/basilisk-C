/**
# Propagation d'onde de marée dans une embouchure


## problème

"La déformation des ondes de marées en pleine mer est assez bien connue. En grand profondeur, cette dynamique peut raisonnablement être modélisée par les équations linéarisées de Saint-Venant.
 Cependant, lors de son entrée sur le plateau continental et encore plus à l'intérieur des estuaires, l'onde se propage en se déformant."
Il s'agit du projet proposé par [wiki.shf-hydro.fr](https://wiki.shf-hydro.fr/index.php?title=Tidal_wave_propagation_in_estuary_/_propagation_d%27une_onde_de_marée_en_estuaire).
"A partir de l'ensemble de ces informations, il s'agira de proposer un modèle simplifié d'estuaire.
Il s'agit ensuite de mettre clairement en évidence les effets spécifiques de chaque terme dans la déformation du signal temporel de marée.

## Topographie sans dimension et paramètres

La carte [https://www.geoportail.gouv.fr/donnees/carte-littorale]() nous donne une longueur $L=89 000 m$, avec 
profondeur $h_0=4 m$. 
A la louche les 4 stations sont réparties en 0 L/4 L/2 3L/4 et L: 1 (Pointe de la Grave), 2 (Richards), 3 (Lamena), 4 (Pauillac), 4 (Bordeaux):

La vitesse de propagation des ondes est donc  
$c_0=\sqrt(9.81*4)= 6.2 m/s$, dans cette version simplifiée, on se donne une seule onde de marée, la M2. La pulsation période 
$T=(12*60*60+25*60+12)$. Sans dimension la hauteur d'entrée est donc, aved $A$ l'amplitude relative de la marée mesurée en $h_0$:
$$h = 1 + A \cos(2 \pi t).$$
Sur les marégraphes, est environ 0.5 (différence haut moins bas: amplitude de la marée 4 m, 
donc deux fois l'amplitude d'un cosinus, profondeur 4 m).


Le domaine fait alors à peu près 0.32 sans dimension car: $c_0 T= 380km$ soit et $L=89km$. 
On suppose que tout sort et continue en amont après Bordeaux.

Au premier ordre en linéarisé on a simplement 
$$\frac{\partial u}{\partial t} = -\frac{\partial h}{\partial x},\;\;\;\;\;\;\frac{\partial h}{\partial t} = -\frac{\partial u}{\partial x}$$
la solution de ce problème est  $h=1 + A \cos (2 \pi (t -x))$ et  $u= A \cos (2 \pi (t -x)).$


## Numerics

An example in simple C (not a Basilisk example). We solve here the  problem with a simple scheme with Rusanov flux. 
$$\dfrac{U_i^{n+1}-U_i^{n}}{\Delta t}+\dfrac{F^n_{i+1/2}-F^n_{i-1/2}}{\Delta x}=0,$$

*/
#include <stdio.h>
#include <stdlib.h> 
#include <math.h> 
#include <string.h>

  double*x=NULL,*h=NULL,*u=NULL,*Q=NULL;
  double t,dt,tmax,dx,Z; 
  int nx;
/** definition of the fluxes, here a simple Rusanov
$$
F_i=f(U_{i-1},U_i)=
\dfrac{F(U_{i-1})+F(U_i)}{2}-c\dfrac{U_{i}-U_{i-1}}{2},
$$
which is coded as
$$f(U_G,U_D)=
\dfrac{F(U_G)+F(U_D)}{2}-c\dfrac{U_D-U_G}{2}
$$
with
$c=\sup\limits_{U=U_G,U_D}({\sup\limits_{j\in\{1,2\}}} |\lambda_j(U)|),$
where $\lambda_1(U)$ and $\lambda_2(U)$ are the eigen values.
*/
double FR1(double ug,double ud,double hg,double hd)
 { double c;
   c=fmax(fabs(ug)+sqrt(hg),fabs(ud)+sqrt(hd));
   return (hg*ug+hd*ud)*0.5-c*(hd-hg)*0.5;
 }
double FR2(double ug,double ud,double hg,double hd)
 { double c;
   c=fmax(fabs(ug)+sqrt(hg),fabs(ud)+sqrt(hd));
   return (ug*ug*hg + hg*hg/2. + ud*ud*hd + hd*hd/2.)*0.5 - c*(hd*ud-hg*ug)*0.5;
 }
  
/*     -------------------------------------------------    */        
int main (int argc, const char *argv[]) {
    int  i,it=0;
    FILE *g;
    double*fp=NULL,*fd=NULL,*un=NULL,*hn=NULL;
    double  q;
    
// parametres --------------------------------------------------------------  
dt=0.001;
tmax=5;
dx=0.0025;
nx=128;
t=0;    

/**
domaine de longueur 2

*/
  
  fprintf(stderr,"  ---------\\                /        \n"); 
  fprintf(stderr,"            \\             /  -->    \n");   
  fprintf(stderr,"               -----------            \n");   
  fprintf(stderr,"  ____________________________________\n");   
/* ----------------------------------------------------------------------*/      
      x= (double*)calloc(nx+1,sizeof(double));
      h= (double*)calloc(nx+1,sizeof(double));
      Q= (double*)calloc(nx+1,sizeof(double));
      u= (double*)calloc(nx+1,sizeof(double));
      fp=(double*)calloc(nx+1,sizeof(double));
      fd=(double*)calloc(nx+1,sizeof(double)); 
      un=(double*)calloc(nx+1,sizeof(double));
      hn=(double*)calloc(nx+1,sizeof(double));
      
 fprintf(stderr,"tmax= %lf \n",tmax);
      
// initialisation cond init ----------------------------      
  for(i=0;i<=nx;i++)
    {  x[i]=(i-nx/2)*dx; 
     Z=0;
     h[i]=1;
     u[i]=0;
     Q[i]=u[i]*h[i];
     }      
// initialisation du fichier de sortie ----------------------      
      g = fopen("solxhQt.OUT", "w");
      fclose(g); 
                  
    while(t<tmax){   // boucle en temps
      t=t+dt;
      it=it+1;      
 
   for(i=1;i<=nx;i++)
    {    
        fp[i]=FR1(u[i-1],u[i],h[i-1],h[i]);
        fd[i]=FR2(u[i-1],u[i],h[i-1],h[i]);  
    }
   
   for(i=1;i<nx;i++)
    {
      hn[i]=h[i]- dt*(fp[i+1]-fp[i])/dx;   //conservation de la masse
  if(h[i]>0.){                             //conservation quantité de mouvement
      q=h[i]*u[i]-dt*(fd[i+1]-fd[i])/dx ;  
      un[i]=q/hn[i];}
      else{
      un[i]=0.;}
    }
/**
   on ajoute un frottement de Manning donc la vitesse est à corriger de :

   $$ \frac{\partial u}{\partial t} =- n^2\frac{c_0gT}{ h_0^{4/3}}(\frac{ u^2}{ h^{4/3}})$$
   A.N.: 
$cf=sqrt(9.81*4)*9.81*(12*60*60+25*60+12)/(4**(4./3))$

 donc pour fond avec n=0.02 (graviers), $cf=173$ 
 (cf [table](http://www.fsl.orst.edu/geowater/FX3/help/8_Hydraulic_Reference/Mannings_n_Tables.htm))

*/
    
   double cf =173;
   for(i=1;i<nx;i++) { if(h[i]>0.) un[i]= un[i]/(1+cf*dt*fabs(un[i])/pow(hn[i],4./3)); } 
/**

 On diminue l'amplitude A car il n'y a aucune dissipation 

 on écrit que les ondes entrent et ne sortent pas $u +2 c =cst$ donc $u + 2 sqrt(h)$ est constant. 

en vives eaux 2.5m d'amplitude, (5 m de la basse à la haute)) A=2.5/4
en mortes eaux 1.5 (3 m de bas à haut)
*/ 
    double A=1.5/4,omega=2*3.141516;
    hn[0]=-hn[1] + (1+ A*cos(omega*t)) *2 ;
    un[0]=un[1]-2*(sqrt(hn[1])-sqrt(hn[0]));  
    hn[nx]=hn[nx-1];
    un[nx]=un[nx-1];  
   
   //swap      
  for(i=0;i<=nx;i++)
    { 
      h[i]=hn[i];
      u[i]=un[i];
      Q[i]=hn[i]*un[i];
    }



  if(it%10==0){
    /* Saving the fields */ 
  g = fopen("solxhQt.OUT", "a");
  for (i=0; i<=nx;i++)
  {      fprintf(g,"%lf %lf %lf %lf \n",x[i],h[i],Q[i],t);}
   fprintf(g,"\n");
   fclose(g);
  }
#ifdef gnuX
if(it%20==0){  
  printf("t=%lf\n",t);
  printf("set xlabel \"x\" ; set title \"t=%lf \"\n",t);
     y1=-.25,y2=1.5;   
  printf("h(x)=1+%lf*cos(%lf*(t-x -%lf)) \n",A,omega,nx*dx/2);

   
  printf("p[%lf:%lf][%lf:%lf] h(x) t'h exact','-' u 1:2 t'Q'  w l,'-'  t'h' w lp,'-'  t'Z' w l linec -1\n ",
         -nx*dx/2,nx*dx/2,y1,y2); 
    for (i=0; i<=nx;i++)
    {
    printf("%lf %lf \n",x[i],Q[i]);}
    printf("e \n");
    for (i=0; i<=nx;i++)
    {
    printf("%lf %lf \n",x[i],h[i]+Z);}
    printf("e \n");
    for (i=0; i<=nx;i++)
    {
      printf("%lf %lf \n",x[i],Z);}
    printf("e \n");}
#endif   
      }
      
     free(x);
     free(fp);
     free(fd); 
     free(un);
     free(hn);  
    return 0;
}


/**

## Run

Programme en C simple fichier (en cours on a regardé un fichier plus compliqué). Cas de rupture de barrage flux Rusanov (essayer HLC) compilation
on le compile et on crée l’exécutable `embouch`

~~~bash
cc -O3 -ffast-math -std=c99 -lm embouch.c -o embouch
~~~

pour lancer le programme:

~~~bash
 ./embouch 
~~~ 

Alternatively, we can compile with `gnuX` option which allows to pipe in gnuplot

~~~bash
 cc  -O3 -ffast-math -std=c99 -lm  -DgnuX  embouch.c -o embouch    
~~~

to run the code 

~~~bash
 ./embouch | gnuplot -persist
~~~ 



## Results


un fichier appelé `solxhQt.OUT` avec x h Q et t est créé, pour le tracer, on lance gnuplot:
 
 Tracé à $t=nT+T/5$, $nT + .1 T$, $nT + .2T$, $nT + .3T$,  de la hauteur d'eau en fonction de $x$:
   
~~~gnuplot  
 set ylabel "h"; set xlabel "x";  
  p[:][0:]'solxhQt.OUT' u 1:($4==4.?$2:NaN) t 't=0.0',''u 1:($4==4.1?$2:NaN)t't=0.1',''u 1:($4==4.2?$2:NaN)t't=0.2',''u 1:($4==4.3?$2:NaN)t't=0.3'
 # p[ : ][0:]'solxhQt.OUT' u 1:($4==4.2?$2:NaN) not,''u 1:($4==4.4?$2:NaN) not,'' u 1:($4==4.6?$2:NaN) not,'' u 1:($4==4.8?$2:NaN) not,'' u 1:($4==5?$2:NaN) not
~~~ 

Ce qui donne pour 1 (Pointe de la Grave), 2 (Richards), 3 (Lamena), 4 (Pauillac), 4 (Bordeaux):

~~~gnuplot     
 set ylabel "h"; set xlabel "t";  
 p[2.5:4.5][0:2]'solxhQt.OUT' u 4:($1==-.16?$2:NaN) t'1',''u 4:($1==-.08?$2:NaN) t'2',\
  '' u 4:($1==0?$2:NaN) t'3',''u 4:($1==.08?$2:NaN) t'4',''u 4:($1==.16?$2:NaN) t'5'
~~~ 


# Go further

 add a slope...
    
   
Do it with Basilsik, use the real ... topography

 
# Links

 
 * [wiki.shf-hydro.fr](https://wiki.shf-hydro.fr/index.php?title=Tidal_wave_propagation_in_estuary_/_propagation_d%27une_onde_de_marée_en_estuaire) 
  
 
# Biblio
 
* [http://www.aquitaineonline.com/actualites-en-aquitaine/nature-et-environnement/1457-marees-estuaire-gironde.html]() 
* [https://www.geoportail.gouv.fr/donnees/carte-littorale]()
* [cours PYL](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv)
* [cours PYL](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf) sur numérique
 

 version sept 2018 Abidjan
 PYL
*/