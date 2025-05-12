/**
# Rupture de barrage/ Dam Break


An example in simple C, this is not a `Basilisk` example, 
this is the example of the [lecture](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf). 

We solve here the dambreak problem with a simple scheme with Rusanov flux.
$$\dfrac{U_i^{n+1}-U_i^{n}}{\Delta t}+\dfrac{F^n_{i+1/2}-F^n_{i-1/2}}{\Delta x}=0,$$

*/

/** definition of the fluxes, here a simple Rusanov
$$
F_{i-1}=f(U_{i-1},U_i)=
\dfrac{F(U_{i-1})+F(U_i)}{2}-c\dfrac{U_{i}-U_{i-1}}{2},
$$
which is coded as
$$f(U_G,U_D)=
\dfrac{F(U_G)+F(U_D)}{2}-c\dfrac{U_D-U_G}{2}
$$
with
$c=\sup\limits_{U=U_G,U_D}({\sup\limits_{j\in\{1,2\}}} |\lambda_j(U)|),$
where $\lambda_1(U)$ and $\lambda_2(U)$ are the eigen values.
 
 
 
 A note on the number of points:  `D` refeers to Domain, we compute the values for the `nx` points
 in the domain, but we need the index 0 and nx+1. Those two extra localisations are
 the ghost cells, `G` refers to Ghost.
 
~~~bash
 G   D   D   D   D    D    G
 |   |   |   |   |    |   |
 0 | 1 | 2 | 3 |...| nx| nx+1
 --|---|---|---|---|---|
    x1  x2             L
 <-------------------->
 <-->
 
 nx cells, +2 ghost cells, (nx+2) data
 
 first cell, index 1,  between  `X0` and `X0+Delta`, centred in  `Delta/2`
 ith cell beween `(i-1) Delta` (left) and `i Delta`(right) centered in `(i-1/2)Delta`
 
 `Delta=L0/N`
~~~


A note on the flux 

~~~bash 
 |           |
 |->Fi       |->Fi+1
 |  thetai   |
 xi+Delta/2
 xi   xi+Delta/2
~~~
 
 Fi is the flux across interface between hi-1 and hi  (sometimes Fi-1/2)
 
 Fi+1 is the flux across interface between hi and hi+1  (sometimes Fi+1/2)
 
 
 
 so, in practice for the height, 
 $$ h_i^{n+1}=  h_i^{n} + {\Delta t}\dfrac{F^n_{i+1/2}-F^n_{i-1/2}}{\Delta x},$$
 
 for `h[i]` the first component of $F_{i-1}$ is the `Fi`  of the figure and is 
the subroutine `FR1`, this gives

~~~bash
hn[i]=h[i]- dt*(fp[i+1]-fp[i])/dx;
~~~

for flow rate, we save the velocity and re compute the flow rate with depth:

~~~bash
q=h[i]*u[i]-dt*(fd[i+1]-fd[i])/dx ;  
      un[i]=q/hn[i];
~~~



# Program/ code 

you can extract this part from here
*/

#include <stdio.h>
#include <stdlib.h> 
#include <math.h> 
#include <string.h>

  double*x=NULL,*h=NULL,*u=NULL,*Q=NULL;
  double t,dt,tmax,dx,Z; 
  int nx;

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
dt=0.005;
tmax=10;
dx=0.125;
nx=320;
t=0;    
  
  fprintf(stderr,"  ------------------\\                          \n"); 
  fprintf(stderr,"                   <--\\                        \n");   
  fprintf(stderr,"                        \\--->                  \n");   
  fprintf(stderr,"  ________________________\\____________________\n");   
/* ----------------------------------------------------------------------*/      
      x= (double*)calloc(nx+2,sizeof(double));
      h= (double*)calloc(nx+2,sizeof(double));
      Q= (double*)calloc(nx+2,sizeof(double));
      u= (double*)calloc(nx+2,sizeof(double));
      fp=(double*)calloc(nx+2,sizeof(double));
      fd=(double*)calloc(nx+2,sizeof(double));
      un=(double*)calloc(nx+2,sizeof(double));
      hn=(double*)calloc(nx+2,sizeof(double));
      
 fprintf(stderr,"tmax= %lf \n",tmax);
      
// initialisation cond init ----------------------------      
  for(i=0;i<=nx+1;i++)
    {  x[i]=(i-nx/2)*dx;
   //dambreak
     Z=0;
     h[i]=1*( x[i]<0); 
     u[i]=0;
     Q[i]=u[i]*h[i];
     }
       
// initialisation du fichier de sortie ----------------------      
      g = fopen("solxhQt.OUT", "w");
      fclose(g); 
                  
    while(t<tmax){   // boucle en temps
      t=t+dt;
      it=it+1;      
 
   for(i=1;i<=nx+1;i++)
    {    
        fp[i]=FR1(u[i-1],u[i],h[i-1],h[i]);
        fd[i]=FR2(u[i-1],u[i],h[i-1],h[i]);  
    }
   
   for(i=1;i<nx+1;i++)
    {
      hn[i]=h[i]- dt*(fp[i+1]-fp[i])/dx;   //conservation de la masse
  if(h[i]>0.){                             //conservation qunatité de mouvement
      q=h[i]*u[i]-dt*(fd[i+1]-fd[i])/dx ;  
      un[i]=q/hn[i];}
      else{
      un[i]=0.;}
    }
      
 // put friction here      

      
 // boundary conditions     
 // flux nul en entree sortie 
    hn[0]=hn[1];
    un[0]=un[1];  
    hn[nx+1]=hn[nx];
    un[nx+1]=un[nx];
   
   //swap      
  for(i=0;i<=nx;i++)
    { 
      h[i]=hn[i];
      u[i]=un[i];
      Q[i]=hn[i]*un[i];
    }

  if(it%100==0){
// Saving the fields 
  g = fopen("solxhQt.OUT", "a");
  for (i=0; i<=nx;i++)
  {      fprintf(g,"%lf %lf %lf %lf \n",x[i],h[i],Q[i],t);}
   fprintf(g,"\n");
   fclose(g);
  }
#ifdef gnuX
if(it%10==0){
  double  y1=-.1,y2=1.25;
  printf("t=%lf\n",t);
  printf("set xlabel \"x\" ; set title \"t=%lf \"\n",t);
     y1=-.1,y2=1.25;  
  printf("set label \"+\" at 0,.296296296296296 \n"); 
  printf("h(x)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t)) \n");

   
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

this is the end of the program. 
You copy and paste in your favorite editor (an ASCII editor).
You save as `svdb.c`in a specific directory somewhere in your hard drive.



You open two terminals (`xterm`), one in wich you compile, the other in which you open``gnuplot`. then `gnuplot` will genrate a figure in a fourth windon (the `X11` window). 

You never close the editor, you never close the two terminals. You play and circle between them. 



## Run avec cc

Programme en C simple fichier (en cours on a regardé un fichier plus compliqué). Cas de rupture de barrage flux Rusanov (essayer HLC). On le compile et on crée l’exécutable `db`

~~~bash
cc -O3 -ffast-math -std=c99 -lm svdb.c -o db
~~~

pour lancer le programme:

~~~bash
 ./db 
~~~ 

Alternatively, we can compile with `gnuX` option which allows to pipe in gnuplot

~~~bash
cc -O3 -ffast-math -std=c99 -lm svdb.c -o db -DgnuX
~~~

to run the code 

~~~bash
 ./db | gnuplot -persist
~~~ 

## Run avec `make`
Avec le `makefile`
 
~~~bash
 make svdb.tst
 make svdb/plots
 make svdb.c.html
~~~

## Results


un fichier appelé solxhQt.OUT avec en première colonne `x`, en 
seconde `h`, en troisième  `Q` et en quatrième `t` est créé.
  Pour le tracer, on lance gnuplot:

~~~bash
 p 'solxhQt.OUT'  w l  
~~~

this gives a view i  with $h$ function of $x$ at the different $t$ superposed.
   
   
~~~gnuplot
 set ylabel "h"; set xlabel "x";  
 p 'solxhQt.OUT'  w l  
~~~



par défaut, 1a  première colonne est l'abcisse, et la deuxième l'ordonnée, on rajoute  `Q` troisième variable en fonction de la première

~~~bash
 reset; set xlabel "x";  
 p 'solxhQt.OUT'  w l,'' u 1:3 w l 
~~~


 ~~~gnuplot
reset; set xlabel "x";  
 p 'solxhQt.OUT'  w l,'' u 1:3 w l 
~~~




Si on veut une représntation 3D (h,x) avec le temps t vers l'arrière 

~~~bash
 set ylabel "t"; set xlabel "x"; set hidden3d;
 sp[-20:20][0:10][0:2] 'solxhQt.OUT' u 1:4:2 not w l linec 3
~~~

   this gives a view in 3D with $h$ function of $x$ at the different $t$.
   


    
       
~~~gnuplot     
  set ylabel "t"; set xlabel "x";  set hidden3d;
  splot  'solxhQt.OUT' u 1:4:2 not w l
~~~ 
 
 An other plot, comparison exact and computed at two times
 
~~~gnuplot     
  set ylabel "h"; set xlabel "x";  
  set label "+" at 0,.296296296296296 
  h(x,t)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t))
  p  'solxhQt.OUT' u 1:(($4==5)||($4==10)? $2 :NaN)w p, h(x,5),h(x,10)
~~~ 

for use of `gnuplot` see
[http://basilisk.fr/sandbox/M1EMN/BASIC/gnuplot_examples.c]()

# Go further

i) First you can change the initial boundary condition to mimick a dam break in a river

~~~bash
 h[i]=0.1+.9*( x[i]<0);
~~~

observe the formation of a shock


ii) test the simple tsunami model 

~~~bash
 h[i]=1+.2*exp(-x[i]*x[i]);
~~~
  and see the right/left waves...
 
  
iii) You can add some viscosity in the code to test the example of [viscous collapse](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv)
    
first change for a heap (two dams breaking!): 

~~~bash
h[i]=1*( fabs(x[i])<1);
~~~
                                                                                
then add viscosity  $\frac{\partial u}{\partial t} = - 3*u /h$ with implicit discretization 
$u(t+\Delta t,x)- u(t,x) = - 3 * u(t+\Delta t,x) \Delta t /h(x,t)^2$ hence 
$$ u(t+\Delta t,x)= \frac{u(t,x)}{1+ 3 \Delta t /h(t,x)^2}$$ 

~~~bash
for(i=1;i<nx+1;i++) { if(h[i]>0.) un[i]= un[i]/(1+3*dt/h[i]/h[i]); } 
~~~
   
 

   
   
   
check the selsimilar variable:

~~~bash
p[-5:5]'solxhQt.OUT' u ($1/$4**.2):($4>10?$2*$4**.2:0) w l
~~~

iv) You can put turbulent friction  (this is the Dressler Problem See  Chanson page 357 and see
  [dam break Dressler](http://basilisk.fr/sandbox/M1EMN/Exemples/damb_dressler.c) )

   turbulent  friction  
  $$
  \frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
  $$
  with $C_f=.05$. 
  
~~~bash                       
for(i=1;i<nx+1;i++) { if(h[i]>0.) un[i]= un[i]/(1+0.5*dt*fabs(u[i])/h[i]); }  
~~~

v) You can put granular  friction...

ex [http://basilisk.fr/sandbox/M1EMN/Exemples/savagestaron.c]()

etc.


v) You may change the boundary conditions,
 
 here BC are output BC, cells `u[0]`and `u[nx+1]` are ghost cells
 
no flux, we copy the values
 
 
~~~bash
 hn[0]=hn[1];
 un[0]=un[1];
 hn[nx+1]=hn[nx];
 un[nx+1]=un[nx];
~~~
 
 
It is changed to have no penetration at the walls: a zero velocity at the final face
 
 dirichlet 0: the final face value is the linear extrapolation of the `nx` and the `nx+1`
 $$u= \frac{u_{nx}+u_{nx+1}}{2}$$
 
~~~bash
 hn[0]=hn[1];
 un[0]=-un[1];
 hn[nx+1]=hn[nx];
 un[nx+1]=-un[nx];
~~~
 
 
 
# Links

the same example of dam break with Basilisk

 * [dam break](http://basilisk.fr/sandbox/M1EMN/Exemples/damb.c)
 
 * the same in [python...](https://colab.research.google.com/drive/1960Q9Cgu9anAv6wfB9MblsBFo0BvWdrO)

  
 
# Biblio
 
* A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows
   Emmanuel Audusse, Franccois Bouchut, Marie-Odile Bristeau, Rupert Klein and Benoıt Perthame
*  Olivier Delestre [Thèse](https://hal.archives-ouvertes.fr/tel-00531377/document)
* [cours PYL](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv)
* [cours PYL](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf) sur numérique
* [Hubert Chanson](https://books.google.fr/books?id=VCNmKQI6GiEC&printsec=frontcover&dq=river+flows+shallow+water+chanson&source=gbs_similarbooks_s&hl=fr#v=onepage&q&f=false)
The hydraulics of open channel flow: an introduction ; basic principles"
Elsevier Butterworth-Heinemann, Amsterdam, The Netherlands. 2004

// version init mai 2012 version mars 2013
// PYL
*/

