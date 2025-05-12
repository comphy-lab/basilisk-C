/**
#    Turbulent fluid "dam break problem", 

 
The classical Shallow Water Equations  in 1D with turbulent friction. 
$$
 \left\{\begin{array}{l}
         \partial_t h+\partial_x Q=0\\
	\partial_t Q+ \partial_x \left[\dfrac{Q^2}{h}+g\dfrac{h^2}{2}\right]
= -C_f u^2
        \end{array}\right. 
$$
with at  $t=0$, a lake on the left $h(|x|<0,t=0)=1$ and zero water $h(|x|>0,t=0)=0$ for $x>0$.
As the flow is a bit viscous, the convective term is important.
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"

double tmax;

u.n[left]   = neumann(0);
u.n[right]  = neumann(0);  
/** 
start by a lake on the left and dry on the right, the soil is flat $z_b=0$
*/
event init (i = 0)
{
  foreach(){
    zb[] =  0;
    h[] =  (x<0);
    u.x[] =0;
   }

#ifdef gnuX
  printf("\nset grid\n");
#endif    
}

/**
 
position of domain `X0`, length `L0`, no dimension `G=1`,  
run with 512 points (less is not enough)
*/

int main() {
  X0 = -5.;
  L0 = 10.;
  G = 1;
  N = 512;
  tmax=2;
  DT = .1;
  run();
}
/**
  We use a simple implicit scheme to implement quadratic bottom
  friction by splitting i.e.
  $$
  \frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
  $$
  with $C_f=.05$. */
event friction(i++){
   foreach() {
    double a = h[] < dry ? HUGE : 1. + .05*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
  }
}    
/**
Output in gnuplot if the flag `gnuX` is defined
*/
#ifdef gnuX
event plot (t<tmax;t+=0.01) {
    printf("set title 'rupture de barrage en 1D --- t= %.2lf '\n"
      "h(x)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t));t= %.2lf  ; "
      "p[%g:%g][-.5:2]  '-' u 1:($2+$4) t'free surface' w l lt 3,"
      "'' u 1:($2*$3) t'Q' w l lt 4,\\\n"
      "'' u 1:4 t'topo' w l lt -1,\\\n"
      "'+' u (0):(0.296296296296296) t'theo Qmax ' , h(x) t 'theo h(x)'\n",
           t,t,X0,X0+L0);
    foreach()
    printf ("%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
    printf ("e\n\n");
}
/**
else put the outputs in `out`
*/
#else
event plot(t<tmax;t+=.5 ) {
    foreach()
    printf ("%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
}
#endif
/**
 
## Run
To compile and run with gnuplot:

~~~bash
  qcc -DgnuX=1   -O2   -o damb damb_dressler.c -lm;  ./damb_dressler | gnuplot
~~~

To generate the plot with   `out`

~~~bash
 qcc  -O2   -o damb_dressler damb_dressler.c -lm;  ./damb_dressler > out  
~~~

## Plots
 
Plot of the surface and flux at time $t=1,2$:

comparison of heights with and without friction

~~~gnuplot  result for h, comparison between computation in red and ideal fluid theory blue and green   
set xlabel "x"
h(x,t)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t))
u(x,t) =2./3*(1+x/t)+0.05*F1t(x,t)
F1t(x,t)=t*(-108./7/(2-x/t)**2. +12/(2-x/t)-8./3+8.*sqrt(3)/189*(2-x/t)**1.5)


F2t(x,t)= ((6./5)/(2-x/t)-2./3)*t+4*sqrt(3)/189.*((2-x/t)**(3./2))*t
 
#h(x,t)=((2./3*(1-(x-0)/(2*t))) +0.05*F2t(x,t) )**2.
p [-2:3][0:1] 'out' u 1:(($5==.5)||($5==1)||($5==1.5)? $2 :NaN)w p, h(x,.5),h(x,1.),h(x,1.5) 
 
~~~

using the Dressler 1952 solution for velocities $u$ ans $h$:
$$
u(x,t)= 2./3*(1+x/t)+C_fF_1(x,t)t
$$
$$c= 1./3*(2-x/t) + C_fF_2(x,t)t$$
$$h = c^2$$
with
$$F_1(x,t)=(-\frac{(108./7)}{(2-x/t)^2} +\frac{12}{(2-x/t)}-8./3+8.\sqrt(3)/189*(2-x/t)^{3/2})$$
  $$F_2(x,t)= ((6./5)/(2-x/t)-2./3)t+4(\sqrt(3)/189)((2-x/t)^{(3./2)})$$
  See Dressler, see as well Chanson and  Delestre (who put a 1./135 in $F_2$).
  
~~~gnuplot  result for velocity, comparison between computation in red and Dressler theory blue and green   
set xlabel "x"
h(x,t)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t))
u(x,t) =2./3*(1+x/t)+0.05*F1t(x,t)
F1t(x,t)=t*(-108./7/(2-x/t)**2. +12/(2-x/t)-8./3+8.*sqrt(3)/189*(2-x/t)**1.5)
p[-2:4][0:2] 'out' u 1:(($5==.5)||($5==1)||($5==1.5)? $3 :NaN)w p, u(x,.5),u(x,1.),u(x,1.5) 
 
~~~

 ~~~gnuplot  result for h, comparison between computation in red and Dressler theory blue and green   
set xlabel "x"
h(x,t)=(((x-0)<-t)+((x-0)>-t)*(2./3*(1-(x-0)/(2*t)))**2)*(((x-0)<2*t))
c(x,t) =  sqrt(h(x,t)) +0.05*F2t(x,t)
F1t(x,t)=t*(-108./7/(2-x/t)**2. +12/(2-x/t)-8./3+8.*sqrt(3)/189*(2-x/t)**1.5)
F2t(x,t)=t*(6./5/(2-x/t) -2./3+4.*sqrt(3)/89*(2-x/t)**1.5)
hd(x,t)=c(x,t)*c(x,t)*((x-0)<2*t)
p[-2:4][0:2] 'out' u 1:(($5==.5)||($5==1)||($5==1.5)? $2 :NaN)w p, hd(x,.5),hd(x,1.),hd(x,1.5) 
 
~~~       
    
 
    
#Links
 
see the same
  [with friction](http://basilisk.fr/sandbox/M1EMN/Exemples/damb_dressler.c)
 and  see non viscous dam break with [standard C](http://basilisk.fr/sandbox/M1EMN/Exemples/svb.c) 
and with [Basilisk](http://basilisk.fr/sandbox/M1EMN/Exemples/damb.c)
 
 
#Bibliographie
 
* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC
* [Dressler 1955](http://nvlpubs.nist.gov/nistpubs/jres/049/jresv49n3p217_A1b.pdf)
* [Chanson](https://books.google.fr/books?id=VCNmKQI6GiEC&pg=PA357&dq=dressler+solution+chanson&hl=en&sa=X&ved=0ahUKEwjj4Ne26_zYAhWBRhQKHW_9ArYQ6AEIJzAA#v=onepage&q=dressler%20solution%20chanson&f=false) p 358
*  CHANSON, H. (2006). "Solutions Analytiques de l’Onde de Rupture de Barrage sur Plan Horizontal et Incliné." ('Analytical Solutions of the Dam Break Wave Problem on Horizontal and Inclined Inverts.') Jl La Houille Blanche, No. 3, pp. 76-86 (ISSN 0018-6368)
            [Chanson](https://espace.library.uq.edu.au/data/UQ_7977/UQ7977_OA_preprint.pdf?Expires=1744301419&Key-Pair-Id=APKAJKNBJ4MJBJNC6NLQ&Signature=dqKffRJz74PUeno2wWa0~eRPyVZd0eNfNBsGfcmHm9vXKLdiMkPeuIqqslB4rCVXaIBcY1cZmLooP~6nwxYGV1O2qswXBmLn58bFjauK7TpNuVLlL19RTCWsB-b853YRRp2hPagJJFBslJYKulCyIVRxiTRTOEeuvo~VGoHzv7eBZ2uQh5xNAkANhdwfzIMQG7tMgLn8SL3ey37ZY6B7Wd1-7dl-YMAt6QS~Oe0CBlth5tkMrPNBeLPdulYrJRzcY~~4xmHEpfLonxz2NSvJqtXfnAHdi1pX5QpUhzGJOgVvQf8sZo1vWuyJX3GzX~8NzwNavP~x2ONWWa7Naj6ZLg__)
* [SWASHES](https://hal.archives-ouvertes.fr/file/index/docid/628246/filename/SWASHES.pdf)   page 24
* [thèse Delestre p 246](https://hal.archives-ouvertes.fr/tel-00531377/document)

Version 1: Paris 2018 

*/

   

 
 
 