/**
# Resolution of  the problem of unsaturated-saturated water table recharge problem
 
 
 What is the unsteady flow in an unsaturated porous media, how does the water table emerge from the infiltration of a river?
 
This problem is a mix of Vauclin (saturated/ unsaturated) and Todd Bear (River topology).

 The river full of water imbibes the soil which is at first unsaturated.
 As time evolves, saturation (i.e. $p>0$) arises.
 A water table is generated at the bottom over the impervious soil in $y=0$. The thickness of the porous media is $H_d$.
 
We model this with Darcy unsaturated theory.
 We have to solve a Darcy problem without dimension:
 $$u= - \beta(p) \frac{\partial p }{\partial x }\text{ and }
   v =- \beta(p)  (\frac{\partial p }{\partial y } - 1 )$$
$\beta$ is proportional to the permeability or the conductivity  and is function of the pressure $p$.

Conservation of water is
 $C(p) \frac{\partial h}{\partial t }+\frac{\partial u}{\partial x }+\frac{\partial v}{\partial y }=0$
 so that we have to solve an unsteady non linear diffusion problem:
$$C(p) \frac{\partial h}{\partial t } =
\frac{\partial }{\partial x } (\beta(p)\frac{\partial p } {\partial x })
+
\frac{\partial }{\partial y } (\beta(p)  \frac{\partial p } {\partial y } - 1)$$
with $C(p)$ and $\beta(p)$ given functions, as $p$ is without dimension $p$ is as well $h$.
  The so called 'hydraulic head' is then $H=p+y$ (again $p$ is without dimension).
 Note that  if pressure is negative, the porous media is unsaturated, if pressure is positive the porous media is saturated.
 The $ad$ $hoc$ functions $C(p)$ and $\beta(p)$ are empirical relations:
 $$ C = 1 \text{ for saturated, and for unsaturated: } C= \frac{1}{1+ \alpha_1 (|p|^{n_1})}$$
 $$ \beta  = 1 \text{ for saturated, and for unsaturated: } \beta=\frac{1}{1+ \alpha_2 (|p|^{n_2})}$$
 empirical values $\alpha_1$,  $\alpha_2$, $n_1$ and $n_2$ are provided...
 and the density is 1 by adimensionalisation.
 
~~~gnuplot
 reset
 set style fill transparent solid 0.5 noborder
 set style function filledcurves y1=0
 unset colorbox
 set lmargin 6
 set xlabel "x"
 set label "porous media unsaturated" at 5,4
 set label "water comming from the river" at 8,5
 set label "iso line p(x,y,t)=0, water table" at 5,1.25
 set object 1 rect from 11,5.5 to 12,7  fc rgb "blue"
 set arrow from 11,5.5 to 10,4  lc rgb 'black'
 set arrow from 11,6 to 9.6,6  lc rgb 'black'
 set arrow from 11.8,5.5 to 11.8,4.5  lc rgb 'black'
 plot [0:12][0:7] 1+exp(-.3*(x-12)*(x-12)) fs solid 1.0 lc rgb "blue" not, (x>11)fs solid 1.0 lc rgb "blue" not
~~~
 
*/
#include "run.h"
#include "poisson.h"
#define level 6
 scalar p[],C[],source[],gammad[],H[],wet[];
 double x0=11;
 double Hd=7.;
 double Hf=5.5;
 double tmax,dt;
 face vector beta[];
 mgstats mgp;
 face vector u[];
 /**
 
 ## boundary conditions
 
 Domain of size $L_0$,
 at the bottom  and top hydrostatic pressure gradient, and zero pressure at left and at the right,
In the upper right corner there is a river which imposes a hydrostatic pressure at the boundary.
 */
u.n[right]  =  neumann(0);
u.t[right]  =  neumann(0);
u.n[left]   =  neumann(0);
u.t[left]   =  neumann(0);
u.n[bottom] =  dirichlet(0);
u.t[bottom] =  neumann(0);
u.n[top]    =  neumann(0);
u.t[top]    =  neumann(0);

p[left]  = neumann(0);
p[right] = neumann(0);
p[top]   = neumann(-1);
p[bottom]= neumann(1);

/**
 domain is `L0 x Hd`, the river pressure
 */
bid river;
u.t[river] = neumann(0);
p[river] =dirichlet(Hd-y) ;

int main()
{
    L0=12.;
    Y0=0;
    X0=0;
    init_grid (1 << level);
    tmax = 15;
    run();
}
/**
 gravity and initial values, we start by an unsaturated soil, the pressure is negative.
 */
const face vector g[] = {0.,-1.};

event init (i = 0) {
    mask (y> Hd ? top: none);
    mask ( x> x0 && y>Hf ? river : none);
    
    foreach_face() {
        beta.x[] = 10./11;
    }
    boundary  ((scalar *){beta});
    
    foreach(){
        p[] = -1 ;
        source[] = 0.;
        C[] = 10./11;
    }
    boundary ({p,source,C});
    
    // We may change the default gradient function (used for advection) to minmod-limited (rather than the centered default).
    //  gradient = minmod2;
}

event pressurevelocities (i++; t<tmax)
{
    /**
     evaluate physical quantities $\beta$ permeability and $\rho$ density as function of
     pressure, if pressure is negative, the porous media is unsaturated, if pressure is positive the porous media is saturated
     $$ C = 1 \text{ for saturated, and for unsaturated: } C= \frac{1}{1+ \alpha_1 (|p|^{n_1})}$$
     $$ \beta  = 1 \text{ for saturated, and for unsaturated: } \beta=\frac{1}{1+ \alpha_2 (|p|^{n_2})}$$
     empirical values $\alpha_1$,  $\alpha_2$, $n_1$ and $n_2$ are provided...
     and the density is 1 by adimensionalisation
     */
    foreach(){
        wet[] =(p[]>0?0:1);
        C[] =  10/(10+ pow(fabs(p[]<0?p[]:0),2.90));
    }
    boundary  ({C});
    foreach_face() {
        beta.x[] = 10/(10+ pow(fabs(p[]<0?p[]:0),5));
    }
    boundary  ((scalar *){beta});
    /**
     
     We have to solve
     $$ 0 = -  \overrightarrow{\nabla}   p -   \frac{1}{\beta} \overrightarrow{  u} - \overrightarrow{ e_y}
     \text{  with }  C\frac{\partial p}{\partial t} + \overrightarrow{  \nabla} \cdot \overrightarrow{  u} =0$$
     this gives :
     $$ C\frac{\partial p}{\partial t} =  \nabla \cdot (\beta \nabla p  )   $$
     implicit discretisation gives a Poisson problem
     $$\nabla \cdot (\beta \nabla p^{n+1}  ) - C \frac{p^{n+1}}{\Delta t} = - C \frac{p^{n}}{\Delta t}$$
     time step is $\Delta t =$ 0.025*/
    dt = dtnext (0.025);
    foreach(){
        source[]=-C[]*p[]/dt;
        gammad[] = -C[]/dt;
    }
    boundary  ({source,gammad});
    /**
     We solve Poisson equation
     $\nabla \cdot (\beta \nabla p  ) + \gamma_d \; p = s$
    
     with [http://basilisk.fr/src/poisson.h](http://basilisk.fr/src/poisson.h)
     */

    mgp = poisson (p, source, beta , gammad);

    /** the velocity is then computed from the gradient
     $$ \overrightarrow{  u} = - \beta  \overrightarrow{\nabla} - \beta \overrightarrow{ e_y}$$
     
     (note the `face_gradient_x p=  ((p[i] - p[i-1])/Delta)`) see
     [http://basilisk.fr/src/poisson.h](http://basilisk.fr/src/poisson.h)
     
     */
    foreach_face()
    u.x[] = (- beta.x[]*face_gradient_x (p, 0) + beta.x[]*g.x[]);
    boundary ((scalar *){u});
 /** Hydro head */
    foreach(){
        H[] = p[] + y ;   }
    boundary({H});
}
/**
 error
 */
event logfile (i++)
{
    stats s = statsf (p);
    fprintf (stderr, "%d %g %d %g %g %g\n",
             i, t, mgp.i, s.sum, s.min, s.max);
    fprintf (stderr,"%g \n",dt);
}
/**
 Save in a files... films...
 */
#if 0
event gfsview (i += 10) {
    static FILE * fp = popen ("gfsview2D gfsview.gfv", "w");
    output_gfs (fp, t = t);
}
#endif

event sauve (t+=.1 )
{
    FILE *  fpc = fopen("pressure.txt", "w");
    output_field ({p,u.x,u.y}, fpc, linear = true);
    fclose(fpc);
    
    FILE *  fpq = fopen("qh.txt", "w");
    double dx=.1;
    for (double x = 0 ; x < L0; x += dx){
        for (double y = 0 ; y < Hd; y += dx){
            fprintf (fpq, "%g %g %g %g %g \n",
                     y,interpolate (p, L0-.01, y), interpolate (beta.x, L0-.01, y),
                       interpolate (p, .01, y), interpolate (beta.x,.01, y)
                     );}
        fprintf (fpq,"\n");}
    fclose(fpq);
}
// 'avconv' ne marche plus???
event movie (t+=.02;t<=tmax)
{
    //static FILE * fp1 = popen ("ppm2mpeg > mouille.mpg", "w");
    output_ppm (wet,file="mouille.mp4",min = 0.,max =2,mask = wet, n = 1024,linear = true,
                box = {{0.05,-1},{L0-0.05,Hd}});
    
  //  static FILE * fp2 = popen ("ppm2mpeg > p.mpg", "w");
    output_ppm (p,file="p.mp4",min = 0.,max = 2,mask = wet, n = 1024,linear = true,
                box = {{0.05,-1},{L0-0.05,Hd}});
}
event photo ( i++ ){
    if (i == (int)(1/dt))
      {output_ppm (wet, file = "wet1.png",  min = 0, max = 2, n = 256, linear = true);
       output_ppm (p,  file = "pres1.png",  min = 0, max = 2, n = 256, linear = true);}
    if (i == (int)(2/dt))
      {output_ppm (wet, file = "wet2.png",  min = 0, max = 2, n = 256, linear = true);
       output_ppm (p, file = "pres2.png",  min = 0, max = 2, n = 256, linear = true);}
    if (i == (int)(3/dt))
      {output_ppm (wet, file = "wet3.png",  min = 0, max = 2, n = 256, linear = true);
       output_ppm (p, file = "pres3.png",  min = 0, max = 2, n = 256, linear = true);}
    if (i == (int)(4/dt))
      {output_ppm (wet, file = "wet4.png",  min = 0, max = 2, n = 256, linear = true);
       output_ppm (p, file = "pres4.png",  min = 0, max = 2, n = 256, linear = true);}
    if (i == (int)(5/dt))
       {output_ppm (wet, file = "wet5.png",  min = 0, max = 2, n = 256, linear = true);
        output_ppm (p, file = "pres5.png",  min = 0, max = 2, n = 256, linear = true);}
    if (i == (int)(10/dt))
    {output_ppm (wet, file = "wet10.png",  min = 0, max = 2, n = 256, linear = true);
        output_ppm (p, file = "pres10.png",  min = 0, max = 2, n = 256, linear = true);}
    if (i == (int)(15/dt))
        {output_ppm (wet, file = "wet15.png",  min = 0, max = 2, n = 256, linear = true);
        output_ppm (p, file = "pres15.png",  min = 0, max = 2, n = 256, linear = true);}
}

/**
# Run
 To compile and run:
 
~~~bash
 make vauclin.tst; make vauclin/plots ; make vauclin.c.html;
 source ../Exemples/c2html.sh vauclin
~~~
 
# Results

## Movies/ images
 
Two movies of the infiltration process,
 
 

[![](./vauclin/pres10.png)](./vauclin/p.mp4)


Animation  of the infiltration: pressure,   pressure evolution and  stabilisation.
 (click on image for animation)
 
 
A film [(here)](vauclin/mouille.mp4)  of the infiltration: unsaturated/ saturated development, blue is wet, green is unsaturated (click for animation)

## snapshots
 

### pressure as function time
 
 fot t=1, 2 3 5  10 15
 
 
  ![t=1](./vauclin/pres1.png)
  ![t=2](./vauclin/pres2.png)
  ![t=3](./vauclin/pres3.png)
  ![t=5](./vauclin/pres5.png)
  ![t=5](./vauclin/pres10.png)
  ![t=15](./vauclin/pres15.png)

### progression of the front
 
  fot t=1, 2 3 4 5 10  and 15
 
 
 ![un/saturated](./vauclin/wet1.png)
 ![un/saturated](./vauclin/wet2.png)
 ![un/saturated](./vauclin/wet3.png)
 ![un/saturated](./vauclin/wet4.png)
 ![un/saturated](./vauclin/wet5.png)
 ![un/saturated](./vauclin/wet10.png)
 ![un/saturated](./vauclin/wet15.png)

## Plots
 
 final pressure
 
~~~gnuplot
reset
 set pm3d map
 set palette rgbformulae 22,13,-31;
 unset colorbox
 set xlabel "x  iso p"
 set size (.5*1.6),1
 splot [][ :][-10:15] 'pressure.txt' u 1:2:3   not
 reset
~~~
 
 Figure of iso pressure at final time
 
 
 ~~~gnuplot pressure
 reset
 set view map
 set size (.5*1.6),1
 set label "river" at 11,8.5
 set label "soil surface" at 2,7.5
 set label "porous media" at 2,5
 set label "impervious rock" at 2,0.5
 unset key
 unset surface
 set contour base
 set cntrparam levels incremental -15,1,15
 splot [][:] 'pressure.txt' u 1:2:($3)  w l not
 ~~~
 
 ~~~gnuplot
 reset
 set isosample 100, 100
 set table 'test.dat'
 splot[0:12][0:10][-10:5] 'pressure.txt' u 1:2:(abs($3)<100?$3:NaN)
 unset table
 
 set contour base
 set cntrparam level incremental -9, .25, 9
 unset surface
 
 set table 'cont.dat'
 splot[0:12][0:10] [-10:5]'pressure.txt'  u 1:2:3
 unset table
 
 reset
 unset key
 unset colorbox
 set palette rgbformulae 33,13,10
 p 'test.dat' with image, 'cont.dat' w l lt -1 lw .15
 ~~~
 
 plot of Hydraulic Head  $H=p/(\rho  g)+z$:
 
 ~~~gnuplot Hydraulic Head
 L0=1.
 reset
 set view map
 set size (.5*1.6),1
 unset key
 unset surface
 set contour base
 set cntrparam levels incremental 0,.5,15
 splot [][:] 'pressure.txt' u 1:2:($3+$2)  w l not
 reset
 ~~~
 
Cuts of pressure at both ends of the domain, showing that at $x=0$ saturation is partial (last time step)
 
 ~~~gnuplot
 set ylabel "y"
 set xlabel "p"
 reset
 p[-4:4][:8]'qh.txt' u 2:1 t'p(x=L0,y)' ,''u 4:1 t 'p(x=0,y)',0 not
 reset
~~~
 
## Links
 
 * [https://stackoverflow.com/questions/20977368/filled-contour-plot-with-constant-color-between-contour-lines]()
 
 * see [http://basilisk.fr/src/hele-shaw.h](http://basilisk.fr/src/hele-shaw.h)
 
 * see [http://basilisk.fr/sandbox/M1EMN/Exemples/darcyLambSneddon.c]()
 
 * see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/dupuit2D.c]()
 
 * see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/richards.c]() Richards problem
 
 * see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/toddbear59.c]()

 * see [http://basilisk.fr/sandbox/M1EMN/HYDGEO/vauclin.c]()
 
## Bibliography

 * M. Vauclin , D. Khanji , G. Vachaud
 ["Experimental and numerical study of a transient, twodimensional unsaturatedsaturated water table recharge problem"](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/WR015i005p01089)
 `Vauclin_et_al-1979-Water_Resources_Research.pdf`
 
 * [https://www.researchgate.net/publication/320515457_Numerical_Solution_of_Richards%27_Equation_A_Review_of_Advances_and_Challenges]()
 
 * Joseph Bear "Dynamic of fluids in porous media" p 238
 
 * [Todd Bear 1959](https://escholarship.org/uc/item/1nx0q3dd) River Seepage investigation
 
 * [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv_aquifere.pdf)
 "Ecoulements en milieux naturels: écoulements en milieux souterrains" Cours MSF12, M1 UPMC
 
 OK
 
 

 
*/
