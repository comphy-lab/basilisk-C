/**
 
# Gnuplotting


 This is a small   tutorial for `gnuplot` use with `Basilisk`. The numerical results are generated with heat equation, so
 this is exactly [./heat.c]()
 but with more and more graphs and graphics.
 We show the difference between files and terminal outputs.
 
 usefull for [http://m2.basilisk.dalembert.upmc.fr/sandbox/PYL/README]()
 
## Problem
 
 Any other problem is convinient, but we have to generate data! 
 So we consider 
 Heat with finite volumes, initial temperature a "Dirac", in fact a thin rectangle so that $\int_{-\infty}^{\infty}Tdx=1$.
 explicit step, discrete flux across interface between ` ì-1` and ` ì`:
 $$q_i = -  \frac{ T_{i} -  T_{i-1}}{\Delta }$$

 and then, balance between flux leaving $q_{i+1}$ and flux entering $q_i$ the cell
 $$T_i^{n+1}= T_i^{n }  - (\Delta t ) ( q_{i+1} - q_i)/\Delta$$
 
*/
#include "grid/cartesian1D.h"
#include "run.h"
scalar T[],q[];
double dt;
int main() {
    L0 = 10.;
    X0 = -L0/2;
    N = 100;
    DT = .201*(L0/N)*(L0/N)/2 ;
#define EPS 0.5
    run();
}
/** BC */
T[left]=neumann(0);
T[right]=neumann(0);
q[left]=dirichlet(0);
q[right]=dirichlet(0);
/** Initial "door" of $T$ */
event init (t = 0) {
    foreach()
    T[] =  1./EPS*(fabs(x)<EPS)/2;
    boundary ({T});
}
/** Finite volumes ($T$ at the center, $q$ at the faces (see section `Drawing a ))
$$T^{n+1}_i=T^n_i - \Delta t \frac{q^n_{i+1}-q^n_{i}}{\Delta x}$$ */
event integration (i++) {
    double dt = DT;
    scalar dT[];
    dt = dtnext (dt);
    foreach()
    q[]=-(T[0,0] - T[-1,0])/Delta; 
    boundary ({q});
    foreach()
    dT[] = - ( q[1,0]  - q[0,0] )/Delta;
    foreach()
    T[] += dt*dT[];
    boundary ({T});
}
/**
 That is just to generate physical results.

## Saving the values
 
At the end we process data. The most simple way is to "print" them in the terminal.
We
 print data through standard output `stdout` in `C`.
 The instruction
 `printf ()`  is equivalent to `fprintf (stdout,)`
 here we print all the $x$ and value $T$ with a new line `\n` after each new couple $(x_i, T_i)$.
 */
event printdatastdout (t = 0.5) {
    foreach()
      printf ( "%g %g \n", x, T[]);
}
/**
 
The values are printed on the terminal. By convention, when you compile withe `make` or when
 you run in the web brouser `Basilisk` saves the standard output of the  terminal into a file named `out`.
 
  In fact Basilisk does something like `./gnuplot_examples 2>log 1>out`

 In this file `out` we have
 the $x_i$ in the first column and the field in the second $T_i$, with all the température vzalues
 
 starting from the first $X_0+\Delta/2$ up to the  $N$th value $X_0 + L -\Delta/2$,
 where $\Delta = L/N$
 <hr>
 <table>
 <tr>
 <td>$x_1=X_0+\Delta/2$</td>
 <td>$T_1$</td>
 </tr>
 <tr>
 <td>$x_2=x_1+\Delta$</td>
 <td>$T_2$</td>
 </tr>
 <tr>
 <td>$x_1$</td>
 <td>$T_1$</td>
 </tr>
 <tr>
 <td>$x_{N}$</td>
 <td>$T_{N}$</td>
 </tr>
 </table>
 

The most simple command in `gnuplot` is just `plot`, or simply `p`
~~~gnuplot
 p 'out'
~~~

 the same with line between points and no point marks,  the ranges are imposed (-3 3)(0 1),
 note the labels of axis and the label of the curve

~~~gnuplot
 set xlabel 'x'
 set ylabel 'T'
 p[-3:3][0:1]'out' t'curve'w l
~~~
Another possibility is to create a file and save the data in the file
*/
event printdata (t = 1){
    char s[80];
    sprintf (s, "shape.txt");
    FILE * fp = fopen (s, "w");
    foreach()
    fprintf (fp, "%g %g %g \n", x, q[], T[]);
    fclose(fp);
}
/**
We can plot the two files on the same plot.
 The second file has `q` in second colomn and `T` in the 3rd column. We have to use columns 1 and 3: hence  we ` use ` the command ` u 1:3 ` to plot `T` and compare to the ` out` file.
 
 we plot a third curve, the second colomn `q`. We just put two quotes to say that it is the same file. Then no need to put ` u 1:2 ` as it is implicit. Note that the first curve plotted  is red, the second is green, and the third is blue.  Colors vary corresponding to the different terminal, that is a drawback of gnuplot.
 
~~~gnuplot
 p 'out' u 1:2 w l,'shape.txt' u 1:3,'' w l
~~~


Here the field is saved at several time steps.
 Each different time step is separated by a blank line : carriage return `\\n'.
 We have a new "block" for every time step.
 */
event printdatas (t +=0.05;t<2){
    char s[80];
    sprintf (s, "shapes.txt");
    static FILE * fp = fopen (s, "w");
    foreach()
      fprintf (fp, "%g %g %g \n", x, T[], t);
    fprintf (fp, "\n");
    fflush(fp);
}
/**
 
 `plot` superposes all the curves
 
~~~gnuplot
 plot 'shapes.txt'  w l
~~~

 We can  select some curves with the ` every` command. Examples:
 
~~~bash
 
 every ::3 skip the first 3 lines
 every ::3::5  plot from the 4-th to 6-th lines
 every ::0::0  plot the first line only
 every 2::::6  plot the 1,3,5,7-th lines
 every :2  plot every 2 data block

 
 p'shapes.txt' ev :5::0::35 w l // plot from block 0 to the 35th every 5 block
~~~
 
~~~gnuplot
 p'shapes.txt' ev :5::0::35 w l
~~~
 
~~~bash
every :::5::8 plot from 5-th to 8-th data blocks
~~~
 
~~~gnuplot
p'shapes.txt' every :::5::8  w l
~~~
 
 

 
# computations with the plotted values
 
 
 The self similar analytical solution is
 $$\theta = \frac{1}{2 \sqrt{\pi t} } e^{-x^2/4 t}$$
 the solution can be written in self similar variables: we plot $\theta \sqrt{\pi t} $ as a function of $x/\sqrt(4t)$ and superpose
 $e^{-x^2}$
 
 Note that we use points of type 5, and size 0.5, le width of the line `lw` is 2
 
~~~gnuplot
 reset
 p[-5:5][:1]'shapes.txt' u ($1/sqrt(4*$3)):($2*sqrt($3*pi)*2) t 'numerical' w p pt 5 ps 0.5, exp(-x*x) t 'self similar' lw 2
~~~

 

 There is another output: the error output `stderr`, here we
 print data in the "error" output. Note that when in a terminal they printed with no difference.
 In fact Basilisk does something like `./gnuplot_examples 2>log 1>out`

Note that here we put two newlines '\\n\\n'   to separate the blocks */

event printdatastderr (t += 0.1; t < 1) {
    foreach()
    fprintf (stderr, "%g %g %g\n", x, T[], t);
    fprintf (stderr, "\n\n");
}
/**

 When two newlines '\\n\\n'separate the blocks we can plot each block using   `index`
 
~~~gnuplot
 p 'log'  index 10 u 1:2,  'log' w l
~~~

 Here we even can do a loop
 
~~~gnuplot
 p for [i=2:4] 'log' index i using 1:2 with lines
~~~

 note the colors: red, green, blue...

 
 
# Drawing a graphics

 
 Note that we change the aspect ratio, the number of samples for the function `f`.
 We put text at several places, and draw arrows.
 
~~~gnuplot
 reset
 set size ratio .5
 set samples 9
 f(x)=1.25*exp(-x*x/2.3)
 set label "T i-1" at 1.5,3.1
 set label "T i" at 2.5,2.5
 set label "T i+1" at 3.5,2.25
 set xtics ("i-2" 0.5, "i-1" 1.5, "i" 2.5,"i+1" 3.5,"i+2" 4.5,"i+3" 5.5)
 set arrow from 2,1 to 2.5,1
 set arrow from 3,1 to 3.5,1
 set label "q i" at 2.1,1.25
 set label "q i+1" at 3.1,1.25
 
 set label "x i-1/2" at 1.5,0.25
 set label "x i" at 2.4,0.25
 set label "x i+1/2" at 3.,0.25
 
 set label "x"  at 0.5,1.5+f(0)
 set label "x"  at 1.5,1.5+f(1)
 set label "x"  at 2.5,1.5+f(2)
 set label "x"  at 3.5,1.5+f(3)
 set label "x"  at 4.5,1.5+f(4)
 set label "x"  at 5.5,1.5+f(5)
 p[-1:7][0:4] 1.5+f(x) w steps not,1.5+f(x) w impulse not linec 1,'out' u ($1+.5):($2*3.4+1.5) t'T' w l
 reset
~~~

this graphics corresponds to 
$$T^{n+1}_i=T^n_i - \Delta t \frac{q^n_{i+1}-q^n_{i}}{\Delta x}$$ 

# multi plot
 
 putting several plots on the same figure
 
 ~~~gnuplot
 reset
 set multiplot layout 2,2
 set xlabel 'x'
 plot [][:1]'shapes.txt' every :::0::0 w l lc 2 not
 plot [][:1]'shapes.txt' every :::5::5 w l lc 2 not
 plot [][:1]'shapes.txt' every :::10::10 w l lc 2 not
 plot [][:1]'shapes.txt' every :::20::20 w l lc 2 not
 unset multiplot
 
 ~~~
 
# 3 D plot
 
 plot in 3D, using `splot` or `sp` :
 
~~~gnuplot
 splot 'shapes.txt' u 3:1:2 w l
~~~

 can change the `view` and hide the hidden part:
 
~~~gnuplot
 set hidden3d
 set view 45,45
 splot  'shapes.txt' u 3:1:2 w l
~~~

 draw color maps, I do not like the standard one, I prefer this one in RGB
 
 
~~~gnuplot
 reset
 set pm3d; set palette rgbformulae 22,13,-31;unset surface;
 set ticslevel 0;
 unset border;
 unset xtics;
 unset ytics;
 unset ztics;
 unset colorbox;
 set view 0,0
 sp[0:.5][-3:3][] 'shapes.txt' u 3:1:2 not
~~~
 
If you prefer gray scale
 
~~~gnuplot
 set pm3d map
 set palette gray negative
 unset colorbox
 set tmargin at screen 0.95
 set bmargin at screen 0.15
 set rmargin at screen 0.95
 set lmargin at screen 0.15
 set xlabel "t"
 set ylabel "x"
 unset key
 splot [:1][-2:2][]'shapes.txt' u 3:1:2
~~~

 
 
# Animations

 Gif animated! 
 
 We have several possibilities, at leat 
 two possibilities, one with `every`....
 
~~~gnuplot
 reset
 set xlabel "x"
 !echo "p[:][:1]'shapes.txt' ev :::k::k title sprintf(\"t=%.2f\",k*.05) w l;k=k+1 ;if(k<40) reread " > mov.gnu
 set term gif animate;
 set output 'slump.gif'
 k=0
 l'mov.gnu'
~~~

 Gif animated of the slump (reload to refresh) or click on image for animation:
 
 

or another possibility for animation
using   `index` instead

~~~bash
reset
!echo "print k "> movi.gnu
!echo "p[:][0:1]'log' i k u 1:2t'T(x,t)'w l \nk=k+1 \nif(k<40) reread" >> movi.gnu
set term gif animate;
set output 'slumpi.gif'
k=0
l'movi.gnu'
~~~

etc.

 
# next level:
 
 It is nice to have a plot on the fly of the output using gnuplot.
  To do that, we put 
  a flag `gnuX`
 
 see [../Exemples/damb.c]() and many others
 
 Output in gnuplot if the flag `gnuX` is defined,
 
~~~bash
 qcc  -DgnuX=1  gnuplot_examples.c -lm
 ./a.out | gnuplot -persist
~~~
 
 this reads: 
*/
#ifdef gnuX
event outputgnu (t += .1 ) {
    fprintf (stdout, "p[-4:4][0:2]  '-' u 1:2   not  w   l \n");
    foreach()
     fprintf (stdout,"%g %g \n", x, T[]);
    fprintf (stdout, "e\n");
}
#endif
/**

there are many other possibilities with gnuplot....


 
# Run
 
 To run the code and generate all the graphics
 
 
~~~bash
 make gnuplot_examples.tst;make gnuplot_examples/plots
 make gnuplot_examples.c.html ; open gnuplot_examples.c.html
~~~
 
 
# Links

 * [http://basilisk.fr/sandbox/amansur/Gnuplot]()
 
 * [heat.c]() the explicit heat equation
 
 * [heat_imp.c]() the implicit heat equation
 
 * [http://basilisk.fr/sandbox/M1EMN/Exemples/herschel-column-noSV.c]()
 
 * [../Exemples/damb.c]()
 
 * etc
 
 
 
 
 
 
 
 
 
# Bibliography
 
 * gnuplot book [https://livebook.manning.com/book/gnuplot-in-action-second-edition/index/]()
 
 * gnuplot examples [http://www.gnuplot.info]()
 
 * gnuplot examples [http://www.gnuplotting.org]()
 
 
 
 
 */
