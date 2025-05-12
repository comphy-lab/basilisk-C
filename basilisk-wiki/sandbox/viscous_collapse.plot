/** gnuplot */

!awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < viscous_collapse.log

set term png

set output viscous_collapse1.png

reset
set title 'Self Similar Solution for a viscous collapse, t=200, 400...1000'

 
set xlabel 'x'
set ylabel 'h(x t^(-1/3),t) t^(1/3)'

s=2.5
b=2.18812
h(x) = 0.543217  * (1-x*x/b/b)**(1./3)
bed(x) = 0
t5=200**.2

set key top right

plot [-3:3] \
      '< grep ^p viscous_collapse.out' u 2:3 w l lc 3 t 'Numerical', \
      h(x) lw 1 lc 1 lt 1 t 'self similar',\
      'viscous_collapse.ref' u ($1/t5):($2*t5) t'order 1, HLL t=200 N=512' w l  lc 2
      
      
# pause -1
set output viscous_collapse2.png
reset
set xlabel 'x'
set ylabel 'h(x,t)'
set arrow from 0,.27 to 0,0.08 
set label "t" at 0.2,0.1

set arrow from 3,.05 to 10,0.05 
set label "t" at 10,0.06 
 
plot [][:.3]'eta-0' not  w l lc -1,'eta-1'not  w l lc -1,\
    'eta-2'not  w l lc -1,'eta-3'not  w l lc -1,'eta-4'not  w l lc -1,'eta-5'not  w l lc -1
 
      
 