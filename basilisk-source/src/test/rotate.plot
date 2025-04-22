reset

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)

f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2

fc(x)=ac+bc*x
fit fc(x) 'clog' u (log($1)):(log($4)) via ac,bc
fc2(x)=ac2+bc2*x
fit fc2(x) 'clog' u (log($1)):(log($2)) via ac2,bc2

set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set key bottom left
set logscale
set xrange [16:256]
set xtics 16,2,256
set grid ytics
set cbrange [1:1]
plot 'log' u 1:4 t 'max (adaptive)', exp(f(log(x))) t ftitle(a,b), \
     'clog' u 1:4 t 'max (constant)', exp(fc(log(x))) t ftitle(ac,bc), \
     'log' u 1:2 t 'norm1 (adaptive)', exp(f2(log(x))) t ftitle(a2,b2), \
     'clog' u 1:2 t 'norm1 (constant)', exp(fc2(log(x))) t ftitle(ac2,bc2)

if (batch) set term @PNG; set output "error.png"; else pause -1;
reset
set pm3d
set pm3d map interpolate 1,1
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,\
     0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0,\
     0.875 0.9333 0 0, 1 0.498 0 0 )
set size ratio -1
splot [0.1:0.4][-0.15:0.15]'out' t ""

if (batch) set term @PNG; set output "interface.png"; else pause -1;
reset
set size ratio -1
plot [-0.5:0.5][-0.5:0.5]'interface' w l t ''
