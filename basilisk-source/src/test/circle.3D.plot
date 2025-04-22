reset
set title 'Poisson solution with a circular refined patch'

minlevel = 5
maxlevel = 7

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) '< grep "max error" log' u (log(2**$3)):(log($4)) via a,b
set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set logscale
set cbrange [1:2]
set xrange [2**(minlevel-1):2**(maxlevel+1)]
set xtics 2**(minlevel-1),2,2**(maxlevel+1)
set grid ytics
plot '< grep "max error" log' u (2**$3):4 t '', exp(f(log(x))) t ftitle(a,b)

if (batch) set term @PNG; set output "res.png"; else pause -1;
reset
set cbrange [1:2]
set title 'Poisson solution with a circular refined patch'
set xlabel 'Multigrid iteration'
set ylabel 'Residual'
set logscale y
set grid ytics
plot for [i = minlevel:maxlevel] \
     '< grep "residual '.i.'" log' u 3:4 w lp t 'level '.i

if (batch) set term @PNG; set output "speed.png"; else pause -1;
reset
set cbrange [1:2]
set title 'Poisson solution with a circular refined patch'
set xlabel 'CPU Time'
set ylabel 'Residual'
set logscale
plot for [i = minlevel:maxlevel] \
     '< grep "speed '.i.'" out' u 4:5 w lp t 'level '.i
