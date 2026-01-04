reset
set title 'Poisson solution with a circular refined patch'

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) '< grep "max error" log' u (log(2**$3)):(log($4)) via a,b
set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set logscale
set cbrange [1:2]
set xrange [64:2048]
set xtics 64,2,2048
set grid ytics
plot '< grep "max error" log' u (2**$3):4 t '', exp(f(log(x))) t ftitle(a,b)

if (batch) set term @PNG; set output "res.png"; else pause -1;
reset
set title 'Poisson solution with a circular refined patch'
set xlabel 'Multigrid iteration'
set ylabel 'Residual'
set logscale y
set grid ytics
plot '< grep "residual 7" log' u 3:4 w lp t 'level 7', \
     '< grep "residual 8" log' u 3:4 w lp t 'level 8', \
     '< grep "residual 9" log' u 3:4 w lp t 'level 9', \
     '< grep "residual 10" log' u 3:4 w lp t 'level 10'

if (batch) set term @PNG; set output "speed.png"; else pause -1;
reset
set title 'Poisson solution with a circular refined patch'
set xlabel 'CPU Time'
set ylabel 'Residual'
set logscale
plot '< grep "speed 7" out' u 4:5 w lp t 'level 7', \
     '< grep "speed 8" out' u 4:5 w lp t 'level 8', \
     '< grep "speed 9" out' u 4:5 w lp t 'level 9', \
     '< grep "speed 10" out' u 4:5 w lp t 'level 10'
