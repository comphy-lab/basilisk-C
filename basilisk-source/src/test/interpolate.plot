reset
set title 'Interpolation on a regular grid'

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($2)) via a,b
set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set logscale
set cbrange [1:2]
set xrange [16:512]
set xtics 16,2,512
set grid ytics
plot 'log' u 1:2 t '', exp(f(log(x))) t ftitle(a,b)
