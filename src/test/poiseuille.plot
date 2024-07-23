ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)
f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2
set xlabel 'Resolution'
set ylabel 'Error'
set logscale
set xrange [4:128]
set cbrange [1:2]
set xtics 4,2,128
set grid ytics
set yrange [1e-5:]
plot 'log' u 1:4 t 'max', exp(f(log(x))) t ftitle(a,b), \
     'log' u 1:2 t 'norm1', exp(f2(log(x))) t ftitle(a2,b2)
