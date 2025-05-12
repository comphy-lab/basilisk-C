reset
set terminal pngcairo size 1024,768 enhanced font 'Arial,18'
set xlabel "Wavenumber"
set ylabel "Growth rate"

do for [i=1:50] {
  stats 'amp-'.i.'.dat' every ::::0 using 2 nooutput
  k = STATS_min
  f(x) = a + b*x
  fit f(x) 'amp-'.i.'.dat' u 2:(log($4)) every ::1 via a,b
  set object circle at first k,b radius char 0.5 \
  fillstyle empty border lc rgb '#aa1100' lw 2
}
set samples 1000
set xrange [0:1.2]
plot real(sqrt(0.5*(x**2-x**4))) t 'Theory'
