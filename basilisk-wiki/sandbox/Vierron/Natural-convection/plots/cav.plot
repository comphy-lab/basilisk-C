set terminal qt 0
reset
set xlabel 'log10(surface)'
set ylabel 'Number of droplets'
binstart = -100
binwidth = 0.2
set boxwidth binwidth
set style fill solid 0.5
Delta = 1./2**8 # minimum cell size
minvol = log10(Delta**2) # minimum cell surface
set arrow from minvol, graph 0 to minvol, graph 1 nohead
set label "Minimum cell volume" at minvol - 0.2, graph 0.35 rotate by 90
t = "300."
plot "< awk '{if ($2==".t.") print $0;}' out" u \
     (binwidth*(floor((log10($4)-binstart)/binwidth)+0.5)+binstart):(1.0) \
     smooth freq w boxes t '' 


set terminal qt 1
set xlabel 'temps t'
set ylabel 'Nu'
set title 'Nu=f(t)'
plot 'data' u 11:4 w l title 'Nu-top',\
      '' u 11:5 w l title 'Nu-bot'
