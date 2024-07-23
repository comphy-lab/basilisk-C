# generate results for Curie
cd 'curie'

# generate weak scaling curves
! bash weak.sh > weak

# change default color for line style 6 from "yellow" to "sea-green"
set style line 6 lw 1 lc rgb "sea-green"
set style increment user

set logscale
set grid
set xrange [2:32768]
set xtics 2

set xlabel "# of cores"
set ylabel "Memory/core (GB)"
set output 'memory.png'
minlevel=10
maxlevel=15
plot [][0.01:] for [i=minlevel:maxlevel] \
     '< sh table.sh poisson '.i u 1:($2/$1) t ''.i.' levels' w lp, \
     18/x**0.9

set ylabel 'Time (sec)'

set output 'poisson.png'
plot [][:100] for [i=minlevel:maxlevel] \
     '< sh time.sh poisson '.i u 1:2 t ''.i.' levels' w lp, \
     'weak' u 1:2 w lp t 'weak scaling', \
     600/x**0.95
     
set output 'poisson-mpi.png'
plot [][0.001:10] for [i=minlevel:maxlevel] \
     '< sh time.sh poisson '.i u 1:3 w lp t ''.i.' levels', \
     4.5/x**0.65

set output 'laplacian.png'
plot [][:10] for [i=minlevel:maxlevel] \
     '< sh time.sh laplacian '.i u 1:2 w lp t ''.i.' levels', \
     50/x**0.93

set output 'laplacian-mpi.png'
plot [][1e-5:1] for [i=minlevel:maxlevel] \
     '< sh time.sh laplacian '.i u 1:3 w lp t ''.i.' levels', \
     2./x**0.7

set output 'restriction.png'
plot [][:1] for [i=minlevel:maxlevel] \
     '< sh time.sh restriction '.i u 1:2 w lp t ''.i.' levels', \
     18/x**0.85

set output 'restriction-mpi.png'
plot [][1e-4:1] for [i=minlevel:maxlevel] \
     '< sh time.sh restriction '.i u 1:3 w lp t ''.i.' levels', \
     2.8/x**0.66
