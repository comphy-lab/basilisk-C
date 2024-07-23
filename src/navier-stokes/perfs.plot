reset
unset mouse

mpl_top    = 0.4 #inch  outer top margin, title goes here
mpl_bot    = 1.0 #inch  outer bottom margin, x label goes here
mpl_left   = 1.0 #inch  outer left margin, y label goes here
mpl_right  = 1.0 #inch  outer right margin, y2 label goes here
mpl_height = 3.0 #inch  height of individual plots
mpl_width  = 3.0 #inch  width of individual plots
mpl_dx     = 0.2 #inch  inter-plot horizontal spacing
mpl_dy     = 0.4 #inch  inter-plot vertical spacing
mpl_ny     = 5   #number of rows
mpl_nx     = 2   #number of columns

# calculate full dimensions
xsize = mpl_left+mpl_right+(mpl_width*mpl_nx)+(mpl_nx-1)*mpl_dx
ysize = mpl_top+mpl_bot+(mpl_ny*mpl_height)+(mpl_ny-1)*mpl_dy

# placement functions
#   rows are numbered from bottom to top
bot(n) = (mpl_bot+(n-1)*mpl_height+(n-1)*mpl_dy)/ysize
top(n)  = 1-((mpl_top+(mpl_ny-n)*(mpl_height+mpl_dy))/ysize)
#   columns are numbered from left to right
left(n) = (mpl_left+(n-1)*mpl_width+(n-1)*mpl_dx)/xsize
right(n)  = 1-((mpl_right+(mpl_nx-n)*(mpl_width+mpl_dx))/xsize)

set offsets
set autoscale fix
set size 1,1

# define x-axis settings for all subplots
set xlabel ''
set format x ''

# start plotting
set multiplot

stats "perfs" u 1:9 nooutput
EVERY = ceil(STATS_records/200.)
unset key
set style fill solid

#-----------------------------------------------
# subplot  1-4
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for third row (top)
set tmargin at screen top(5)
set bmargin at screen bot(5)

set ylabel "dt"
set logscale y
plot 'perfs' u 1:2 every EVERY w boxes lc 0
unset logscale y

#-----------------------------------------------
# subplot  2-4
#  set horizontal margins for second column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for third row (top)
set tmargin at screen top(5)
set bmargin at screen bot(5)

if (STATS_min_y == STATS_max_y) set y2range [STATS_max_y-1:STATS_max_y+1]; \
else set y2range [STATS_min_y:STATS_max_y];
unset ytics
set y2tics
unset ylabel
set y2label '# cells'
plot 'perfs' u 1:9 every EVERY w boxes lc 0
set yrange [*:*]

#-----------------------------------------------
# subplot  1-3
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for second row (middle)
set tmargin at screen top(4)
set bmargin at screen bot(4)

unset ylabel
unset y2label
set ytics auto
unset y2tics
colori(i) = (i < 10 ? 65280 : i > 20 ? 16711680 : 16753920)
set ylabel 'mgp.i'
plot [][0:]'perfs' u 1:3:(colori($3)) every EVERY w boxes lc rgbcolor variable

#-----------------------------------------------
# subplot  2-3
#  set horizontal margins for second column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for second row (middle)
set tmargin at screen top(4)
set bmargin at screen bot(4)

stats "perfs" u 1:4 nooutput
if (STATS_min_y == STATS_max_y) set y2range [STATS_max_y-1:STATS_max_y+1]; \
else set y2range [STATS_min_y:STATS_max_y];
unset ytics
set y2tics auto
unset ylabel
set y2label 'mgp.nrelax'
plot [][0:]'perfs' u 1:4:(colori($4)) every EVERY w boxes lc rgbcolor variable

#-----------------------------------------------
# subplot  1-2
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for second row (middle)
set tmargin at screen top(3)
set bmargin at screen bot(3)

stats "perfs" u 1:5 nooutput
if (STATS_min_y == STATS_max_y) set yrange [0:STATS_max_y+1]; \
else set yrange [STATS_min_y:STATS_max_y];
unset ylabel
unset y2label
set ytics auto
unset y2tics
set ylabel 'mgpf.i'
plot [][0:]'perfs' u 1:5:(colori($5)) every EVERY w boxes lc rgbcolor variable

#-----------------------------------------------
# subplot  2-2
#  set horizontal margins for second column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for second row (middle)
set tmargin at screen top(3)
set bmargin at screen bot(3)

stats "perfs" u 1:6 nooutput
if (STATS_min_y == STATS_max_y) set y2range [0:STATS_max_y+1]; \
else set y2range [STATS_min_y:STATS_max_y];
unset ytics
set y2tics auto
unset ylabel
set y2label 'mgpf.nrelax'
plot [][0:]'perfs' u 1:6:(colori($6)) every EVERY w boxes lc rgbcolor variable
set yrange [*:*]

#-----------------------------------------------
# subplot  1-1
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for second row (middle)
set tmargin at screen top(2)
set bmargin at screen bot(2)

stats "perfs" u 1:7 nooutput
if (STATS_min_y == STATS_max_y) set yrange [0:STATS_max_y+1]; \
else set y2range [STATS_min_y:STATS_max_y];
unset ylabel
unset y2label
set ytics auto
unset y2tics
set ylabel 'mgu.i'
plot [][0:]'perfs' u 1:7:(colori($7)) every EVERY w boxes lc rgbcolor variable
set yrange [*:*]

#-----------------------------------------------
# subplot  2-1
#  set horizontal margins for second column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for second row (middle)
set tmargin at screen top(2)
set bmargin at screen bot(2)

stats "perfs" u 1:8 nooutput
if (STATS_min_y == STATS_max_y) set y2range [0:STATS_max_y+1]; \
else set y2range [STATS_min_y:STATS_max_y];
unset ytics
set y2tics auto
unset ylabel
set y2label 'mgu.nrelax'
plot [][0:]'perfs' u 1:8:(colori($8)) every EVERY w boxes lc rgbcolor variable
set yrange [*:*]

#-----------------------------------------------
# subplot  1-0
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for first row (bottom)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set xlabel "time"
unset y2label
set ylabel "Wall-clock time"
set xtics auto
set format x '% g'
unset y2tics
set ytics auto
plot 'perfs' u 1:10 every EVERY w boxes lc 0

#-----------------------------------------------
# subplot  2-0
#  set horizontal margins for second column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for first row (bottom)
set tmargin at screen top(1)
set bmargin at screen bot(1)

stats "perfs" u 1:($11/$12) nooutput
if (STATS_min_y == STATS_max_y) set y2range [STATS_max_y-1:STATS_max_y+1]; \
else set y2range [STATS_min_y:STATS_max_y];
unset ytics
set y2tics auto
unset ylabel
set y2label 'points.step/sec/core'
set format y '%e'
plot 'perfs' u 1:($11/$12) every EVERY w boxes lc 0

unset multiplot

pause 10
reread
