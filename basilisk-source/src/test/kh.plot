unset surface
set pm3d map
set contour base
set cntrparam levels 10
set cntrlabel onecolor
set zrange [-0.49999:0.49999]
unset ytics
unset key
splot 'out' u 1:2:3 w l lc rgbcolor "black"
