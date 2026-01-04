set xlabel 'y'
plot 'log' u 2:3 t "", \
     -0.51*(x - 0.5) + x**2/2. - 0.125 t "hydrostatic pressure", \
     0.01 + 0.5 - x t "density"

set output 'velocity.png'

reset
set size ratio -1
unset key
s=100
plot 'log' u 1:2:($4*s):($5*s) w vec, 'out' w l
