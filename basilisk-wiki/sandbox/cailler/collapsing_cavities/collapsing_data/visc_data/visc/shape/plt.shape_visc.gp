set term qt enhanced persist font "Helvetica Neue, 14"

set size ratio 1

path = "visc-zenon_gas_refined-shape_s230_lvl12_t29.8_th120_m50_tinv15.dat"
path2 = "invisc-zenon-shape_s230_lvl12_t29.8_th120_m50_tinv15.dat"

set xrange [-40:50]
set yrange [0:90]

set xlabel 'z'
set ylabel 'r'

# set title "L_0 = 500 ; N_{max} = 9 ; θ_0 = 120° ; D = 10^{-3} ; t = 200 ; μ_0 = -4\n \
# AMR Erf func, RADIUS = 500"

p path u ($1):($2) w l lt -1 lw 1 t 'visc., t = 29.8', \
  path2 u ($1):($2) w l lw 1 t 'invisc, t = 29.8'