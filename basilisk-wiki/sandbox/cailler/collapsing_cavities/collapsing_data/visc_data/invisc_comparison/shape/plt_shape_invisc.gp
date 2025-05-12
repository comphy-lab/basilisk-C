set term qt enhanced persist font "Helvetica Neue, 14"

set size ratio 1

path = "invisc-zenon-shape_s230_lvl12_t35_th120_m50_tinv15.dat"

set xrange [-0.4*230:0.6*230]
set yrange [0:230]

set xlabel 'z'
set ylabel 'r'

# set title "L_0 = 500 ; N_{max} = 9 ; θ_0 = 120° ; D = 10^{-3} ; t = 200 ; μ_0 = -4\n \
# AMR Erf func, RADIUS = 500"

p path u ($1):($2) w l lt -1 lw 1 notitle