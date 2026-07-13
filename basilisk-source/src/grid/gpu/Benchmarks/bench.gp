set multiplot layout 1,2

set title "Speed in grid points x timesteps / second"
plot for [i=2:np] \
    sprintf("< awk -v findex=3 -v minlevel=7 -f advection.awk %s", input) \
    using i:xtic(1) title columnhead lt i

set title "Speedup relative to 8 x Intel Core i7"
plot for [i=3:np] \
    sprintf("< awk -v findex=3 -v minlevel=7 -f advection.awk %s", input) \
    using (column(i)/$2):xtic(1) title columnhead lt i

unset multiplot
