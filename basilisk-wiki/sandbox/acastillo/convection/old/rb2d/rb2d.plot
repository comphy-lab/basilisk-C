
set terminal pngcairo size 500,300
set output "rb2d.png"
set xlabel "Time"
set ylabel "Energy"

p "rb2d.asc" u 1:2 w l title "Kinetic", "" u 1:3  w l title "Potential", "" u 1:($3-$4)  w l title "Available potential"
