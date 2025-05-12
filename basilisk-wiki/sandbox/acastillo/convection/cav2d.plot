set macros
bin='binary u 2:1:3'
color=' binary u 2:1:3 w ima notitle'
lines=' w l notitle'

file1="./cav2d_Psi.bin"    ; file2="./cav2d_T.bin"

set table 'table1.dat'
set contour ; unset surface ; set view map ; set cntrparam levels 20
sp file1 @bin
unset table

set table 'table2.dat'
set contour ; unset surface ; set view map ; set cntrparam levels 20
sp file2 @bin
unset table

set terminal pngcairo size 1000,360
set output "cav2d.png"
set xtics 0.5
set ytics 0.5
set xlabel "x"
set ylabel "y"
set view map
set size ratio 1
set multiplot layout 1,3
set title "Stream-function field"
set cbtics 0.01
p file1 @color, 'table1.dat' @lines
set palette defined (-1 "#F1A340", 0 "white", 1 "#998ec3")
set title "Temperature field"
set cbtics 0.2
set cbrange [-0.5:0.5]
p file2 @color, 'table2.dat' @lines

set title "Numerical convergence"
set size ratio 0.75
set xlabel "h" ; set xtics auto ; set logscale x
set ylabel "%Error(Nu)" ; set ytics auto ; set logscale y
set key bottom right
p "cav2d.asc" u (1./$3):(($4-1.118)/1.118*100.) title "Left", "" u (1./$3):(($5-1.118)/1.118*100.) title "Right"
