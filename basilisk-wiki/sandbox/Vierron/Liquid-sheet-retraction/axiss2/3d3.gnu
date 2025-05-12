#!/usr/bin/gnuplot
### rotation of a curve
#Oh=0.1
#pas = 1./10
#tend = 30
filename(j) =  sprintf("allaxi%g%g.dat",j*pas,Oh)
stats filename((tend-pas)/pas) nooutput
zmax = STATS_max_x
do for [j=0:1/pas*tend]{
	reset
	set print $stat
	stats filename(j) nooutput
	print sprintf("%g %g", STATS_max_y, STATS_pos_max_y)
	print sprintf("%g %g", STATS_max_y, -STATS_pos_max_y)
	unset print
	set angle degrees
	set table $Rotation
	    array A[1]   # dummy array for plotting a single empty line
	    do for [i=0:360:10] {
		plot filename(j) u ($2*cos(i)):($2*sin(i)):1 w table
		plot filename(j) u ($2*cos(i)):($2*sin(i)):(-$1) w table
		plot $stat u ($1*cos(i)):($1*sin(i)):2 w table
		plot A u ("") w table
	    }
	unset table
	unset key
	set term pngcairo size 1024, 768
	#set output sprintf("liquid%g.png",j+1) uncomment without pipe
	set border
	set title sprintf("t=%g, Oh=%g, R0/H=10",j*pas,Oh)
	#set style fill solid 1.00 border -1
	set view 62, 8, 1.245, 1.0
	set palette defined (0 'yellow',1 'orange', 2 'red',3 'blue', 4 'navy',5 'blue', 6 'red',7 'orange', 8 'yellow')
	#set palette defined (0 'navy',1 'blue', 2 'red',3 'orange', 4 'yellow') non-symetrique
	set zrange [-3.5:3.5]
	#set cbrange [0:zmax] non-symetrique
	set cbrange [-zmax:zmax]
	set cblabel "Thickness"
	set ticslevel 0
	set pm3d depthorder interpolate 4,4 lighting specular 0.6 at s
	splot $Rotation u 1:2:3 w pm3d lt -2
}
# Video parameters
#system("ffmpeg -r 15 -f image2 -i liquid%d.png -c:v libx264 -y movie.mp4")

### end of code
