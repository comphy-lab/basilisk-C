/**
# Flashlube data

~~~gnuplot Flashlube verbruik met schoef verdraaingen. De grijze lijnen geven het beoogde verbruik weer (1mL op 10km).
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5
set pointintervalbox 3
set offsets graph 0.05, 100, 0 ,0
set xlabel 'ODO [km]'
set ylabel 'Flashlube volume [mL]'
set key off
plot 'out' u 1:($2)+15:3 w labels, 'out' u 1:2 w linespoints ls 1 ,\
  for [i=30:40] -0.1*x + 10*i + 68600*0.1 lt rgb 'grey'

~~~
*/

#define COLS (3)
#define ENTRIES (7)

int main() {
  double data[ENTRIES][COLS] = {{68598, 400,  0},
                                {68619, 350, -1},
                                {68652, 310, -1.5},
                                {69021, 285,  0},
                                {69852, 220,  0.25},
                                {70543, 165, 0},
                                {71417, 125, 0.5}
                               }; 
  printf ("#km\tV\tturns\n");
  for (int i = 0; i < ENTRIES; i++) { 
    for (int j = 0; j < COLS; j++) 
      printf ("%g\t", data[i][j]);
    printf ("\n");
  }
}
