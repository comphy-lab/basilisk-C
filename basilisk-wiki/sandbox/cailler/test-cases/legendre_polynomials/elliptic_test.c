#include "elliptic.h"

/**
~~~gnuplot Legendre functions of fractional orders
reset session

set style line 1  lc rgb '#0025ad' lt 1 lw 3.5 # --- blue
set style line 2  lc rgb '#0042ad' lt 1 lw 3.5 #      .
set style line 3  lc rgb '#0060ad' lt 1 lw 3.5 #      .
set style line 4  lc rgb '#007cad' lt 1 lw 3.5 #      .
set style line 5  lc rgb '#0099ad' lt 1 lw 3.5 #      .
set style line 6  lc rgb '#00ada4' lt 1 lw 1.5 #      .
set style line 7  lc rgb '#00ad88' lt 1 lw 1.5 #      .
set style line 8  lc rgb '#00ad6b' lt 1 lw 1.5 #      .
set style line 9  lc rgb '#00ad4e' lt 1 lw 1.5 #      .
set style line 10 lc rgb '#00ad31' lt 1 lw 1.5 #      .
set style line 11 lc rgb '#00ad14' lt 1 lw 1.5 #      .
set style line 12 lc rgb '#09ad00' lt 1 lw 1.5 # --- green

set term pngcairo enhanced size 1000,500
set output 'legendre_half.png'

set size ratio 1
set multiplot layout 1,2

# Plot Legendre P polynomials
set title "Legendre functions of the first kind"
set xlabel "x"  
set xrange [0:pi]
set yrange [-2.3:2.3]
# set xtics (0, 'π/2' pi/2. ,'π' pi)
set xzeroaxis

p "legendre_half_order_profile" u 1:2 w l ls 3 t "P_1_/_2", \
  "legendre_three_half_order_profile" u 1:2 w l ls 6 t "P_3_/_2", \
  "legendre_five_half_order_profile" u 1:2 w l ls 9 t "P_5_/_2"

# Plot Legendre Q polynomials
set title "Legendre functions of the second kind"
set xlabel "x"  
set xrange [0:pi]
set yrange [-2.3:2.3]
# set xtics (0, 'π/2' pi/2. ,'π' pi)
set xzeroaxis

p "legendreQ_half_order_profile" u 1:2 w l ls 3 t "Q_1_/_2", \
  "legendreQ_three_half_order_profile" u 1:2 w l ls 6 t "Q_3_/_2", \
  "legendreQ_five_half_order_profile" u 1:2 w l ls 9 t "Q_5_/_2"

unset multiplot
~~~
*/



int main(){
  char filename1[80];
  char filename2[80];
  char filename3[80];
  char filename4[80];
  char filename5[80];
  char filename6[80];
  sprintf(filename1, "legendre_half_order_profile");
  sprintf(filename2, "legendreQ_half_order_profile");
  sprintf(filename3, "legendre_three_half_order_profile");
  sprintf(filename4, "legendreQ_three_half_order_profile");
  sprintf(filename5, "legendre_five_half_order_profile");
  sprintf(filename6, "legendreQ_five_half_order_profile");
  FILE * fp1 = fopen(filename1, "w");
  FILE * fp2 = fopen(filename2, "w");
  FILE * fp3 = fopen(filename3, "w");
  FILE * fp4 = fopen(filename4, "w");
  FILE * fp5 = fopen(filename5, "w");
  FILE * fp6 = fopen(filename6, "w");
  for (double x = 0.; x <= pi; x += 0.01){
    double P_half = legendreP_half(x);
    double Q_half = legendreQ_half(x);
    double P_three_half = legendreP_three_half(x);
    double Q_three_half = legendreQ_three_half(x);
    double P_five_half = legendreP_five_half(x);
    double Q_five_half = legendreQ_five_half(x);
    fprintf (fp1, "%g %g \n", x, P_half);
    fprintf (fp2, "%g %g \n", x, Q_half);
    fprintf (fp3, "%g %g \n", x, P_three_half);
    fprintf (fp4, "%g %g \n", x, Q_three_half);
    fprintf (fp5, "%g %g \n", x, P_five_half);
    fprintf (fp6, "%g %g \n", x, Q_five_half);
  }
  fclose (fp1);
  fclose (fp2);
  fclose (fp3);
  fclose (fp4);
  fclose (fp5);
  fclose (fp6);
}
