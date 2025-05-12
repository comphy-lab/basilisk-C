/**
# COVID-19 Data

From the [European Center for Disease Prevention and
Control](https://www.ecdc.europa.eu).

~~~gnuplot 
set term svg font ',11' size 800,600

# get the data if it's not already there
# note that ssconvert is required and can be installed with
# sudo apt install gnumeric

! test -f covid-`date +%Y-%m-%d`.xslx || \
  (wget https://www.ecdc.europa.eu/sites/default/files/documents/COVID-19-geographic-disbtribution-worldwide-`date +%Y-%m-%d`.xlsx \
   -O covid-`date +%Y-%m-%d`.xslx && \
   ssconvert covid-`date +%Y-%m-%d`.xslx covid.csv) || \
   rm -f covid-`date +%Y-%m-%d`.xslx

set key top left
set xlabel 'Days (after the hundredth case)'
set grid

# choose which countries to display
countries="United_States_of_America Italy Spain France Germany United_Kingdom China Japan South_Korea Netherlands Sweden Argentina Australia New_Zealand"

replaceu(x)=system("echo ".x."| sed 's/_/ /g'")

# Cases

item(n) = word(countries,n)
igap(n) = word(gap,n)


load "../average.gp"

set ylabel 'Number of cases, 7 days moving averaging'
plot [0:60] for [i=1:words(countries)] sum = 0 '< awk -f covid.awk < covid.csv' \
index item(i)  u ($1):(avg_n($3)) title replaceu(item(i)) w lp

~~~

~~~gnuplot
pop = "33000 6000 4700 6700 8300 6600 140000 12600 5100 1728 1023 4500 2546 485"
ipop(n) = word(pop,n)
set ylabel 'Number of cases, 7 days moving averaging / 10000 people'
plot [0:60][0:2] for [i=1:words(countries)] sum = 0 '< awk -f covid.awk < covid.csv' \
index item(i)  u ($1):(avg_n($3)/ipop(i)) title replaceu(item(i)) w lp, \
 sum = 0 '< awk -f covid.awk < covid.csv' index 'France' \
 u ($1):(avg_n($3)/ipop(4)) t "" lt -1 w l
~~~

~~~gnuplot
title = "Zoom"
set ylabel 'Number of cases, 7 days moving averaging / 10000 people'
plot [-5:60][0:0.5] for [i=1:words(countries)] sum = 0 '< awk -f covid.awk < covid.csv' \
index item(i)  u ($1):(avg_n($3)/ipop(i)) title replaceu(item(i)) w lp
~~~

~~~gnuplot

set xlabel 'Days (after the hundredth case)'
set ylabel 'Number of cases'

plot [0:50][0:10000] for [i=1:words(countries)] sum = 0 '< awk -f covid.awk < covid.csv' \
index item(i)  using ($1):(sum=sum+$3) title replaceu(item(i)) w lp
~~~

~~~gnuplot

set xlabel 'Days (after the hundredth case)'
set ylabel 'Number of deaths'

plot [0:30][0:600] for [i=1:words(countries)] sum = 0 '< awk -f covid.awk < covid.csv' \
index item(i)  using ($1):(sum=sum+$4) title replaceu(item(i)) w lp
~~~

How many days left in France?
~~~gnuplot

xp = 45
set xtics ("1/03" -4, "16/03 \n confinement" 11,"8/04" 34, "1/05" 57, "11/05" 67,"1/06" 88, "1/07" 118)
set arrow from xp,800 to xp,100
set xlabel 'Date'
set label strftime("%b %d %Y", time(0)-24*60*60)  at xp,1100 center
gap = "0 8 0 2"

plot [-10:160][:10000] for [i=2:4] 	\
     '< awk -f covid.awk < covid.csv' index item(i)		\
     using ($1-igap(i)):(avg_n($3)) w lp t replaceu(item(i))
     
set xtics auto
unset label
~~~



~~~gnuplot
alatina = "Brazil Chile Bolivia Argentina Peru Uruguay Paraguay"
set ylabel 'Cases'
set xlabel 'Cases from case =  100'
set log y
plot [0.1:40][100:25000] for [country in alatina] sum=0,		\
     '< awk -f covid.awk < covid.csv' index country		\
     using 1:(sum=sum+$3) w lp t replaceu(country)	
   
unset log y
~~~

The length is around 50 days, we are in the 27th-28th day, so
we have still

- 50 - 28 = 22 days
- security time (contagion) 10 - 15 days
- between 32 and 37 days starting from 10 April 2020...


A dummy program. Everything is done by gnuplot.
*/

int main() {
  fprintf (stderr, "Please ignore this error. It is just to force updates.\n");
  exit (1); // always fails (forces rerun)
}

