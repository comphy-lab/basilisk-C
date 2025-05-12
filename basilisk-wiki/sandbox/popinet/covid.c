/**
# COVID-19 Data

From the [European Center for Disease Prevention and
Control](https://www.ecdc.europa.eu).

~~~gnuplot Number of cases. The grey curve $f(x)$ is a fit of a [sigmoid function](https://en.wikipedia.org/wiki/Sigmoid_function).
set term svg font ',11' size 800,600

# get the data if it's not already there
# note that ssconvert is required and can be installed with
# sudo apt install gnumeric

! test -f covid-`date +%Y-%m-%d`.xslx || \
  (wget https://www.ecdc.europa.eu/sites/default/files/documents/COVID-19-geographic-disbtribution-worldwide-`date +%Y-%m-%d`.xlsx \
   -O covid-`date +%Y-%m-%d`.xslx && \
   ssconvert covid-`date +%Y-%m-%d`.xslx covid.csv) || \
   rm -f covid-`date +%Y-%m-%d`.xslx

set key above
set xlabel 'Days (after the hundredth case)'
set ylabel 'Number of cases'
set grid

# choose which countries to display
countries="United_States_of_America Italy Spain France Germany United_Kingdom China Japan Brazil Sweden Argentina New_Zealand Australia India"

replaceu(x)=system("echo ".x."| sed 's/_/ /g'")

# Try a sigmoid function
f(x)=a*exp((x-c)/b)/(exp((x-c)/b) + 1.)
a = 50000.
b = 5.
c = 20.
sum = 0.
fit [0:]f(x) '< awk -f covid.awk < covid.csv' index 'Italy' \
  using ($1):(sum=sum+$3) via a,b,c

set label 1 at 17,50000 "double every 2 days" rotate by 57
set label 2 at 30,80000" double every 3 days" rotate by 45
set label 3 at 30,21000" double every 4 days" rotate by 37

tend = 300

set logscale y
plot [0:tend][100:1e7]for [country in countries] sum=0,		\
     '< awk -f covid.awk < covid.csv' index country		\
     using 1:(sum=sum+$3) w lp t replaceu(country),		\
     100.*exp(log(2)/2.*x)  not w l linec -1,                   \
     100.*exp(log(2)/3.*x)   not w l linec -1,  		\
     100.*exp(log(2)/4.*x) not w l linec -1,                    \
     f(x) lc rgb 'grey'
unset label 1
unset label 2
unset label 3
~~~

~~~gnuplot Relative number of cases
set ylabel 'Relative number of cases (%)'
plot [0:tend][1e-2:] for [country in countries] sum=0,	       \
      '< awk -f covid.awk < covid.csv' index country	       \
      using 1:((sum=sum+$3)/$5*100.) w lp t replaceu(country)
~~~

~~~gnuplot Number of deaths
g(x)=a1*exp((x-c1)/b1)/(exp((x-c1)/b1) + 1.)
a1 = 50000.
b1 = 5.
c1 = 20.
sum = 0.
fit [0:]g(x) '< awk -f covid.awk < covid.csv' index 'Italy' \
  using 1:(sum=sum+$4) via a1,b1,c1

set ylabel 'Number of deaths'
plot [0:tend][1:3e5] for [country in countries] sum=0,	       \
      '< awk -f covid.awk < covid.csv' index country	       \
      using 1:(sum=sum+$4) w lp t replaceu(country),                     \
      g(x) lc rgb 'grey'
~~~

~~~gnuplot Relative number of deaths
set ylabel 'Relative number of deaths (%)'
plot [0:tend][1e-4:] for [country in countries] sum=0,	       \
      '< awk -f covid.awk < covid.csv' index country	       \
      using 1:((sum=sum+$4)/$5*100.) w lp t replaceu(country)
~~~

~~~gnuplot Mortality (%)
unset logscale
set ylabel 'Mortality (%)'
plot [0:tend][0:20]for [country in countries] sum=0, sum1=0,		\
       '< awk -f covid.awk < covid.csv' index country		        \
       using 1:(100.*(sum1=sum1+$4)/(sum=sum+$3)) w lp t replaceu(country),       \
       100.*g(x)/f(x) lc rgb 'grey'
~~~

~~~gnuplot New cases per day
set ylabel 'New cases / day'
df(x)=(a*exp((x-c)/b))/(b*(exp((x-c)/b)+1))-(a*exp((2*(x-c))/b))/(b*(exp((x-c)/b)+1)**2)
set logscale y
plot [0:tend][10:1e5]for [country in countries] '< awk -f covid.awk < covid.csv' \
    index country using 1:3 w p t replaceu(country), \
    df(x) lc rgb 'grey'
~~~

~~~gnuplot Relative number of new cases per day
set ylabel 'Relative # new cases / day (%)'
plot [0:tend][1e-4:]for [country in countries] '< awk -f covid.awk < covid.csv' \
    index country using 1:($3/$5*100.) w p t replaceu(country)
~~~

~~~gnuplot Deaths per day
set ylabel 'Deaths / day'
dg(x)=(a1*exp((x-c1)/b1))/(b1*(exp((x-c1)/b1)+1))-(a1*exp((2*(x-c1))/b1))/(b1*(exp((x-c1)/b1)+1)**2)
plot [0:tend][1:]for [country in countries] '< awk -f covid.awk < covid.csv' \
    index country using 1:4 w p t replaceu(country), \
    dg(x) lc rgb 'grey'
~~~

~~~gnuplot Relative # deaths per day
set ylabel 'Relative # deaths / day (%)'
plot [0:tend][1e-6:]for [country in countries] '< awk -f covid.awk < covid.csv' \
    index country using 1:($4/$5*100.) w p t replaceu(country)
~~~

~~~gnuplot New cases against total nr. of cases
set ylabel 'New cases / day'
set xlabel 'Number of Cases'
set logscale xy
plot [10:1e7][1:100000]                                 \
     for [country in countries] sum = 0, '< awk -f covid.awk < covid.csv' \
    index country using (sum=sum+$3):3 w p t replaceu(country), \
    x/5
~~~

If we assume an exponential growth $N(t)\propto e^{t/b}$, the ratio of the total number of cases to the number of new cases per day is
$$
N(t)/\frac{dN(t)}{dt} = b
$$
which gives roughly $b=5$ days, according to the graph below. This means an average "doubling time" of $\log(2)\times 5\approx 3.5$ days.

~~~gnuplot Total number of cases / New cases
set ylabel 'Total number of cases / New cases'
set xlabel 'Number of Cases'
unset logscale
set logscale xy
plot [10:1e7][:]                                 \
     for [country in countries] sum = 0, sum1 = 0, '< awk -f covid.awk < covid.csv' \
    index country using (sum=sum+$3):((sum1=sum1+$3)/$3) w p t replaceu(country), 5
~~~

## Rescaling with a sigmoid

The [sigmoid](https://en.wikipedia.org/wiki/Sigmoid_function) or [logistic](https://en.wikipedia.org/wiki/Logistic_function) function used to fit all data is
$$
f(t) = \frac{N}{1 + e^{- (t - t_d)/t_g}}
$$
where $N$ is the maximum, $t_g$ is the "growth timescale" and $t_d$ is
(half) the total duration.

Note that the choice of this function can be justified based on the same arguments used by [Verhulst](https://en.wikipedia.org/wiki/Pierre_Fran%C3%A7ois_Verhulst) who first introduced it as a model of population growth (in particular to counter the arguments advanced by [Malthus](https://en.wikipedia.org/wiki/Thomas_Robert_Malthus)): i.e. the initial exponential growth is eventually bounded by a "crowding effect". In the case of epidemics, this can be interpreted as the increasing difficulty of the pathogen to spread as the density of already infected individuals increases.

Note that this function 
$f(t)$ is solution of the following differential equation :
$$
\frac{d f(t)}{dt}= \frac{f(t)}{t_g} (1 - \frac{f(t)}{N})
$$

The fit for the total number of cases for all countries is convincing,
although beware of the ["von Neuman fitting
effect"](https://en.wikiquote.org/wiki/John_von_Neumann): *With four
parameters I can fit an elephant, and with five I can make him wiggle
his trunk.*

~~~gnuplot Rescaled number of cases. The grey curve is the [sigmoid function]().

nc = 0; do for [c in countries] { nc = nc + 1; }

array country[nc]
array N[nc]
array Ne[nc]
array tg[nc]
array td[nc]
array tde[nc]

i = 1; do for [c in countries] { country[i] = c; i = i + 1; }

# Try a sigmoid function
f(x,a,b,c)=a*exp((x-c)/b)/(exp((x-c)/b) + 1.)
do for [i=1:nc] {
  a = 500000.
  b = 4.5
  c = 100.
  sum = 0.
  set fit errorvariables
  set fit prescale
  print country[i]
  fit [0:]f(x,a,b,c) '< awk -f covid.awk < covid.csv' index country[i] \
     using ($1):(sum=sum+$3) via a,b,c
  N[i] = a
  Ne[i] = a_err/a
  tg[i] = b
  td[i] = c
  tde[i] = c_err/c
}

set xlabel '(t - t_d)/t_g'
set ylabel 'Number of cases / N'
set grid
unset logscale

plot [-5:10][*:*]							\
     for [i=1:nc] sum=0,					        \
     '< awk -f covid.awk < covid.csv' index country[i]			\
     using (($1-td[i])/tg[i]):((sum=sum+$3)/N[i]) w p t replaceu(country[i]), \
     f(x,1.,1.,0.) lc rgb 'grey' t 'sigmoid'
~~~

And also for the number of deaths (with different fitting parameters).

~~~gnuplot Rescaled number of deaths
array N1[nc]
array N1e[nc]
array t1g[nc]
array t1d[nc]

do for [i=1:nc] {
  a = 100000.
  b = 5.
  c = 80.
  sum = 0.
  print country[i]
  fit [0:]f(x,a,b,c) '< awk -f covid.awk < covid.csv' index country[i]	\
  using 1:(sum=sum+$4) via a,b,c
  N1[i] = a
  N1e[i] = a_err/a
  t1g[i] = b
  t1d[i] = c
}

set xlabel '(t - t1_d)/t1_g'
set ylabel 'Number of deaths / N1'
plot [-5:10][*:*]							\
     for [i=1:nc] sum=0,					        \
     '< awk -f covid.awk < covid.csv' index country[i]			\
     using (($1-t1d[i])/t1g[i]):((sum=sum+$4)/N1[i]) w p t replaceu(country[i]), \
        f(x,1.,1.,0.) lc rgb 'grey' t 'sigmoid'
~~~

Things are more noisy for the slopes.

~~~gnuplot Rescaled new cases per day
phi(x)=1./sqrt(2.*pi)*exp(-x*x/2.)
Phi(x)=(1.+erf(x/sqrt(2.)))/2.

set xlabel '(t - t_d)/t_g'
set ylabel 'New cases per day x t_g / N_{max}'
df(x,a,b,c)=(a*exp((x-c)/b))/(b*(exp((x-c)/b)+1))-(a*exp((2*(x-c))/b))/(b*(exp((x-c)/b)+1)**2)
set logscale y
plot [-6:10][1e-3:1]							\
  for [i=1:nc]								\
     '< awk -f covid.awk < covid.csv' index country[i]			\
     using (($1-td[i])/tg[i]):($3*tg[i]/N[i]) w p t replaceu(country[i]), \
     df(x,1.,1.,0.) lc rgb 'grey' t 'derivative of sigmoid'
~~~

~~~gnuplot Rescaled deaths per day
set xlabel '(t - t1_d)/t1_g'
set ylabel 'Deaths per day x t1_g / N1_{max}'
plot [-6:10][1e-3:1]							\
     for [i=1:nc]							    \
     '< awk -f covid.awk < covid.csv' index country[i]		 	    \
     using (($1-t1d[i])/t1g[i]):($4*t1g[i]/N1[i]) w p t replaceu(country[i]), \
     df(x,1.,1.,0.) lc rgb 'grey' t 'derivative of sigmoid'
~~~

The evolutions for each country can be summarised by these three
parameters, as represented below. Only countries for which the fitting error on $N$ or $N1$ is smaller than 20% are represented.

~~~gnuplot Fitting parameters for cases. The percentages give the relative fitting error on $N$.
set xlabel 't_d x 2 (days)'
set ylabel 't_g x log(2) (days)'

# select only countries for which the fitting error is less than 20%
select(i,j) = Ne[i] < 0.2 ? j : 1e100

unset key
set xrange [35:120]
set yrange [2:10]
unset logscale
scale=1000.
plot									\
"< echo '40 4.5 100000'" u 1:2:(sqrt($3/scale)) pt 7 ps variable lc 12,	      \
"< echo '40 4.5'" u 1:2:("100,000") w labels,				\
td using (select($1,$2*2)):(tg[$1]*log(2)):(sqrt(N[$1]/scale)):1	    \
pt 7 ps variable lc variable,						\
td using (select($1,$2*2)):(tg[$1]*log(2)):(sprintf("%s\n%.1f\%", replaceu(country[$1]), Ne[$1]*100.)) w labels
~~~

~~~gnuplot Fitting parameters for deaths. The percentages give the relative fitting error on $N$.
set xlabel 't1_d x 2 (days)'
set ylabel 't1_g x log(2) (days)'

# select only countries for which the fitting error is less than 20%
select1(i,j) = N1e[i] < 0.2 ? j : 1e100
scale=100.
plot									\
  "< echo '40 4.5 10000'" u 1:2:(sqrt($3/scale)) pt 7 ps variable lc 12,      \
  "< echo '40 4.5'" u 1:2:("10,000") w labels,				\
  t1d using (select1($1,$2*2)):(t1g[$1]*log(2)):(sqrt(N1[$1]/scale)):1	      \
     pt 7 ps variable lc variable,					\
  t1d using (select1($1,$2*2)):(t1g[$1]*log(2)):(sprintf("%s\n%.1f\%", replaceu(country[$1]), N1e[$1]*100.)) w labels
~~~

A dummy program. Everything is done by gnuplot.
*/

int main() {
  fprintf (stderr, "Please ignore this error. It is just to force updates.\n");
  exit (1); // always fails (forces rerun)
}
