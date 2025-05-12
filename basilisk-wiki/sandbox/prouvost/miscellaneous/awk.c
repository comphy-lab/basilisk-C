
int main() {}

/**
I share here some useful examples of what we can do with awk when post-treating data.
You can do the same with other programming languages such as python, but awk is fun, it does the job in only one line and it's easy to write awk commands to use inside gnuplot plots.
It's also very efficient when used in bash programs.
Note that I don't pretend to write the most efficient awk programs.


# Elements of AWK

Awk is a program which reads an input file line by line and allows to perform operations based on the data written in these lines.

## Extract lines containing a particular value

Lets assume you made a parametric study/series of experiments and you registered in a single file all the results.
For example, in the first column you have the radius of a bubble, and in the second one the pressure in the bubble at a given time.
You can simply extract the results for a particular radius (lets says R=1) writing
*/
/*
awk '$1==1' data_file
*/
/**
where $1 describes the values in the first columns of the data file.

~~~gnuplot Example extract lines containing a particular value
set term pngcairo enhanced size 500,500
set output 'ex1.png'

p "./../data1" w p pt 7 ps 2 lc rgb 'blue' t 'original data',\
"< awk '$1==1' ./../data1" w p pt 7 lc rgb 'red' t 'data R=1'
~~~

*/


/**
## Extract max value

Lets assume you want to extract the maximum pressure value in the previously described experiment and the corresponding bubble radius.
Then,
*/
/*
awk '{ if ($2>tmp_p) {tmp_p=$2;tmp_R=$1} } END{print tmp_R, tmp_p}' data_file
 */
/**
will do the job. 
Note that the variables declared inside awk are automatically initialized with the value 0.

~~~gnuplot Example extract maximum pressure and corresponding radius
set term pngcairo enhanced size 500,500
set output 'ex2.png'

p "./../data1" w p pt 7 ps 2 lc rgb 'blue' t 'original data',\
"< awk '{ if ($2>tmp_p) {tmp_p=$2;tmp_R=$1} } END{print tmp_R, tmp_p}' ./../data1" w p pt 7 lc rgb 'red' t 'max pressure'
~~~

## Compute moving average

Lets assume you registered a noisy pressure signal over time and you want to filter it using a moving average on 100 points.
In the first column is the time and in the second the pressure.
You can use
*/
/*
awk -v n=100 '{sum -= mem[NR%n]; mem[NR%n]=$2; sum += $2} (NR>(n-1)){print $1,sum/n}' data_file
 */
/**
~~~gnuplot Example moving average
set term pngcairo enhanced size 500,500
set output 'ex3.png'

p "./../data2" u 2:5 w p lc rgb 'blue' t 'original data',\
"< awk -v n=100 '{sum -= mem[NR%n]; mem[NR%n]=$5; sum += $5} (NR>(n-1)){print $2,sum/n}' ./../data2" w l lw 2 lc rgb 'red' t 'moving average'
~~~


## Integration

Last example for now, lets assume you want to compute the impulse 
$$ I = \int_t (p-p_0) dt $$
with $p_0=1$
and you only have registered the time in the first column and the pressure in the second one. 
You can use a "trapeze" integral
 */
/*
awk  '{ if ($1!=tmpt) {s+=($2+tmpy)/2*($1-tmpt)} ; print $1,s ; tmpy=$2 ; tmpt=$1  }'  data_file
 */
/**
~~~gnuplot Example integrate
set term pngcairo enhanced size 500,500
set output 'ex4.png'

set ytics nomirror tc "blue"
set y2tics nomirror tc "red"
set ylabel "p" tc "blue"
set y2label "I = \int (p-1) dt" tc "red"

p "./../data3" u 1:3 w p lc rgb 'blue' axis x1y1 not,\
"< awk 'BEGIN{s=0}{ if ($1!=tmpt) {s+=(($3+tmpy)/2-1)*($1-tmpt)} ; print $1,s ; tmpy=$3 ; tmpt=$1 }' ./../data3" w p lc rgb 'red' axis x1y2 not
~~~
*/





/**

# Structure/options I often use

The option -v allows to use variables computed outside of the awk command.
The keyword BEGIN is for operation to do at the begining of the program, END at the end (after reading all lines), (NR>1) is for all line number strictly greater than 1.

*/

/*
awk -v var=42 -v var2=3 'BEGIN{...} (NR>1){...} END{...}' filename
*/

/**

Note: the index of the array defined in awk can be everything: you can even define non integer indices such as a[0.1] (ex: awk '{a[0.1]=1 ; print a[0.1]}' data1 works...)

*/

/**

# User-defined Functions

awk does not contain all "standard" mathematical functions as built-in functions. It contains [these ones](https://www.gnu.org/software/gawk/manual/html_node/Numeric-Functions.html).
The "missing" functions (such as arcsin, arctanh, ...) and variables ($\pi$) can be deduced and created from the built-in functions.

Example:

*/

/*

awk 'function asin(x) { return atan2(x, sqrt(1.-x*x)) } BEGIN {pi = atan2(0, -1) ; print pi,asin(pi)} '

*/

/**

arcsin : function asin(x) { return atan2(x, sqrt(1.-x*x)) }

arctanh : function atanh(x) {return log((1.+x)/(1.-x))/2}

$\pi$ : pi = atan2(0, -1)

Note: beware of the domain of definition, and problematic values with these definitions ($x$ close to 1) for arctanh for example...)

Note: It is possible to write awk scripts and not all the commands in one line

*/

/**

# Miscellaneous

## awk as a calculator in bash:
*/
/*
a=$(awk "BEGIN {toto = 1+1 ; print toto}")
*/


/**
## print a blank line

It may be usefull to reorganise a file to plot it with line-points with gnuplot 
*/
/*
awk '($1>=13){print ; c++ ;  if (!(c%3)) print Blank} '  data_file
*/











