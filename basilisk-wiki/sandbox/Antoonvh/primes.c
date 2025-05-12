/** 
![A $4^2$ Morton curve. Basilisk employs a diagnoally mirrored version, i.e. (capital) N-order curve to be precise. Image courtesy of [Asger Hoedt](http://asgerhoedt.dk/?p=276)](http://asgerhoedt.dk/wp-content/uploads/2012/10/MortonCurve-4x4.png)

# Prime numbers along a Z-order space-filling curve
[Stanislaw Ulam](https://en.wikipedia.org/wiki/Stanislaw_Ulam) once dicided that is was a good idea to order the prime numbers (primes) along a spiralling space-filling curve. Doing so, he discovered what is now known as the [*Ulam's spiral*](https://www.youtube.com/watch?v=iFuR97YcSLM). It is attractive to use the Basilisk toolbox to study the behaviour of primes along a Z-order space-filling curve. If you are curious what role this curve plays within Baslisk, have a look [here](The_Tree_Grid_Structure_in_Basilsik#indexing-and-mpi-load-balancing).   
*/

#include "grid/quadtree.h" //<- For it's 2D Z-order indexing iterator
#include "utils.h" //<- For visualization purposes
#include "tag.h" //<- For finding connected regions. 
int i=9;
/**
## Find primes
We need to find all primes upto *n* and store them in an array *b* of length *n*. This is done by using [Eratosthenes' sieve](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes), an ancient algorithm that was not developed with computational efficientcy nor paralellization in mind. 
*/
void getprimes(int b[],int n){
  b[0]=0; //zero is not a prime
  b[1]=0; //one is not a prime
  for (int j=2;j<n;j++){
    b[j]=j;
  }
  int j=2;
  while(j<=ceil(sqrt(n))){
    int num = b[j];
    if (b[j]!=0){
      int ind = 2*num;
      while (ind<=n){
        b[ind]=0;
        ind+=num;
      }
    }
    j++;
  }
}
/**
## The loop
The Z-order (or (capital)N-order) space-filling curve seems to be most suitable to study $N\times N$ grids where $N$ is a power of 2 (e.g. 2,4,128 etc.). In order to study the behaviour of the prime-number locations in the Z-order indexed grid we perform our analysis on an increasingly larger grid.   
*/
int main(){
  char name[100];
  FILE * fp2 = fopen("connectedregions","w");
  FILE * fp1 = popen ("ppm2gif --delay 200 > g.gif ", "w");
  /** 
  Loop over increasingly larger grids
  */
  for (int maxlevel=1;maxlevel<=i;maxlevel++){  
    init_grid(1<<maxlevel);
    int d[1<<(maxlevel*dimension)];
    getprimes(d,1<<(maxlevel*dimension));
    int m = 1;
    /**
    Mark cells at prime locations along the space-filling curve:
    */
    scalar field[];
    foreach(){
      field[]=d[m++];
    }
    /**
    We can view the result by using this line of output.
    */
    output_ppm (field, fp1,n=pow(2,i),min=0,max=1,map = gray);
    /**
    Here it is, for all the iterations:
    
    ![White indicates locations of the prime numbers](primes/g.gif)
    
    ## Connected regions
    We see that there exist connected regions, these could potentially have interesting properties. Hence, we tag the connected regions with a unique *tag*. 
    */
    int am = tag(field);
    int regsize[am];
    for (m=0;m<am;m++)
      regsize[m]=0;
    /**
    And we store the length (i.e. size) of each region.
    */
    foreach(){
      if (field[]>0)
        regsize[(int)field[]]++;
    }
    sprintf(name,"prime%d.dat",maxlevel);
    FILE * fp = fopen (name, "w");
    for (m=1;m<am;m++)
      fprintf(fp,"%d\n",regsize[m]);
    /**
    For the largest grid (i.e. $512 \times 512$ points), we plot some statistics on the lengths of the connected regions.
   
   ~~~gnuplot
   
reset
set key box 
set border 3
file='prime9.dat'
binwidth=1
bin(x,width)=width*floor(x/width)

set table "data"

plot file using (bin($1,binwidth)):(1.0) smooth freq with boxes


unset table

set boxwidth 1
set logscale y
set yrange [0.1:20000]
set xlabel 'Length of connected region' 
set ylabel 'Frequency' 
plot "data" index 0 using ($1):($2) with boxes title 'Freqency',\
 100000*exp(-2*x) t 'e^{-2x} scaling ' lw 3
~~~

A truly remarkable feature is that the frequency ($f(i)$) of the region's lengths ($i$) appears to scale with the region's length, according to;

$$f(i) \propto e^{-2i}.$$

We also log the number of connected regions for each grid-refinement iteration. 
 */
    fprintf(fp2,"%d\t%d\n",maxlevel,am+1);
  }
  
 /**
 The result of this procedure are plotted below:
 ~~~gnuplot
 reset
 set xr[0.5:9.5]
 set xlabel 'Refinement level'
 set ylabel 'Number of connected regeons'
 set logscale y
 set key box bottom right
 set size square
 plot 'connectedregions' using 1:2 title 'Number of connected regions',\
 0.1*2**(1.9*x) with lines title '1.9-Dimensional scaling' lw 3
 ~~~
 
 Again we obtain a result that I did not expect, but I am not an expert. 
 */
}

/**
## The next step
The next step is to increase the dimensionality of our Z-order-indexing curve and do the same analysis in 3 dimensions. The results are presented [here](primes3D.c).
*/
