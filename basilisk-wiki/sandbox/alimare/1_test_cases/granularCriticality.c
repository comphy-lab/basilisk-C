/**
# Granular flow, sand heap formation

Freely adapted from the [Bak, tang, Wiesenfeld](#bak1988self) paper for a fortran project
by I. Delbende. It is a cellular automaton model.

![Animation of the formation of the sand](granularCriticality/movie_normal.mp4)

Notice that we do not use any grid.
*/


#define MAXLEVEL 7
#define nt_grains 4
#define hmax 50

/**
## Noise function

Random number between 0 and 1.
*/

static inline double mynoise(){
  return ((noise()+1.)/2.);
}

/**
## Output function
For outputing infos
*/

void output(FILE * fp, int stack[1<<MAXLEVEL], int j, int fileNumber){
  fprintf (fp,
     "set term png notransparent truecolor\n"
     "set yrange [0:%d]\n"
     "set xrange [0:%d]\n"
     "set output 'gnuplot/plot-%06d.png'\n",hmax,(int)(1.5*hmax),fileNumber);
     // "stats 'out' nooutput\n", 
  fprintf (fp,
     "set title 't = %06d'\n"
     "plot '-' u 1:2 w l t sprintf('t=%d')\n",
     j,j);
  for (int k = 0; k < ((1<<MAXLEVEL) - 2 ); ++k){
    fprintf (fp, "%d %d\n", k, stack[k]);
  }
  fprintf (fp, "e\n");
  fflush (fp);
}

int main(){

/**
## The CA model

We model the presence of sand by $r_{max} = 2^{\text{MAXLEVEL}}$ "columns"
denoted
*stack*. The cellular automaton for grain displacement is the following:

* if $stack_{i} - stack_{i+1} \geq 2$ then transfer $n_{g}$ grains from $i$ to
$i+1$ where
$$
n_{g} =  1+ \{0.5*(2 + stack_{i} - stack_{i+1} )\xi \}
$$



*/

  int stack[1<<MAXLEVEL];

  for (int i = 0; i < (1<<MAXLEVEL); ++i){
    stack[i] = 0;
  }

  // temporal loop
  int j=0;
  system ("rm -r gnuplot");
  system ("mkdir gnuplot");
  FILE * fp = popen ("gnuplot 2> /dev/null", "w");

  fprintf (fp,
   "unset key\n"
   "set xlabel 'x'\n"
   "set ylabel 'Column height'\n"
   );
  int mycount = 0;
/**
There is a flow of grains on the first column, we add a grain every $n_{t}$
iterations. The calculation is stopped when the first column reaches $h_{max}$
height.
*/

  while(stack[0]< hmax){
    j++; if(j%nt_grains ==0)stack[0]++;
    fprintf(stderr, "##%d\n", stack[0]);

    for (int k = 0; k < ((1<<MAXLEVEL) - 2 ); ++k){
      int diff = stack[k]-stack[k+1];
      if(diff >= 2){ // grain transfer
        int ng = 1+ floor(0.5*(2.+diff)*mynoise());
        stack[k]   -=ng;
        stack[k+1] +=ng;
        mycount++;
        // output(fp,stack,j,mycount);
      }
    }
    mycount++;
    if(j%10==0)output(fp,stack,j,mycount);
  }


  system ("for f in gnuplot/plot-*.png; do"
    " convert $f ppm:- && rm -f $f; done | "
    "ppm2mp4 movie_normal.mp4");
  fprintf (stderr, "\n\nDone\n");

  return 0;
}

/**
~~~bib

@article{bak1988self,
  title={Self-organized criticality},
  author={Bak, Per and Tang, Chao and Wiesenfeld, Kurt},
  journal={Physical review A},
  volume={38},
  number={1},
  pages={364},
  year={1988},
  publisher={APS}
}

~~~
*/