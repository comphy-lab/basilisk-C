/**
![[This 5 Megapixel image](https://upload.wikimedia.org/wikipedia/commons/2/21/Mandel_zoom_00_mandelbrot_set.jpg) (click for full resolution) provided via Wikipedia does not statisfy our thirst for resolution. Image courtesy of Wolfgang Beyer.](https://upload.wikimedia.org/wikipedia/commons/thumb/2/21/Mandel_zoom_00_mandelbrot_set.jpg/220px-Mandel_zoom_00_mandelbrot_set.jpg)

# A High-Resultion Visualization of the Mandelbrotset Using a Quadtree Grid. 

For a given complex number $c$ one can define the sequence of numbers $z_{n}$ indexed with n as,

$$z_{n+1}=z_n^2+c,$$

starting from $z_0=c$. It can be shown that the sequence diverges in
an absolute sense if $z_nz_n^* > 2$. The mandelbrot set can be
visualized using color coding for the amount of iterations it takes to
get $z_n$ such that $z_nz_n^* \geq 2$. Furthermore, there are
sequenses that do not diverge (eg. the sequence corresponding to
$c=-1$). Therefore we choose to stop checking for divergence after
$n_{max}=2000$ iterations.

On this page we have a look at the scaling of the number of grid cells required to make a smooth looking image of the set. 
We will produce an image at the maximum resolution of 16384 x 16384 pixels (=268 Megapixels), or 14 levels of quadtree refinement. A link will be provided to the result of a 16 Gigapixel redering, using 17 levels of refinement.   

Notice that this (uncompressed) PPM image will be approx. 800 Megabytes in size. Outputting directly to a more sensible format causes crashes when using the inbuild converter.

For the most part we use a similar methods as presented on to the [zoom-in-of-the-Mandelbrot-set page](mandelbrot.c). 
*/
#include "grid/quadtree.h"
#include "utils.h"

void costumcmap (double cmap[NCMAP][3]){
  for (int i = 0; i < NCMAP; i++){
    cmap[i][0]=sq(sin((double)i*M_PI/64.));
    cmap[i][1]=sq(sin(((double)i+20.)*M_PI/64.));
    cmap[i][2]=sq(sin(((double)i+40.)*M_PI/64.));           
  }
  cmap[NCMAP-1][0]=0;
  cmap[NCMAP-1][1]=0;
  cmap[NCMAP-1][2]=0;
}

double xz=-0.75;
double yz=0;
scalar cc[],m[],c[];

double jmax=2000.;
int maxlevel = 14;
int minlevel = 3;
double f=1.4;

/**
   We define a function that calculates how many iterations it takes before the aforementioned sequence diverges. If it does not diverge within *jmax* iterations, it will return *jmax*. Also we track the total number of iterations made for each snapshot.
*/

double determine_iterations(double x,double y){
  double j=0.;
  double xx=xz+(x/f); 
  double yy=yz+(y/f);
  double a=xx;
  double b=yy;
  double c;
  double ab=(sq(a)+sq(b));
  while (j<jmax && ab < 4){
    c=(a*a)-(b*b)+xx;// Real part 
    b=(2*a*b)+yy; //Imag part
    a=c;
    j+=1.;
    ab=((a*a)+(b*b));
  }
  return (double) j;
}

int main(){
  double l10 = log(10.);
  X0=-2.;
  Y0=-2.;
  L0=4.;
  FILE * fp2 = fopen("frac.dat","w");
  init_grid(1<<minlevel);
  foreach(){
    cc[]=determine_iterations(x,y);
    c[]=log(cc[]+1)/l10;
    if (cc[]<jmax)
      m[]=nodata;
    else//mask
      m[]=-1.;
  }
  int n;
  /**
  We iteratively increase the maximum level of refinement. 
  */
  for(int l = minlevel; l<=maxlevel; l++){
    fprintf(ferr,"l=%d\n",l);
    boundary({c});
    /**
    We update the grid where the color-coded field $c$ is not sufficiently smooth.
    */
    while(adapt_wavelet({c},(double[]){3.2/256.},l,minlevel).nf)
      boundary({c});
    n=0;
    foreach(){
      n++;
      cc[]=determine_iterations(x,y);
      c[]=log(cc[]+1.)/l10;
    }
    /**
    The used number of grid cells are stored in a file.
    */
    fprintf(fp2,"%d\t%d\n",l,n);
    fflush(fp2);
  }
  /**
  At the maximum resolution, two PPM images are generated. One with the color-coded field, and one displaying the used resolution. 
  */
  scalar lev[];
  foreach(){
    lev[]=level;
    if (cc[]<jmax)
      m[]=nodata;
    else //mask
      m[]=-1.;
  }
  n=(1<<maxlevel);
  printf("outputting PPM\n");
  FILE * fp = fopen ("nja.ppm","w");
  boundary({c,m});
  output_ppm(c, fp, n=n, linear=true, mask=m, map=costumcmap, min = 0, max=3.2);
  fclose(fp);
  FILE * fpl = fopen ("lev.ppm","w");
  output_ppm(lev, fpl, n=n, min = minlevel, max=depth());
}
/**
# Results

~~~gnuplot
set xr [2.5:14.5]
set yr [10 : 30000000]
set logscale y
set xlabel 'Level of Refinement' font "Hershey/Complex_Roman,12"
set ylabel 'Used Grid Cells' font "Times-Roman,12"
set key top left

plot "frac.dat" t 'Data' ,\
        2**(1.5*x) t '1.5 D scaling'
~~~

The fractal dimension for the color-coded field is approx 1.5. Notice that formally the Mandelbrot set itself (the black boundary) has a Haussdorf dimension of 2. I *think* we do not find that value for the scaling because we have set a maximum number of iterations and the refinement is focussed upon the color-coded field. Notice that the quadtree grid for the computations contained 10 Million grid cells instead of 268 Million for a regular cartesian approach. The corresponding ratio is likely to increase further when generating higher resultion images (1.5D scaling instead of 2D). In fact we have extended the methods to 17 levels of refinement and the 1.5D scaling was maintained. A link to the resulting 16 Giga pixel image is presented below. Notice this required the use of a single core for a few hours using approx. 200GB or RAM memory, a rather silly hardware layout.  

For the 17 levels of refinement run:  

~~~bash
  *** MTRACE: max traced memory size: 204907454271 bytes (tracing overhead 0.002%)  
~~~

## Visuals:

The 268 Mega pixel images are quite a lot to digest and the 16 Giga pixels even more. Therefore, one may have a look at the images via the very nice "Easy-Zoom" website.

* For the Color-coded field, [follow think link to Easy-Zoom](https://easyzoom.com/image/109963) (900 MB, PNG formatted image).  
* The Resolution field upto 14 levels of refinement, [follow another link to Easy-Zoom](https://easyzoom.com/image/109845) (2 MB, PNG formatted image).

For the output of the images with more than one Gigapixel a modified version of the `output_PPM` function was used. 

# Is this record breaking?
[No](http://www.fractalforums.com/mandelbrot-and-julia-set/largest-image-of-the-mandelbrot-ever-created/), but it's quite close actually, and it appears the algorithm presented here is quite efficient (this scripts takes minutes rather than hours and scales well when it is run in parallel). 

Altough, in the end, the image in a 2D object, and the file IO is really hampering performance and the benefit of adaptivity is lost a bit.
*/