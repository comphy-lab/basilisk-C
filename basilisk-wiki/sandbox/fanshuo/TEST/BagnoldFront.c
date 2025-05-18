
#define ML 1
#define HYDRO 1
#define MUI 1
#define BAGNOLDDRY 1


#include "grid/multigrid1D.h"
#if !ML
# include "saint-venant.h"
#else // ML
# include "./hydroMT.h"
# define phi q
# if !HYDRO
#   include "layered/nh.h"
# endif
# include "layered/remap.h"
# include "layered/perfs.h"
#endif // ML


const double NU = 0.1, T0 = 240;
const double HR = 1 [1];
const double Ltas = 20 [1];


double Ub(Point point)
{
  double Ia = 0.; 
  double zc = 0.;
  double ans;
  for (int l = - point.l; l < nl - point.l; l++) {
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;

  Ia = (tan(slope)-mu0)*I0/(mu0+deltamu-tan(slope));
 
  ans = 2./3.*Ia/dg*sqrt(G*pow(HR,3)*cos(slope))*(1-pow((1-zc/HR),1.5));
  //ans = sqrt(G*pow(HR,3))*0.1/HR/dg;
  return ans;

}


double Ub2(double zc)
{
  double Ia = 0.; 
  double ans;
  
  Ia = (tan(slope)-mu0)*I0/(mu0+deltamu-tan(slope));
  //Ia = (tan(slope)-mu0)*I0/deltamu;
  ans = 2./3.*Ia/dg*sqrt(G*pow(HR,3)*cos(slope))*(1-pow((1-zc/HR),1.5));

  return ans;

}

double Sb(double zc)
{
  double Ia = 0.; 
  double ans;
  
  Ia = (tan(slope)-mu0)*I0/(mu0+deltamu-tan(slope));
 
  ans = Ia/dg*sqrt(G*pow(HR,1)*cos(slope))*pow((1-zc/HR),0.5);

  return ans;


}


/*- - - - - - - - - - - - - - - - - - - - - - - - - */
FILE *file1;
FILE *file2;
FILE *file3;
int main()
{
 
  // slope = 21 degree
  slope = 0.383864[0];  
  L0 = 400.[1];
  G = 1. [1,-2];
  N =256; 
  
  nu = NU;
  nl =64; 
#if ML
#if NOMETRIC
  max_slope = 0.;
#endif
#if !HYDRO  
  NITERMIN = 2;
#endif
#endif

  I0 = 0.4;
  mu0 = 0.35;
  deltamu = 0.21;
  dg = 0.04;
  rho = 1.0;

  file1 = fopen("profil.dat", "w");
  fclose(file1);
  file2 = fopen("profilfront.dat", "w");
  fclose(file2);
  file3 = fopen("profilfactor.dat", "w");
  fclose(file3);

  DT = 1;
  run();
}



/**
We initialise the topography and the initial thickness of each layer *h*. */

event init (i = 0)
{
  foreach() {
    zb[] = 0.;
#if !ML
    h[] = HR - zb[];
#else
    foreach_layer()
     // h[] = (HR)/nl;
		 h[] = dry;
#endif
  }

 // redefine the h distribution, considering the variation in x and in y

 foreach(){
	 double zz = 0.;
	 double hh = HR/nl;
	 
  zz = -HR/Ltas*x+HR;
   foreach_layer(){
		 if(x<Ltas){
			 if(point.l<((int)(zz/hh))){
         h[] = hh;
			 }
			 else if(point.l==((int)(zz/hh))){
         h[] = zz-(int)(zz/hh)*hh;
			 }
			 else{
         h[] = dry;
			 }
		 }
		// else{
    // eta[] = HR*(1-x/Ltas);
		//	 h[] = eta[]/nl;
		// }
		} 
 }
 


  //periodic(right);
	u.n[left] = neumann(0);
  u.n[right] = neumann(0);
  h[left] = neumann(0);
  h[right] = neumann(0);

//Mettre une vitesse initiale dans le domaine 
 foreach (serial) {
	 foreach_layer(){
          u.x[] = (x<Ltas) ? Ub(point):0;
					//u.x[] = Ub(point);
	}
    }

 boundary ({u.x,h});

 }

event acc(++i){
  foreach (serial) {
	 foreach_layer(){
          u.x[] = u.x[] + G*sin(slope)*dt;
     }
  }
}


#if HYDRO
scalar w = {-1};

event update_eta (i++)
{
  if (w.i < 0)
    w = new scalar[nl];
  vertical_velocity (w, hu, hf);

  /**
  The layer interface values are averaged at the center of each
  layer. */
  
  foreach() {
    double wm = 0.;
    foreach_layer() {
      double w1 = w[];
      w[] = (w1 + wm)/2.;
      wm = w1;
    }
  }
}
#endif // HYDRO

event outputfront (t +=5;t<=T0)
{
  file2 = fopen ("profilfront.dat", "a");
  foreach(serial){
    fprintf(file2,"%g %g \n", x, eta[]);
   }
    
  fprintf (file2,"\n\n");   

  fclose(file2);
  
}


event outputprofilfactor (t +=5;t<=T0)
{ 
	double alpha  = 0.;
	double sum = 0.;
	double carresum = 0.;
  file3 = fopen ("profilfactor.dat", "a");
  foreach(serial){
    		if(eta[]>0.001){
    foreach_layer(){
		    carresum = carresum + h[]*u.x[]*u.x[];
		sum = sum + h[]*u.x[];
		}
    alpha = eta[]*carresum/(sum*sum);
    sum = 0.;
		carresum = 0.;
		}
    fprintf(file3,"%g %g\n",x,alpha);
		   }
    
  fprintf (file3,"\n \n");   

  fclose(file3);
  
}

/**
##Front 

~~~gnuplot (Nx = 256, nl = 50)
 set xlabel "x"
 set ylabel "y/h0" 
 set key font ",15"
 p "profilfront.dat" u 1:($2/0.956) w l t"t=0,5,10,15,20,..."
~~~
*/
/**
# A theoretical solution for the front shape
If we compare the front with a theoretical prediction (G. Saingier et al 2016) 
$$X =H-\frac{\frac{2}{3}\log(1-\sqrt{H})-\frac{1}{3}\log(H+\sqrt{H}+1)+\frac{2}{\sqrt{3}}\arctan{(\frac{2\sqrt{H}+1}{\sqrt{3}})-(\alpha-1)Fr^2\ln(\frac{H^d}{(1-H^{3/2})^{2/3}})}}{d-1}$$  

With Froude number $Fr=\frac{u_0}{\sqrt{gh_0\cos{\theta}}}$, $u_0$ the mean velocity in the whole flow packing; $H=h/h_0$; $X=x(\tan{\theta}-\mu_0)/h_0$; $d=(\tan{\theta}-\mu_0)/\Delta \mu$.

The comparaison gives : 
~~~gnuplot front propagates as a wave (Nx= = 256, nl = 64)
 set xlabel "X-Xf"
 set ylabel "H"
 p "profilfront.dat" u ($1*0.0339-9.75):($2/0.96) i 35 w l t"num solution","../anafront_theta21.dat" u ($1-0.9):($2) w l t"ana solution"
~~~
*/
/**
# shape factor $\alpha$
the shape factor is defined from flow velocity profile
$$
\alpha = \frac{\frac{1}{h}\int_0^hu^2(z)dz}{(\frac{1}{h}\int_0^hu(z)dz)^2}
$$
The variation of $\alpha$ along $x$ can be extracted from simulations.
~~~gnuplot 
 set xlabel "x-xf" font ",15"
 set key font ",15"
 set ylabel "y/h_0" font ",15"
 set y2label "alpha" font ",15"
 set y2range [1:2]
 set yrange [0:1.1]

 set y2tics nomirror 0.1
 set ytics nomirror
 p "profilfactor.dat" u ($1-280):2 i 35 w l axes x1y2 t"profil factor alpha", "profilfront.dat" u ($1-280):($2/0.956) i 35 w l axes x1y1 t"num sol front"
~~~
We notice that the shape factor goes from 1.25 (Bagnold flow) to a unit value (plug flow) as approching the front of the avalanche.
*/
