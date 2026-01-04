/**
Test for the remap function in C in comparison with the remap that was previously used in Fortran. We validate our implementation by using 3 test functions : a polynom of order 2, a piecewise constant function and a sum of gaussians. */

#include "ppr/ppr.h"
#include "layered_perso/remap-util.h"


/**
Test functions */

double polynom_a(double z) {
    return sq(z);
}

double integral_polynom(double za,double zb) {
    int N=50;
    double integral=0.;
    double zi[N];
    for (int n=0;n<N;n++) {
        zi[n]=za+n*(zb-za)/(N-1.);
    }
    for (int n=0;n<N-1;n++) {
        integral+=(polynom_a(zi[n+1])+polynom_a(zi[n]))/2.*(zi[n+1]-zi[n]);
    }
    integral/=(zb-za);
    return integral;
}

double sum_gaussian(double z) {
    return exp(-sq(z+6.))+3./4.*exp(-sq(z+3.))+2./3.*exp(-sq(z))+1./2.*exp(-sq(z-3.))+1./3.*exp(-sq(z-6.));
}

double integral_sum_gaussian(double za,double zb) {
    int N=50;
    double integral=0.;
    double zi[N];
    for (int n=0;n<N;n++) {
        zi[n]=za+n*(zb-za)/(N-1.);
    }
    for (int n=0;n<N-1;n++) {
        integral+=(sum_gaussian(zi[n+1])+sum_gaussian(zi[n]))/2.*(zi[n+1]-zi[n]);
    }
    integral/=(zb-za);
    return integral;
}

double sharp(double z) {
    if ((z>=-7.)&&(z<-3.)) {
        return 4./10.;
    }
    else if ((z>=-3.)&&(z<1.)) {
        return 12./10.;
    }
    else if ((z>=1.)&&(z<4.)) {
        return 8./10.;
    }
    else {
        return exp(-1./2.*sq(z-9.));
    }
}

double integral_sharp(double za,double zb) {
    int N=50;
    double integral=0.;
    double zi[N];
    for (int n=0;n<N;n++) {
        zi[n]=za+n*(zb-za)/(N-1.);
    }
    for (int n=0;n<N-1;n++) {
        integral+=(sharp(zi[n+1])+sharp(zi[n]))/2.*(zi[n+1]-zi[n]);
    }
    integral/=(zb-za);
    return integral;
}



int main() {

    int nposinit=60;
    int nposend=54;
    int ndof=1;
    int nvar=1;

    double zinit[60];
    double zend[54];

    for (int i=0;i<nposinit;i++) {
        zinit[i]=-10.+i*(10.-(-10.))/(nposinit-1.);
    }
    zend[0]=-10.;
    zend[nposend-1]=10.;
    for (int i=1;i<nposend-1;i++) {
        zend[i]=-10.+i*(10-(-10.))/(nposend-1.)+0.1*noise();
    }
  
  /**
  Parameter for the Fortran code */
  int edge_meth=p3e_method;
  int cell_meth=ppm_method;
  int cell_lim=null_limit;//mono_limit;

  int bc=0.;
  
  /**
  Parameter for the C code */
    double f_b=100.;
    double lambda_b=0.;
    double f_t=100.;
    double lambda_t=0.;

  int Nremap=100;  //number of remap
  
  
  /**
  Polynom test */
    double init_pc[59]={0.};  //number of point=npos-1 because those are mean values on each interval
  double init_pf[59]={0.}; 
    double init2_p[59]={0.};
    double end_pc[53]={0.};
  double end_pf[53]={0.};

    for (int i=0;i<nposinit-1;i++) {
        init_pc[i]=integral_polynom(zinit[i],zinit[i+1]);
      init_pf[i]=integral_polynom(zinit[i],zinit[i+1]);
        init2_p[i]=init_pc[i];
    }


    for (int n=0;n<Nremap;n++) {
        my_remap_perso_perso_C (&nposinit, &nposend, &ndof, &nvar, zinit, zend, init_pc, end_pc,
		    &edge_meth, &cell_meth, &cell_lim, &f_b,&lambda_b,&f_t,&lambda_t);
        my_remap_perso_perso_C (&nposend, &nposinit, &ndof, &nvar, zend, zinit, end_pc, init_pc,
		    &edge_meth, &cell_meth, &cell_lim, &f_b,&lambda_b,&f_t,&lambda_t);
        my_remap (&nposinit, &nposend, &ndof, &nvar, zinit, zend, init_pf, end_pf,
		    &edge_meth, &cell_meth, &cell_lim);
        my_remap (&nposend, &nposinit, &ndof, &nvar, zend, zinit, end_pf, init_pf,
		    &edge_meth, &cell_meth, &cell_lim);
    }
                        
  FILE * test_polynom=fopen("test_polynom.txt","w");
    for (int i=0;i<nposinit-1;i++) {
        fprintf(test_polynom,"%d %g %g %g %g\n",i,zinit[i],init2_p[i],init_pc[i],init_pf[i]);
        fflush(test_polynom);
    }
 
  /**
~~~gnuplot Polynom
plot "test_polynom.txt" using 2:3 title "init" with boxes fc 'gray', "test_polynom.txt" using 2:4 title "remap C++" lc 'red', "test_polynom.txt" using 2:5 title "remap Fortran" lc 'black'
~~~
*/
    /**
  Piecewise constant test */
  
  f_b=0.;
    lambda_b=0.;
    f_t=0.;
    lambda_t=HUGE;
  
  cell_lim=mono_limit;
  
    double init_sc[59]={0.};  //number of point=npos-1 because those are mean values on each interval
  double init_sf[59]={0.}; 
    double init2_s[59]={0.};
    double end_sc[53]={0.};
  double end_sf[53]={0.};

    for (int i=0;i<nposinit-1;i++) {
        init_sc[i]=integral_sharp(zinit[i],zinit[i+1]);
      init_sf[i]=integral_sharp(zinit[i],zinit[i+1]);
        init2_s[i]=init_sc[i];
    }


    for (int n=0;n<Nremap;n++) {
        my_remap_perso_perso_C (&nposinit, &nposend, &ndof, &nvar, zinit, zend, init_sc, end_sc,
		    &edge_meth, &cell_meth, &cell_lim, &f_b,&lambda_b,&f_t,&lambda_t);
        my_remap_perso_perso_C (&nposend, &nposinit, &ndof, &nvar, zend, zinit, end_sc, init_sc,
		    &edge_meth, &cell_meth, &cell_lim, &f_b,&lambda_b,&f_t,&lambda_t);
        my_remap(&nposinit, &nposend, &ndof, &nvar, zinit, zend, init_sf, end_sf,
		    &edge_meth, &cell_meth, &cell_lim);
        my_remap(&nposend, &nposinit, &ndof, &nvar, zend, zinit, end_sf, init_sf,
		    &edge_meth, &cell_meth, &cell_lim);
    }
  
FILE * test_sharp=fopen("test_sharp.txt","w");
    for (int i=0;i<nposinit-1;i++) {
        fprintf(test_sharp,"%d %g %g %g %g\n",i,zinit[i],init2_s[i],init_sc[i],init_sf[i]);
        fflush(test_sharp);
    }
 
  /**
~~~gnuplot Piecewise constant
plot "test_sharp.txt" using 2:3 title "init" with boxes fc 'gray', "test_sharp.txt" using 2:4 title "remap C++" lc 'red', "test_sharp.txt" using 2:5 title "remap Fortran" lc 'black'
~~~
*/
      /**
  Sum of gaussian test */
  
  f_b=0.;
    lambda_b=0.;
    f_t=0.;
    lambda_t=HUGE;
  
  cell_lim=mono_limit;
  
    double init_gc[59]={0.};  //number of point=npos-1 because those are mean values on each interval
  double init_gf[59]={0.}; 
    double init2_g[59]={0.};
    double end_gc[53]={0.};
  double end_gf[53]={0.};

    for (int i=0;i<nposinit-1;i++) {
        init_gc[i]=integral_sum_gaussian(zinit[i],zinit[i+1]);
      init_gf[i]=integral_sum_gaussian(zinit[i],zinit[i+1]);
        init2_g[i]=init_gc[i];
    }


    for (int n=0;n<Nremap;n++) {
        my_remap_perso_perso_C (&nposinit, &nposend, &ndof, &nvar, zinit, zend, init_gc, end_gc,
		    &edge_meth, &cell_meth, &cell_lim, &f_b,&lambda_b,&f_t,&lambda_t);
        my_remap_perso_perso_C (&nposend, &nposinit, &ndof, &nvar, zend, zinit, end_gc, init_gc,
		    &edge_meth, &cell_meth, &cell_lim, &f_b,&lambda_b,&f_t,&lambda_t);
        my_remap (&nposinit, &nposend, &ndof, &nvar, zinit, zend, init_gf, end_gf,
		    &edge_meth, &cell_meth, &cell_lim);
        my_remap (&nposend, &nposinit, &ndof, &nvar, zend, zinit, end_gf, init_gf,
		    &edge_meth, &cell_meth, &cell_lim);
    }
  
FILE * test_gaussian=fopen("test_gaussian.txt","w");
    for (int i=0;i<nposinit-1;i++) {
        fprintf(test_gaussian,"%d %g %g %g %g\n",i,zinit[i],init2_g[i],init_gc[i],init_gf[i]);
        fflush(test_gaussian);
    }
 /** 
~~~gnuplot Sum of gaussian
plot "test_gaussian.txt" using 2:3 title "init" with boxes fc 'gray', "test_gaussian.txt" using 2:4 title "remap C++" lc 'red', "test_gaussian.txt" using 2:5 title "remap Fortran" lc 'black'
~~~
  */
 

    
}