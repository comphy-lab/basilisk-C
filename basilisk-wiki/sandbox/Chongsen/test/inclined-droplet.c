/**
# Droplet spreading on an inclined plane using EBM_VOF.

![Sessile droplet spreading on an inclined plane](inclined-droplet/movie.mp4)
*/

#include "Chongsen/src/EBM_VOF/myembed.h"
#include "navier-stokes/centered.h"
#include "Chongsen/src/EBM_VOF/embed_contact.h"
#include "Chongsen/src/EBM_VOF/embed_two-phase.h"
#include "Chongsen/src/EBM_VOF/embed_tension.h"

vector tmp_h[], o_interface[], ncc[], hnew1[];

double thetac =  70;
double vinitial = 0;
#define MAXLEVEL 6

#define r0        0.2
#define l0        10*r0
#define ed        l0/pow(2,MAXLEVEL)
#define oxy       ed/4.

#define ooox      (oxy)
#define oooy      (oxy)

#define rho01    1.
#define rho02    1.
#define mu01     75./100000.
#define mu02     75./100000.
#define sigma0   0.1

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

int main()
{
  p.nodump = 0;
  DT = .1;
  size(l0);
  origin(-l0/2,-l0/2);
  rho1 = rho01;
  rho2 = rho02;
   mu1 = mu01;
   mu2 = mu02;
  f.sigma = sigma0;
  tmp_c.height = tmp_h;
  tmp_c.hnew1 = hnew1;
  tmp_c.oxyi = o_interface;
  tmp_c.nc = ncc;
  N = 1 << MAXLEVEL;
  run();
}

event init (t = 0)
{
  //embedded_solid
  solid(cs,fs, x + 2.*oxy + y);
  
  //droplet
  fraction(f, - (sq(y + oooy) + sq(x + ooox) - sq(r0)));

  double v=0;
  foreach(){
    f[]*=cs[];
    
    //set contact angle
    contact_angle[] = thetac;
    
    if (dv()>0)
      v += dv()/cm[]*f[];
  }
  vinitial = v;
}

#define tend 20.

event radiou1 (i+=10; t<=tend){

  static FILE * fpshape = fopen("shape.dat", "a");
  fprintf(fpshape,"%.4f ",t);
  
  //v
  double v = 0;
  foreach()
    if (dv()>0)
      v += sq(Delta)*f[];
  
  fprintf(fpshape,"%.16f %.16f %.16f ",vinitial,v,fabs(v-vinitial)/vinitial);
  
  //r
  //Calculate the droplet radius through contact points.
  //The data of r is in "../Chongsen/test/r.dat"
  vector ms[];
  scalar alphacs[];
  reconstruction_cs(cs, fs, ms, alphacs);
  double xx1=0,yy1=0,xx2=0,yy2=0;
  double dr = 100, dl =100;
  foreach(){
    if ((y < -oooy) && (x > -ooox) && (mark[]==4 || mark[] ==5)) {
      coord n, ns;
      coord p_mof[2]={{10,10},{10,10}};
      coord pp[5]={{10,10},{10,10},{10,10},{10,10},{10,10}};
      for (scalar tmp_c in tmp_interfaces) 
        n = interface_normal (point, tmp_c);
      ns.x = ms.x[];
      ns.y = ms.y[];
      polygon_alpha(f[], n, ns, alphacs[], p_mof, pp);
      for(int i=0;i<=1;i++){
        double d = fabs(x + p_mof[i].x*Delta + y + p_mof[i].y*Delta + 2.*oxy)/sqrt(2);
        if (d < dr){
            xx1 = x + p_mof[i].x*Delta;
            yy1 = y + p_mof[i].y*Delta;
            dr = d;
        }
      }
    }
    else if ((y > -oooy) && (x<-ooox) && (mark[]==4 || mark[] ==5)) {
      coord n,ns;
      coord p_mof[2]={{10,10},{10,10}};
      coord pp[5]={{10,10},{10,10},{10,10},{10,10},{10,10}};
      for (scalar tmp_c in tmp_interfaces) 
        n = interface_normal (point, tmp_c);
      ns.x = ms.x[];
      ns.y = ms.y[];
      polygon_alpha(f[], n, ns, alphacs[], p_mof, pp);
      for(int i=0;i<=1;i++){
        double d = fabs(x + p_mof[i].x*Delta + y + p_mof[i].y*Delta + 2.*oxy)/sqrt(2);
        if (d < dl){
            xx2 = x + p_mof[i].x*Delta;
            yy2 = y + p_mof[i].y*Delta;
            dl = d;
        }
      }
    }
  }

double theor = sqrt(v/(thetac*pi/180-cos(thetac*pi/180)*sin(thetac*pi/180)))*sin(thetac*pi/180);

  if (dr<10 && dl<10)
    fprintf(fpshape,"%.16f %.16f %.16f ",theor,sqrt(pow(xx1-xx2,2)+pow(yy1-yy2,2))/2,(sqrt(pow(xx1-xx2,2)+pow(yy1-yy2,2))/2-theor)/theor);
  else
    fprintf(fpshape,"\n");

  //h
  double hhh = 0;
  foreach(){
      if (mark[]==6) {
      coord n;
      for (scalar tmp_c in tmp_interfaces) 
        n = interface_normal (point, tmp_c);
      double alpha = plane_alpha (f[], n);
      coord segment[2];
      if (facets (n, alpha, segment) == 2){
        for (int i = 0 ;i <=1;i++)
          hhh = max(hhh, fabs(x + segment[i].x*Delta + y + segment[i].y*Delta + 2.*oxy)/sqrt(2.));
      }
    }
  }
  double theoh = theor*(1/sin(thetac*pi/180)-1/tan(thetac*pi/180));
  fprintf (fpshape,"%.16f %.16f %.16f\n",theoh,hhh,(hhh-theoh)/theoh);
  fflush(fpshape);
}


#include "view.h"
event movie(t += tend/300){
  view (quat = {0.000, 0.000, 0.000, 1.000},
        fov = 30, near = 0.01, far = 1000,
        tx = -0.00, ty = -0.00, tz = -1.0,
        width = 500, height = 500);
  squares (color = "f");
  draw_vof ("cs", fc = {0.6,0,0}, lc = {1,0,0});
  draw_vof ("tmp_c");
  save ("movie.mp4");
}

#if TREE
event adapt (i++) {
    scalar sf1[];
    foreach() {
    sf1[] = (8. * tmp_c[] +
        4. * (tmp_c[-1] + tmp_c[1] + tmp_c[0, 1] + tmp_c[0, -1] + tmp_c[0, 0, 1] + tmp_c[0, 0, -1]) +
        2. * (tmp_c[-1, 1] + tmp_c[-1, 0, 1] + tmp_c[-1, 0, -1] + tmp_c[-1, -1] +
            tmp_c[0, 1, 1] + tmp_c[0, 1, -1] + tmp_c[0, -1, 1] + tmp_c[0, -1, -1] +
            tmp_c[1, 1] + tmp_c[1, 0, 1] + tmp_c[1, -1] + tmp_c[1, 0, -1]) +
        tmp_c[1, -1, 1] + tmp_c[-1, 1, 1] + tmp_c[-1, 1, -1] + tmp_c[1, 1, 1] +
        tmp_c[1, 1, -1] + tmp_c[-1, -1, -1] + tmp_c[1, -1, -1] + tmp_c[-1, -1, 1]) / 64.;
    sf1[]+=cs[];
    }
  adapt_wavelet ({sf1}, (double[]){1e-3}, minlevel = 3, maxlevel = MAXLEVEL);
}
#endif

/**

## References

~~~bib
~~~
*/