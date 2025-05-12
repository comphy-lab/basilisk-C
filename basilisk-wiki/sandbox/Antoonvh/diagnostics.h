/**
# Diagnostics header file for the lumped parameter dirunalcycle case

This file uses many of the difinitions from the case set-up for *the*
diurnal cycle case of Van Hooft et al (2019). Therefore this header
file should be included after all the relevant declarations. For some
reason "view.h" needs to be included before those declarations in the
set-up file.
*/
//#include "view.h"    //<- Should already be included in the .c file
#include "profile5b.h" //profile5c.h is better
#include "profflux.h"
#include "domainquan.h"
#include "slicer.h"
#include "lambda2.h"

/**
   A movie is rendered with a frame each half minute. Furthermore, the
   movie lasts approx 2 min and shows $\lambda_2$-isosurfaces that are
   coloured with the local buoyancy. Also two slices of the buoyancy
   gradient and the grid structure at the `front`$=$`back` boundary
   are displayed. Finally, the surface buoyancy is visualized to
   reveal the near surface structures if it is not obscured by the
   $\Lambda_2$ surfaces.
 */
# if (dimension == 3)
float bgc[3];
void bgcsetter (double t){
  double day[3] = {126./256., 192./256., 238./256.};
  double sun[3] = {253./256., 94./256., 83./256.};
  double night[3] = {19./256., 24./256., 98./256.};
  double m = max(sin( 2.*M_PI*t/(24*60)), -0.3);
  if (m > 0.){
    m = sqrt(m);
    for (int l = 0; l<3; l++)
      bgc[l] = m*day[l] + (1-m)*sun[l];
  }else{
    m /= -0.3;
    for (int l = 0; l<3; l++)
      bgc[l] = m*night[l] + (1-m)*sun[l];
  }
}

event bviewer3D(t+=0.5;t<=T){
  double bmax = zi*sq(Nb)*1.0;
  double bmin = B1/Lambda;
  scalar gb[],lambda[];
  lambda2(u,lambda);
  // The value of the $\lambfa_2$ isosurface varies with time. 
  double isoval = min(pow(5./3.,(double)maxlevel-6.)*-0.4*sin(t*2*M_PI/T),-0.05);
  boundary({b});
  foreach(){
    gb[]=0;
    foreach_dimension()
      gb[]+=sq((b[1]-b[-1])/(2*Delta));
    gb[]=log10(sqrt(gb[])+1.);
  }
  boundary({gb});
  double avgb= 0;
  foreach_boundary(bottom reduction(+:avgb))
    avgb += BSURF*sq(Delta);
  avgb/=sq(L0);
  double minbr = avgb - 0.1;
  double maxbr = avgb + 0.1;
  char title[100];
  char des[100];
  sprintf(title,"Time Since Sunrise = %02d:%02d (HH:MM)",	\
	  (int)(t/60), (int)floor(fmod(t,60)));
  sprintf(des, "Lambda2 isosurface: %.2g a.u.", isoval);
  double shift = L0/2.;
  double mshift = -L0/2.;
  coord ny = {0, 1, 0};
  //The camera rotates to show a variable perspective.
  double th = 1.57*sin(t*2.*M_PI/360.); // 1 cycle per 6 "hours"
  double ph = 0.3 + (0.6*cube(max(sin(t*2.*M_PI/360.),0.))); //""
  bgcsetter(t);
  clear();
  view(fov = 17, phi = ph,theta = th,  bg = {bgc[0],bgc[1],bgc[2]},
       width = 1920, height = 1080); 
  //view(fov = 17, phi = ph, theta = th, bg = {bgc[0],bgc[1],bgc[2]},
  //     width = 600, height = 400);
  //cells(n=ny,alpha=0.0011);
  squares("b", n = ny,alpha = 0.001, min = minbr,max = maxbr, linear = true);
  cells(alpha = shift, lw = 2);
  squares("gb", alpha = mshift, min = 0, max=1, linear=true);
  isosurface("lambda", isoval, "b", min = bmin , bmax);
  double tc[3];
  if (t < 13*60){
    tc[0] = 0., tc[1] = 0., tc[2] = 0.;  
  }else
    tc[0] = 1., tc[1] = 1., tc[2] = 1.;
draw_string (title, pos = 3, size = 40 , lc = {tc[0],tc[1],tc[2]}, lw = 5);
draw_string(des, pos = 1, size = 60, lc = {tc[0],tc[1],tc[2]}, lw = 3);
  save ("diurnal.mp4");

 
}
#endif
/**
   Each 48-th part of the simulation run, vertical profiles are
   outputted for some fields that have our interest. Also we dump the
   solution for later use to restore and analyse the solution structure
   in more detail with possibly new knowledge/questions. Furthermore,
   horizontal slices at some heights are outputted.
*/

event profs_and_sclices(t+=T/48.){
  scalar lev[],e[],lb[],uh[],sgsbflx[];
  foreach(){
    sgsbflx[]=Kh.y[]*(b[0,1]-b[])/(Delta);
    uh[]=pow(sq(u.x[])+sq(u.z[]),0.5);
    lev[]=level;
    e[]=0;
    lb[]=0;
    foreach_dimension(){
      e[]+=sq(u.x[]);
      lb[]+=sq((b[1]-b[-1])/(2.*Delta));
    }
    if (lb[]>0)
      lb[]=sqrt(lb[]);
  }
  char fn[100];
  sprintf(fn,"proft=%g",t);
  profile((scalar*){b,e,u,uh,lb,Evis,sgsbflx,lev},fn);
  sprintf(fn,"profvart=%g",t);
  profile_var_flx((scalar*){b,u,uh},fn);
  for (double yp = Y0; yp<=zi*1.1;yp+=zi/5.){
    sprintf(fn,"slicebXZt=%gy=%g",t,yp);
    sliceXZ(fn,b,yp,maxlevel,true);
    sprintf(fn,"sliceuyXZt=%gy=%g",t,yp);
    sliceXZ(fn,u.y,yp,maxlevel,true);
  }
  sprintf(fn,"slicebXYt=%gy=%g",t,0.);
  sliceXY(fn,b,0,maxlevel,false);
  sprintf(fn,"sliceuyXYt=%gy=%g",t,0.);
  sliceXY(fn,u.y,0,maxlevel,false);
}
/**
   Additionally the simulation is dumped 48 times during the simulation run
*/
event dumper(t+=T/48){
  char fn[100];
  sprintf(fn,"dumpt=%g",t);
  dump(fn);
}

/**
   Each time unit (a physical minute) we output some global (i.e. integrated) quantities to
   track the evolution of the solution. The focus is a bit on the energy budgets.
*/

event timeseries(t+=1.){
  static FILE * fpt = fopen("timeseries","w");
  double H = 0, G = 0, Qs = 0,diss = 0, en = 0, bint = 0;
  foreach(reduction(+:H) reduction(+:G) reduction(+:Qs) reduction(+:diss) \
	  reduction(+:en) reduction(+:bint)){
    double cv =1.;
    foreach_dimension()
      cv*=Delta;
    foreach_dimension(){
      en += cv * sq(u.x[]);
      diss+=Evis[]*cv*(sq(u.x[1] - u.x[-1]) +
		       sq(u.x[0,1] - u.x[0,-1]) +
		       sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
    bint += (b[]-sq(Nb)*y)*cv;
    if (y<Delta){
      double m = sq(Delta);
      Qs+=Qn*m;
      G+=Gbflx*m;
      H+=(Qn+Gbflx)*m;
    }
  }
  double endm=e(u);
  double totb=0;
  foreach(reduction(+:totb))
    totb+=b[]*cube(Delta);
  double totbndm=sndom(b);
  if (pid()==0){
    if (t==0)
      fprintf(fpt,"t\ti\tH\tG\tQn\te\ten\tdiss\ttotbndm\ttotb\tbint\tLc\tUc\n");
    fprintf(fpt,"%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", \
	    t,i,H,G,Qs,endm,en,diss,totbndm,totb,bint,Lc,Uc);
    fflush(fpt);
  }
}

/**
   There is also a tall virual observation tower in our domain that
probes the boundary layer at given locations. Instantanious
meassurements are taken each 20 timesteps. The location is *almost* in
the middle of the domain. The format is not very convinient, and as
follows;

Header: t1   i1  
       "y"   s1.name    s2.name ... sn.name             
data:   y1   s_1(y1,t1) ..      ... s_n(y1,t1)      
        y2   s_1(y2,t1) ..      ... s_n(y2,t1)  
        ..   ..         ..      ...  
        yn   s_1(yn,t1) ..      ... s_n(yn,t1)  
        t2   i2   
	"y"   s1.name   s2.name ... sn.name  
        y1   s_1(y1,t2) ..      ... s_n(y1,t2)      
        |    |          |       |   |  
        yn   s_1(yn,t2) ..      ... s_n(yn,t2)  
        t3   i3  
        |    |  
        tn   in  
        |    |  
        yn   s_1(yn,tn) ..      ... s_n(yn,tn)  
   

where `\t` and `\n` are used for tabs and linebreaks, respectively.
*/
event tower(i += 20){
  int No  = 301;
  double maxy = 1.5*zi;
  double miny = Y0;
  double dy = (maxy-miny)/((double)(No - 1));
  scalar * list = {u.x, u.y, u.z, b};
  boundary(list);
  int len = list_len(list);
  double data[len][No];
  static FILE * fpt = fopen("towerdata","w");
  double xp = ((M_PI*sqrt(2.)*3./5.)/(L0*100.)) + X0 + L0/2;
  double zp = ((M_PI*sqrt(2.)*3./5.)/(L0*100.)) + Z0 + L0/2;
  int jj = 0;
  fprintf(fpt, "%g\t%d\ny", t, i);
  for (scalar s in list)
    fprintf(fpt, "\t%s", s.name);
  fprintf(fpt, "\n");
  for (double yp = miny; yp <= maxy; yp += dy){
    Point point = locate(xp, yp, zp);
    int ii = 0;
    for (scalar s in list){
      data[ii][jj] = point.level >= 0 ? interpolate(s, xp, yp, zp): 0.;
      ii++;
    }
    jj++;
  }
#if _MPI
  if (pid() == 0){ //Favorite worker
    MPI_Reduce(MPI_IN_PLACE, &data, No*len, MPI_DOUBLE, MPI_SUM, 0,
		  MPI_COMM_WORLD);
  }else
    MPI_Reduce (&data[0], NULL, No*len, MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
#endif
  jj = 0;
  if (pid() == 0){
    for (double yp = miny; yp <= maxy; yp += dy){
      int ii = 0;
      fprintf(fpt, "%g", yp);
      for (scalar s in list){
	fprintf(fpt, "\t%g", data[ii][jj]);
	ii++;
      }
      fprintf(fpt, "\n");
      jj++;
    }
  }
  fflush(fpt);
}	      

/**
   Each 25 time steps we output some characteristics of the solver. This
   aims to trace the solver's performance.
*/

int gc = 0;
double tc = 0;
event logfile(i+=25){
  fprintf(ferr,"%d\t%g\n",i,t);
  static FILE * fpp = fopen ("perfs", "w");
  if (i == 0)
    fprintf (fpp,
  	     "i t dt mgp.i mgp.nrelax mgu.i mgu.nrelax"
  	     "grid->tn perf.t perf.speed npe\n");
  fprintf (fpp, "%d %g %g %d %d %d %d %ld %g %g %d\n",
  	   i,t, dt, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax,
  	   grid->tn, perf.t, perf.speed, npe());
  fflush (fpp);
#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#if TRACE
  static FILE * fptrace = fopen ("tracing", "w");
  trace_print (fptrace, 1);
#endif
}
