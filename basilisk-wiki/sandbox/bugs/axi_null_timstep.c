
/**

The boundary() are not imposed on cm[] at the end of the initialisation in axi.h.
As a result, if someone tries to compute a timestep before applying boundary((scalar *){cm}), then the obtained timestep is null.

According to S. Popinet in [https://groups.google.com/g/basilisk-fr/c/ZrykBQI6K5I/m/aH6Qzb5GAAAJ](https://groups.google.com/g/basilisk-fr/c/ZrykBQI6K5I/m/aH6Qzb5GAAAJ) the problem may be linked with the "non-symmetrical way" cm[] is called in timestep.h.
However, we see below that make it symmetric doesn't solves the issue.

After applying the BC on cm, the timestep is computed as expected.

*/

#include "axi.h"
#include "run.h"

#include "timestep.h"


double symmetric_timestep (const face vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt = Delta/fabs(u.x[]);
#if EMBED
      assert (fm.x[]);
      dt *= fm.x[];
#else
      dt *= (cm[]+cm[-1])/2.;   // "symmetric call" of cm[]
#endif
      if (dt < dtmax) dtmax = dt;
    }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}


int main()
{
  N = 2;
  run();
}



event plot_cm(i=0){

  face vector uf[];
  foreach_face()
    uf.x[]=fm.x[];   // 1.*fm

  double my_dt = timestep(uf,1.);  // the boundary on cm[] is not done before the loop, thus the timestep is null !
  fprintf(stderr,"%.10g\n", my_dt);


  my_dt = symmetric_timestep(uf,1.); // doesn't work either with a "symmetric epression" for cm in the timstep computation: the timestep is still null
  fprintf(stderr,"%.10g\n", my_dt);
  
  boundary((scalar *){cm}); // After applying BC on cm[], the timestep is no longer null, which solves the issue.
  my_dt = timestep(uf,1.);
  fprintf(stderr,"%.10g\n", my_dt);

}




