/* A constant level of fluid defined by (phi[]=y-2.;) gives the following arithmetic error :
Program received signal SIGFPE, Arithmetic exception. Location of error fractions.h:266 */

event init (t = 0) {
  vertex scalar phi[];
  foreach_vertex(){
  phi[]=y-2.;
  }
  fractions (phi, f);
}
