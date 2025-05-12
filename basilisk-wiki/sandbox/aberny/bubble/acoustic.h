#include <stdlib.h>

event acoust (i++){
  double xLoc1 = -13.6;
  double xLoc2 = -13.8;
  double xLoc3 = -13.4;

  double yLoc1 = 0.001;
  double yLoc2 = 0.005;
  double yLoc3 = 0.01;

  scalar a[];

  foreach(){
    a[] = p[];
  }

  boundary({a});

  double p1 = interpolate(a, xLoc1, yLoc1);
  double p2 = interpolate(p, xLoc1, yLoc2);
  double p3 = interpolate(p, xLoc1, yLoc3);

  double p4 = interpolate(p, xLoc2, yLoc1);
  double p5 = interpolate(p, xLoc2, yLoc2);
  double p6 = interpolate(p, xLoc2, yLoc3);

  double p7 = interpolate(p, xLoc3, yLoc1);
  double p8 = interpolate(p, xLoc3, yLoc2);
  double p9 = interpolate(p, xLoc3, yLoc3);

  static FILE * fp1 = fopen("pressure1.dat", "w");
  if (i == 0){
    fprintf(fp1, "i  t/tc  x/Rb  y/Rb  p\n");    
  }
  fprintf(fp1, "%d %g %g %g %g\n", i, t, xLoc1, yLoc1, p1);

  static FILE * fp2 = fopen("pressure2.dat", "w");
  if (i == 0){
    fprintf(fp2, "i  t/tc  x/Rb  y/Rb  p\n");    
  }
  fprintf(fp2, "%d %g %g %g %g\n", i, t, xLoc1, yLoc2, p2);

  static FILE * fp3 = fopen("pressure3.dat", "w");
  if (i == 0){
    fprintf(fp3, "i  t/tc  x/Rb  y/Rb  p\n");    
  }
  fprintf(fp3, "%d %g %g %g %g\n", i, t, xLoc1, yLoc3, p3);

  static FILE * fp4 = fopen("pressure4.dat", "w");
  if (i == 0){
    fprintf(fp4, "i  t/tc  x/Rb  y/Rb  p\n");    
  }
  fprintf(fp4, "%d %g %g %g %g\n", i, t, xLoc2, yLoc1, p4);
  
  static FILE * fp5 = fopen("pressure5.dat", "w");
  if (i == 0){
    fprintf(fp5, "i  t/tc  x/Rb  y/Rb  p\n");    
  }
  fprintf(fp5, "%d %g %g %g %g\n", i, t, xLoc2, yLoc2, p5);
  
  static FILE * fp6 = fopen("pressure6.dat", "w");
  if (i == 0){
    fprintf(fp6, "i  t/tc  x/Rb  y/Rb  p\n");    
  }
  fprintf(fp6, "%d %g %g %g %g\n", i, t, xLoc2, yLoc3, p6);
  
  static FILE * fp7 = fopen("pressure7.dat", "w");
  if (i == 0){
    fprintf(fp7, "i  t/tc  x/Rb  y/Rb  p\n");    
  }
  fprintf(fp7, "%d %g %g %g %g\n", i, t, xLoc3, yLoc1, p7);
  
  static FILE * fp8 = fopen("pressure8.dat", "w");
  if (i == 0){
    fprintf(fp8, "i  t/tc  x/Rb  y/Rb  p\n");    
  }
  fprintf(fp8, "%d %g %g %g %g\n", i, t, xLoc3, yLoc2, p8);
  
  static FILE * fp9 = fopen("pressure9.dat", "w");
  if (i == 0){
    fprintf(fp8, "i  t/tc  x/Rb  y/Rb  p\n");    
  }
  fprintf(fp9, "%d %g %g %g %g\n", i, t, xLoc3, yLoc3, p9);

  fflush(fp1);
  fflush(fp2);
  fflush(fp3);
  fflush(fp4);
  fflush(fp5);
  fflush(fp6);
  fflush(fp7);
  fflush(fp8);
  fflush(fp9);




#if Bview
  clear();
  view (fov = 20, quat = {0,0,-0.707107,0.707107}, ty = 0.3, tx = 0, bg = {1,1,1}, width = 1280, height = 720, samples = 4);
  draw_vof("f");
  squares("a", min = -5, max = 5, linear = true, map = cool_warm);
  mirror({0,1,0},0) {
    draw_vof("f");
    squares("a", min = -5, max = 5, linear = true, map = cool_warm);
  }


  static FILE * fpPressure = popen ("ppm2mp4 pressure.mp4", "w");
  save(fp = fpPressure);
#endif
  // dump("pression");

}