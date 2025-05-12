/**test_restart.c
The restore() function is restting the value of N to 1. 
Make your favourite dump file called test_restart 
and then try and read it in here.
A possible code for making a dump file is:
//********************************************
#include "saint-venant.h"

#define MAXLEVEL 7

int main (){
    init_grid(1<<MAXLEVEL);
    run();
}

event init(i=0){
    foreach(){
      zb[] = 1.-x * x - y * y;
      h[] = max(0.0, 0.5 - zb[]);
    }
}

// Output a dump file
event output(t = 0.){
        dump(file = "test_restart");
	output_ppm(zb, file="bathymetry0.ppm");
	output_ppm(h, file="h0.ppm");
}
//***********************************************

*/
#include "saint-venant.h"

#define MAXLEVEL 7

int main (){
    init_grid( 1 << MAXLEVEL );
    run();
}

event init(i=0){
  printf("Before restore, N = %d\n",N);
  restore(file = "test_restart");
  printf("After restore, N= %d\n",N);
  /* If you uncomment the next line it resets the value of N to what it should be 
  and the .ppm files are output correctly
  */
  //  N=1<<MAXLEVEL;
}
event output(i=1){
  output_ppm(zb, file="bathymetryRS.ppm");
  output_ppm(h,file="hRS.ppm");
}
