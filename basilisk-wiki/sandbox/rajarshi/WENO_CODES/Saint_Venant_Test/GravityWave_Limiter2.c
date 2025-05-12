/**
#GRAVITY WAVE TEST CASE - LIMITER2 OPTION
*/

#include "grid/multigrid1D.h"
#include "Reconstruction.h"
#define pi 3.14159265359

double amp;
scalar KE[], PE[], TE[];

int main()
{ 
  system ("rm -f Time-Period.dat");
  system ("rm -f Dissipation_Coefficients.dat");
  L0 = 600.;
  X0 = -L0/2.;
  G = 10.;
  periodic(left);
  for(N=32; N<=128; N *=2){
     system ("rm -f KE_Location.dat");
     system ("rm -f PE_Location.dat");
     run();
   }
}


event init (i = 0)
{
  foreach()
    h[] = 1.5 + 0.001*sin(pi*x/75.); 
}


event logfile (i++) {
  foreach(){
    KE[] = 0.5*h[]*u.x[]*u.x[];
    PE[] = G*h[]*h[]/2.;
    TE[] = PE[] + KE[];
   }
  static FILE * fp1 = fopen("KE_Location.dat","w");
  static FILE * fp2 = fopen("PE_Location.dat","w");
  static FILE * fp3 = fopen("Heights.dat","w");
  foreach(){
     if(x == (L0/16. +Delta/2.)){
        fprintf(fp1, "%.10g %.13g \n", t, KE[]);
        fprintf(fp2, "%.10g %.13g \n", t, PE[]);
        fprintf(fp3, "%.10g %.13g \n", t, h[]);
      }
   }

  if(N==64){
    static FILE * fp4 = fopen("KE_Location-Display.dat","w");
    static FILE * fp5 = fopen("PE_Location-Display.dat","w");
    static FILE * fp6 = fopen("TE_Location-Display.dat","w");
    foreach(){
     if(x == (L0/16. +Delta/2.)){
        fprintf(fp4, "%.10g %.13g \n", t, KE[]);
        fprintf(fp5, "%.10g %.13g \n", t, PE[]);
        fprintf(fp6, "%.10g %.13g \n", t, TE[]);
      }
    }
  }
}


event output (t = 110) {

 FILE *fp2 = fopen("Numerical_Solution.dat","w");
 foreach()
    fprintf (fp2, "%g %g %g\n", x, h[], u.x[]);
 fclose(fp2);

 double tn=0, KEn=0;
 double tmid = 0.;
 double KE_Mid = -1.;
 double KE_Lower = -2.;
 char str1[14], str2[14];
 
 FILE *fp3 = fopen("KE_Location.dat","r");
 FILE *fp4 = fopen("Dissipation.dat","w");
 
 while (fscanf(fp3, "%s %s\n", &str1, &str2) != EOF){ 
    tn = atof(str1);
    KEn = atof(str2);
    if( KEn <= KE_Mid && KE_Mid > KE_Lower)
       fprintf(fp4,"%g %g \n",tmid,KE_Mid);
     KE_Lower = KE_Mid;
     KE_Mid = KEn;
     tmid = tn;
   }
 fclose(fp3);
 fclose(fp4);

 double A,B,Aavg=0.,Bavg=0.;
 int ctr=0;
 fp3 = fopen("Dissipation.dat","r");
 fp4 = fopen("Dissipation_Coefficients.dat","a");
 fscanf(fp3, "%s %s\n", &str1, &str2);
 tmid = atof(str1);
 KE_Mid = atof(str2);
 while (fscanf(fp3, "%s %s\n", &str1, &str2) != EOF){
   tn = atof(str1);
   KEn = atof(str2);
   ctr++;
   B = log(KE_Mid/KEn)/(tn-tmid);
   A = KE_Mid/exp(-1.*B*tmid);
   Aavg += A;
   Bavg += B;  
   tmid = tn;
   KE_Mid = KEn;
  }
 Aavg /= ctr;
 Bavg /= ctr;
 fprintf(fp4,"%g %g %g \n",L0/N,Aavg,Bavg); 
 fclose(fp3);
 fclose(fp4);

 fp3 = fopen("Time-Period.dat","a");
 fp4 = fopen("Heights.dat","r");
 tn=0, tmid = 0;
 double PEn = 0;
 double PE_Mid = 100.;
 double PE_Lower = 200.;
 ctr=0;
 while (fscanf(fp4, "%s %s\n", &str1, &str2) != EOF && ctr==0){ 
    tn = atof(str1);
    PEn = atof(str2);
    if(tn>=2){
    if( PEn <= PE_Mid && PE_Mid > PE_Lower){
       printf("\n\n %.10lf %.12lf",tn,PEn);
       fprintf(fp3,"%g %g %g %g \n",L0/N,tmid,5.*sqrt(L0/G),fabs(tmid - 5.*sqrt(L0/G)));
       ctr=1;
     }
     }
     PE_Lower = PE_Mid;
     PE_Mid = PEn;
     tmid = tn;
   }
 fclose(fp3);
 fclose(fp4);
}

/**
## Results

~~~gnuplot Velocity
set output 'Velocity.png'
set xlabel 'X'
set ylabel 'U(X)'
set grid
plot 'Numerical_Solution.dat' u 1:3 w p t 'Numerical Velocity', 'Analytical_Solution.dat' u 1:3 w p t 'Analytical Velocity'
~~~

~~~gnuplot Kinetic Energy
set output 'Kinetic_Energy-LIMITER2.png'
set xlabel 'time'
set ylabel 'KE(time)'
set grid
plot 'KE_Location-Display.dat' u 1:2 w p t 'KINETIC ENERGY LIMITER2', 'KE_Location-Display.dat' u 1:(1.91565e-07*exp(-1*0.00394411*($1))) w l t 'Envelope-Function'
~~~

~~~gnuplot Potential Energy
set output 'Potential_Energy-LIMITER2.png'
set xlabel 'time'
set ylabel 'PE(time)'
set grid
plot 'PE_Location-Display.dat' u 1:2 w p t 'POTENTIAL ENERGY LIMITER2'
~~~

~~~gnuplot Total Energy
set output 'Total_Energy-LIMITER2.png'
set xlabel 'time'
set ylabel 'TE(time)'
set grid
plot 'TE_Location-Display.dat' u 1:2 w p t 'TOTAL ENERGY LIMITER2'
~~~

~~~gnuplot DAMPING COEFFICIENT (B) -> exp(-B*time)
set output 'DampingCoefficient_Convergence.png'
set xlabel 'Resoulution'
set ylabel 'B(Delta)'
set grid
plot 'Dissipation_Coefficients.dat' u 1:3 w p t 'Damping-Coefficient-Convergence'
~~~

~~~gnuplot TIME PERIOD ERROR
set output 'TIME-PERIOD-ERROR-Convergence.png'
set xlabel 'Resoulution'
set ylabel '||T(Delta) - T_analytical(Delta)||'
set grid
plot 'Time-Period.dat' u 1:4 w p t 'TIME-PERIOD-ERROR-Convergence'
~~~

*/
