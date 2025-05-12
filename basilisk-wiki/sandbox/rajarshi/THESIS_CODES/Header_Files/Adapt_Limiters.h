/**
# HIGHER-ORDER ADAPTIVITY WITH LIMITING.

This header-file includes the implementation of the full-limited version of the prolongation function. It is computationally expensive compared to 
the non-limited higher-order adaptive implementation.

## NOTES FOR BASILISK USERS

1. Include this header-file only when your fields have discontinuity. In case of continuous fields, use the header-file Adapt_No_Limiters.h
2. This code is functional for 1D/2D/3D cases.
3. This code has implementations for volume-averaged functions only. Hence it can be used for scalar and vector fields, but for face-vector fields (especially
   face-velocity fields, which need to implement divergence-free condition), further development is needed to build appropriate higher-order functions. 
4. For detailed understanding of the Prolongation function implementation, see Rajarshi-PhD_Thesis - section 2.10.4 on Prolongation operator.
5. The fourth-order Poisson solver uses this file (or its non-limited version by default) to borrow the prolongation function for multigrid implementation.
6. How to implement in Basilisk :-   
   scalar s[];
   s.refine = refine_order5;
   s.prolongation = refine_order5;
   adapt_wavelet ({s},(double []){cmax},maxlevel=MAX);
*/

#include "utils.h"
#define epsilon_P 1e-06
static inline void Prolongation_Quadrature_Order5 (double * T, int marker, double * q) {

  int i, j, k;
  static const double weight[4][3][3] = {

                                         {
                                             {-0.011903542149445, -0.03254374839074,  1.0444472905402   }, 
                                             {-0.068254374839074, 1.0801579169885  , -0.011903542149445 }, 
                                             {0.8753947924713   , 0.19285958236778 , -0.068254374839074 } 
                                         },

                                         {
                                             {-0.011903542149445, -0.03254374839074,  1.0444472905402   }, 
                                             {-0.068254374839074, 1.0801579169885  , -0.011903542149445 }, 
                                             {0.8753947924713   , 0.19285958236778 , -0.068254374839074 } 
                                         }, 
 
                                         {
                                             { 0.11458333333333, -0.47916666666667,  1.3645833333333  }, 
                                             {-0.13541666666667,  1.0208333333333 ,  0.11458333333333 }, 
                                             { 0.61458333333333,  0.52083333333333, -0.13541666666667 }
                                         },
 
                                         { 
                                             { 0.27857020881611, -1.0007895849426 ,  1.7222193761265  }, 
                                             {-0.16507895849426,  0.88650874967815,  0.27857020881611 }, 
                                             { 0.39127187419537,  0.77380708429889, -0.16507895849426 }                                                
                                         } 
                                       };

  static const double gamma[4][3] = {
                                       { 0.18862794688089, 0.80268704050847, 0.0086850126106366  }, 
                                       { 0.4818440346213 , 0.51260956890584, 0.0055463964728521  },  
                                       {  0.22414772727273,  0.60013111888112,  0.17572115384615 }, 
                                       {  0.11823527923684,  0.60955211778235,  0.27221260298081 } 
                                    };


  double IS[3];
  IS[0] = 13.*sq( T[0] - 2.*T[1] + T[2] )/12. + sq(    T[0] - 4.*T[1] + 3.*T[2] )/4.;
  IS[1] = 13.*sq( T[1] - 2.*T[2] + T[3] )/12. + sq(    T[1]       -        T[3] )/4.;
  IS[2] = 13.*sq( T[2] - 2.*T[3] + T[4] )/12. + sq( 3.*T[2] - 4.*T[3] +    T[4] )/4.;
    
  double w[4][3];
  double sum;
  double temp[4];
  if (marker == 1) {
      for ( i = 0; i <= 3; i++ ) {
          for ( j = 0; j <= 2; j++ )  
              w[i][j] = gamma[i][j] / sq(epsilon_P + IS[j] );
          sum = w[i][0] + w[i][1] + w[i][2];
          for ( j = 0; j <= 2; j++ )
              w[i][j] /= sum;
       
          temp[i] = 0.;
          for(j=0; j<=2; j++)
             for(k=0; k<=2; k++)
                temp[i] += w[i][j]*weight[i][j][k]*T[j+k]; 
      }
      q[0] = 4.607061858653*temp[0] - 3.607061858653*temp[1]; 
      q[1] = temp[2];
      q[2] = temp[3];                   
  }
      
  else if (marker == -1) {
      for ( i = 0; i <= 3; i++ ) {
          for ( j = 0; j <= 2; j++ )  
              w[i][j] = gamma[3-i][2-j] / sq(epsilon_P + IS[j] );
          sum = w[i][0] + w[i][1] + w[i][2];
          for ( j = 0; j <= 2; j++ )
              w[i][j] /= sum; 
          temp[i] = 0.;
          for(j=0; j<=2; j++)
             for(k=0; k<=2; k++)
                temp[i] += w[i][j]*weight[3-i][2-j][2-k]*T[j+k];
      }
      q[0] = temp[0];
      q[1] = temp[1];
      q[2] = 4.607061858653*temp[3] - 3.607061858653*temp[2];
  }
}


#if dimension == 1

   static inline double order5_refine (Point point, scalar s) {
   
     double T_hold[5];
     double QuadV[3];
  
     for (int i = -2; i <= 2; i++)
         T_hold[i+2] = coarse(s,i);
  
     Prolongation_Quadrature_Order5 (T_hold, child.x, QuadV);
     return (5.*QuadV[0] + 8.*QuadV[1] + 5.*QuadV[2])/18.;
}

#elif dimension == 2

   static inline double order5_refine (Point point, scalar s) {
   
     double T_hold[5];
     double QuadPtr[3];
     double QuadArr1[5][3];
     double QuadArr2[3];

     for (int j = -2; j <= 2; j++) {
        for (int i = -2; i <= 2; i++)
            T_hold[i+2] = coarse(s,i,j);
        Prolongation_Quadrature_Order5 (T_hold, child.x, QuadPtr);
        for (int i = 0; i <= 2; i++)
            QuadArr1[j+2][i] = QuadPtr[i];
     }

     for (int j = 0; j <= 2; j++) {
         for (int i = 0; i <= 4; i++)
             T_hold[i] = QuadArr1[i][j];
         Prolongation_Quadrature_Order5 (T_hold, child.y, QuadPtr);
         QuadArr2[j] = (5.*QuadPtr[0] + 8.*QuadPtr[1] + 5.*QuadPtr[2])/18.;
      }

      return (5.*QuadArr2[0] + 8.*QuadArr2[1] + 5.*QuadArr2[2])/18.;
}

#elif dimension == 3
  
   static inline double order5_refine (Point point, scalar s) {
     
     double T_hold[5];
     double QuadPtr[3];
     double QuadArr1[5][5][3];
     double QuadArr2[5][3][3];
     double QuadArr3[3][3][3];

     for (int k = -2; k <= 2; k++)
         for (int j = -2; j <= 2; j++) {           
             for (int i = -2; i <= 2; i++)
	         T_hold[i+2] = coarse(s,i,j,k);
             Prolongation_Quadrature_Order5 (T_hold, child.x, QuadPtr);
             for (int i = 0; i <= 2; i++)
          	 QuadArr1[j+2][k+2][i] = QuadPtr[i];
         }

     for (int i = 0; i <= 2; i++)
         for (int j = 0; j <= 4; j++) {          
             for (int k = 0; k <= 4; k++)
	         T_hold[k] = QuadArr1[j][k][i];
             Prolongation_Quadrature_Order5 (T_hold, child.z, QuadPtr);
             for (int k = 0; k <= 2; k++) 
	         QuadArr2[j][k][i] = QuadPtr[k]; 
         }
             
     for (int i = 0; i <= 2; i++)
         for (int k = 0; k <= 2; k++) {
             for (int j = 0;j <= 4; j++)
	         T_hold[j] = QuadArr2[j][k][i];
             Prolongation_Quadrature_Order5 (T_hold, child.y, QuadPtr); 
             for (int j = 0; j <= 2; j++)
	         QuadArr3[j][k][i] = QuadPtr[j];
         }

     double QSum[3] = {5./18., 8./18., 5./18. };
     double QuadSummation = 0.;
     for (int k = 0; k <= 2; k++)
         for (int j = 0; j <= 2; j++)
             for (int i = 0; i <= 2; i++)
	         QuadSummation += QSum[j]*QSum[k]*QSum[i]*QuadArr3[j][k][i];   
  
     return QuadSummation;
}

#endif


static inline void refine_order5 (Point point, scalar s)
{
  foreach_child()
    s[] = order5_refine (point, s);
}

