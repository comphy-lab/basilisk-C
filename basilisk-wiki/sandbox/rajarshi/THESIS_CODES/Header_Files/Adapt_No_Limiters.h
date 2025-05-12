/**
# HIGHER-ORDER ADAPTIVITY WITHOUT LIMITING.

This header-file includes the implementation of the non-limited version of the prolongation function. It is computationally cheaper compared to 
the full-limited higher-order adaptive implementation.

## NOTES FOR BASILISK USERS

1. Include this header-file only when your fields are continuous. In case of fields with discontinuity, use the header-file Adapt_Limiters.h
2. This code is functional for 1D/2D/3D cases.
3. This code has implementations for volume-averaged functions only. Hence it can be used for scalar and vector fields, but for face-vector fields (especially
   face-velocity fields, which need to implement divergence-free condition), further development is needed to build appropriate higher-order functions. 
4. For detailed understanding of the Prolongation function implementation, see Rajarshi-PhD_Thesis - section 2.10.4 on Prolongation operator.
5. The fourth-order Poisson solver uses this file (or its full-limited version) to borrow the prolongation function for multigrid implementation.
6. How to implement in Basilisk :-   
   scalar s[];
   s.refine = refine_order5;
   s.prolongation = refine_order5;
   adapt_wavelet ({s},(double []){cmax},maxlevel=MAX);
*/

static inline void Prolongation_Quadrature_Order5 (double * T, int marker,
						   double * q)
{
  static const double weight[3][5] = {
    {0.0103444235735617, -0.09792213521550909, 1.107094656676454,
     -0.01815143469025667, -0.001365510344249945},
    {0.02568359375, -0.188671875, 1.026497395833333,
     0.1602864583333333, -0.02379557291666666},
    {0.0329368264264383, -0.2189528647844909, 0.8505095099902126,
     0.3804431013569233, -0.04493657298908339}
  };
  if (marker == 1)
    for (int i = 0; i <= 2; i++) {
      q[i] = 0.;
      for (int j = 0; j <= 4; j++)
	q[i] += weight[i][j]*T[j];
    }
  else if (marker == -1)
    for (int i = 0; i <= 2; i++) {
      q[i] = 0.;
      for(int j = 0; j <= 4; j++)
	q[i] += weight[2-i][4-j]*T[j];
    }
}

#if dimension == 1
double order5_refine (Point point, scalar s)
{
  double T_hold[5];
  double QuadV[3];
  
  for (int i = -2; i <= 2; i++)
    T_hold[i+2] = coarse(s,i);
  
  Prolongation_Quadrature_Order5 (T_hold, child.x, QuadV);
  return (5.*QuadV[0] + 8.*QuadV[1] + 5.*QuadV[2])/18.;
}
#elif dimension == 2
static inline double order5_refine (Point point, scalar s)
{
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
static inline double order5_refine (Point point, scalar s)
{
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
