/**
#HEADER FILE - All functions necessary for flux computation (Prolongation Operator included)
*/

#define epsilon 1E-06

double * Prolongation_Quadrature(double *T, int marker){

     int i,j,k;
     double sum,temp;
     static double T_Quad[3];
     double stencil_weights[3][3][3] = 
         {
            {
            { 2./15. - 3.*sqrt(15)/80., -31./60. +    sqrt(15)/8. , 83./60. - 7.*sqrt(15)/80.},
            {-7./60. +    sqrt(15)/80.,  59./60. +    sqrt(15)/40.,  2./15. - 3.*sqrt(15)/80.},
            {19./30. +    sqrt(15)/16.,  29./60. - 3.*sqrt(15)/40., -7./60. +    sqrt(15)/80.}
            },
            {
            { 11./96.,-23./48.,131./96.},
            {-13./96., 49./48., 11./96.},
            { 59./96., 25./48., -13./96.}
            },
            {
            { 2./15. + 3.*sqrt(15)/80., -31./60. -    sqrt(15)/8. , 83./60. + 7.*sqrt(15)/80.},
            {-7./60. -    sqrt(15)/80.,  59./60. -    sqrt(15)/40.,  2./15. + 3.*sqrt(15)/80.},
            {19./30. -    sqrt(15)/16.,  29./60. + 3.*sqrt(15)/40., -7./60. -    sqrt(15)/80.}
            },
        };

     double gamma[3][3],Weights[3][3],Beta[3];
     gamma[0][0] = (277./12800. - 7.*sqrt(15.)/2400.)/(2./15. - 3.*sqrt(15.)/80.);
     gamma[0][2] = (-889./38400. + 9.*sqrt(15.)/1600.)/(-7./60. + (1./80.)*sqrt(15.));
     gamma[0][1] = 1. - gamma[0][0] - gamma[0][2];
     gamma[1][0] = 789./3520.;
     gamma[1][1] = 13731./22880.;
     gamma[1][2] = 731./4160.;
     gamma[2][0] = (277./12800. + 7.*sqrt(15.)/2400.)/(2./15. + 3.*sqrt(15.)/80.);
     gamma[2][2] = (-889./38400. - 9.*sqrt(15.)/1600.)/(-7./60. - (1./80.)*sqrt(15.));
     gamma[2][1] = 1. - gamma[2][0] - gamma[2][2];

     if( marker == -1 ){
     
       for(i=0; i<=1; i++){          
          for(j=0; j<=2; j++){ 
                for(k=0; k<=2; k++){
                   if((i==0) || ((i==1) && (j+k<=2) && j!=2)){             
                      temp = stencil_weights[i][j][k];
                      stencil_weights[i][j][k] = stencil_weights[2-i][2-j][2-k];
                      stencil_weights[2-i][2-j][2-k] = temp;
                   }
                }
             }
          }
        
       for(i=0; i<=1; i++){ 
           for(j=0; j<=2; j++){
               if(i+j<=2 && i!=2){             
                 temp = gamma[i][j];
                 gamma[i][j] = gamma[2-i][2-j];
                 gamma[2-i][2-j] = temp;
                } 
              }
           }
       }

     Beta[0]      = 13.*sq(*T     - 2.*(*(T+1)) + *(T+2))/12. + sq( *T         - 4.*(*(T+1))      + 3.*(*(T+2)))/4.;    
     Beta[1]      = 13.*sq(*(T+1) - 2.*(*(T+2)) + *(T+3))/12. + sq(*(T+1)      -   *(T+3)  )/4.;
     Beta[2]      = 13.*sq(*(T+2) - 2.*(*(T+3)) + *(T+4))/12. + sq(3.*(*(T+2)) - 4.*(*(T+3))      +    *(T+4))/4.;

     for(i=0; i<=2; i++){
       for(j=0; j<=2; j++)
          Weights[i][j] = gamma[i][j]/sq(epsilon + Beta[j]);
       sum = Weights[i][0] + Weights[i][1] + Weights[i][2];
       for(j=0; j<=2; j++)
          Weights[i][j] /= sum;
     }
         
     for(i=0; i<=2; i++){
       T_Quad[i] = 0.;
       for(j=0; j<=2; j++)
         for(k=0; k<=2; k++)
            T_Quad[i] += gamma[i][j]*stencil_weights[i][j][k]*(*(T+j+k));
      }
    
     return (T_Quad);  
}

#if dimension == 1
static inline double WenoRefine1D (Point point, scalar X)
{
     double T_hold[5];
     double * QuadV;

     for(int i = -2; i<=2; i++)
        T_hold[i+2] = coarse(X,i);

     QuadV = Prolongation_Quadrature(T_hold, child.x);
     return(5.*(*QuadV)/18. + 8.*(*(QuadV+1))/18. + 5.*(*(QuadV+2))/18.);
}

#elif dimension == 2
static inline double WenoRefine2D (Point point, scalar X)
{
     double T_hold[5];
     double * QuadPtr;
     double QuadArr1[5][3];
     double QuadArr2[3];
     int i,j;

     for(j = -2; j<=2; j++){
        for(i = -2; i<=2; i++)
            T_hold[i+2] = coarse(X,i,j);
        QuadPtr = Prolongation_Quadrature(T_hold, child.x);
        
        for(i = 0; i<=2; i++)
          QuadArr1[j+2][i] = *(QuadPtr+i);         
     }

     for(j = 0; j<=2; j++){
        for(i = 0; i<=4; i++)
            T_hold[i] = QuadArr1[i][j];
        QuadPtr = Prolongation_Quadrature(T_hold, child.y);
        QuadArr2[j] = 5.*(*(QuadPtr))/18. + 8.*(*(QuadPtr+1))/18. + 5.*(*(QuadPtr+2))/18.;
     }
  
    return(5.*QuadArr2[0]/18. + 8.*QuadArr2[1]/18. + 5.*QuadArr2[2]/18. );
}
#endif

static inline void refine_weno (Point point, scalar s)
{
  #if dimension == 1
    foreach_child()
       s[] = WenoRefine1D (point, s);
  #elif dimension == 2
    foreach_child()
       s[] = WenoRefine2D (point, s);
  #endif
}


double WENO5_Reconstruction_LeftY(Point point, scalar *Xl, double gradL){

  scalar X = Xl[0];
  double T_Stencil[3], Beta[3], Weights[3];
  int i;
  double sum;

  double gammaL[3] = {1./10.,3./5.,3./10.};
  T_Stencil[0]  =   1.*(X[0,-1] - 2.*Delta*gradL)/3. - 7.*X[0,-2]/6. + 11.*X[0,-1]/6.;
  Beta[0]       =  13.*sq(X[0,-1] - 2.*Delta*gradL - 2.*X[0,-2] + X[0,-1])/12. + sq(X[0,-1] - 2.*Delta*gradL - 4.*X[0,-2] + 3.*X[0,-1])/4.;
  T_Stencil[1]  =  -1.*X[0,-2]/6. + 5.*X[0,-1]/6. + 1.*X[]/3.;
  Beta[1]       =  13.*sq(X[0,-2] - 2.*X[0,-1] + X[])/12. + sq(X[0,-2]-X[])/4.;
  T_Stencil[2]  =   1.*X[0,-1]/3. + 5.*X[]/6. - 1.*X[0,1]/6.;
  Beta[2]       =  13.*sq(X[0,-1] - 2.*X[] + X[0,1])/12. + sq(3.*X[0,-1] - 4.*X[] + X[0,1])/4.;
     
  for(i=0; i<=2; i++)
    Weights[i] = gammaL[i]/sq(epsilon + Beta[i]);
  
  sum = Weights[0] + Weights[1] + Weights[2];
  for(i=0; i<=2; i++)
    Weights[i] /= sum;

  sum = 0.;
  for(i=0;i<=2;i++)
   sum += T_Stencil[i]*Weights[i]; 
 
  return(sum);
}


double WENO5_Reconstruction_RightY(Point point, scalar *Xl, double gradR){

 scalar X = Xl[0];
 double T_Stencil[3], Beta[3], Weights[3];
 int i;
 double sum;

 double gammaR[3] = {3./10.,3./5.,1./10.};
 T_Stencil[0] = -1.*X[0,-2]/6. + 5.*X[0,-1]/6. + 1.*X[]/3.;
 Beta[0]      = 13.*sq(X[0,-2] - 2.*X[0,-1] + X[])/12. + sq(X[0,-2] - 4.*X[0,-1] + 3.*X[])/4.;
 T_Stencil[1] =  1.*X[0,-1]/3. + 5.*X[]/6. - 1.*X[0,1]/6.;
 Beta[1]      = 13.*sq(X[0,-1] - 2.*X[] + X[0,1])/12. + sq(X[0,-1]-X[0,1])/4.;
 T_Stencil[2] = 11.*X[]/6. - 7.*X[0,1]/6. + (X[]+2.*Delta*gradR)/3.;
 Beta[2]      = 13.*sq(2.*X[] - 2.*X[0,1] + 2.*Delta*gradR)/12. + sq(4.*X[] - 4.*X[0,1] + 2.*Delta*gradR)/4.;
 for(i=0; i<=2; i++)
   Weights[i] = gammaR[i]/sq(epsilon + Beta[i]);
  
 sum = Weights[0] + Weights[1] + Weights[2];
 for(i=0; i<=2; i++)
   Weights[i] /= sum;

 sum = 0;
 for(i=0;i<=2;i++)
   sum += T_Stencil[i]*Weights[i]; 
 
 return(sum);
}



double WENO5_Reconstruction_LeftX(Point point, scalar *Xl, double gradL){

  scalar X = Xl[0];
  double T_Stencil[3], Beta[3], Weights[3];
  int i;
  double sum;

  double gammaL[3] = {1./10.,3./5.,3./10.};
  T_Stencil[0]  =   1.*(X[-1,0] - 2.*Delta*gradL)/3. - 7.*X[-2,0]/6. + 11.*X[-1,0]/6.;
  Beta[0]       =  13.*sq(X[-1,0] - 2.*Delta*gradL - 2.*X[-2,0] + X[-1,0])/12. + sq(X[-1,0] - 2.*Delta*gradL - 4.*X[-2,0] + 3.*X[-1,0])/4.;
  T_Stencil[1]  =  -1.*X[-2,0]/6. + 5.*X[-1,0]/6. + 1.*X[]/3.;
  Beta[1]       =  13.*sq(X[-2,0] - 2.*X[-1,0] + X[])/12. + sq(X[-2,0]-X[])/4.;
  T_Stencil[2]  =   1.*X[-1,0]/3. + 5.*X[]/6. - 1.*X[1,0]/6.;
  Beta[2]       =  13.*sq(X[-1,0] - 2.*X[] + X[1,0])/12. + sq(3.*X[-1,0] - 4.*X[] + X[1,0])/4.;
     
  for(i=0; i<=2; i++)
    Weights[i] = gammaL[i]/sq(epsilon + Beta[i]);
  
  sum = Weights[0] + Weights[1] + Weights[2];
  for(i=0; i<=2; i++)
    Weights[i] /= sum;

  sum = 0.;
  for(i=0;i<=2;i++)
   sum += T_Stencil[i]*Weights[i]; 
 
  return(sum);
}


double WENO5_Reconstruction_RightX(Point point, scalar *Xl, double gradR){

 scalar X = Xl[0];
 double T_Stencil[3], Beta[3], Weights[3];
 int i;
 double sum;

 double gammaR[3] = {3./10.,3./5.,1./10.};
 T_Stencil[0] = -1.*X[-2,0]/6. + 5.*X[-1,0]/6. + 1.*X[]/3.;
 Beta[0]      = 13.*sq(X[-2,0] - 2.*X[-1,0] + X[])/12. + sq(X[-2,0] - 4.*X[-1,0] + 3.*X[])/4.;
 T_Stencil[1] =  1.*X[-1,0]/3. + 5.*X[]/6. - 1.*X[1,0]/6.;
 Beta[1]      = 13.*sq(X[-1,0] - 2.*X[] + X[1,0])/12. + sq(X[-1,0]-X[1,0])/4.;
 T_Stencil[2] = 11.*X[]/6. - 7.*X[1,0]/6. + (X[]+2.*Delta*gradR)/3.;
 Beta[2]      = 13.*sq(2.*X[] - 2.*X[1,0] + 2.*Delta*gradR)/12. + sq(4.*X[] - 4.*X[1,0] + 2.*Delta*gradR)/4.;
 for(i=0; i<=2; i++)
   Weights[i] = gammaR[i]/sq(epsilon + Beta[i]);
  
 sum = Weights[0] + Weights[1] + Weights[2];
 for(i=0; i<=2; i++)
   Weights[i] /= sum;

 sum = 0;
 for(i=0;i<=2;i++)
   sum += T_Stencil[i]*Weights[i]; 
 
 return(sum);
}
   


double * QuadratureValues2DX(Point point, vector * WENO_AvgV ){

 face vector WENO_Avg = WENO_AvgV[0];
 int i,j,k;
 double U_Quad_Stencil[4][4],Weights[4][3];
 double Beta[3];
 double sum;
 
 double StencilQuadWeights[4][3][3]=   
        {
            {
            {(2.  - 3.*sqrt(15))/60., (-4. + 12.*sqrt(15))/60., (62. - 9.*sqrt(15))/60.},
            {(2.  + 3.*sqrt(15))/60.,         56./60.         , ( 2. - 3.*sqrt(15))/60.},
            {(62. + 9.*sqrt(15))/60., (-4. - 12.*sqrt(15))/60., ( 2. + 3.*sqrt(15))/60.}
            },
            {
            {-1./24. , 1./12.  , 23./24. },
            {-1./24. , 26./24. , -1./24. },
            {23./24. ,  2./24. , -1./24. }
            },
            {
            {-1./24. , 1./12.  , 23./24. },
            {-1./24. , 26./24. , -1./24. },
            {23./24. ,  2./24. , -1./24. }
            },
            {
            {(2.  + 3.*sqrt(15))/60., (-4. - 12.*sqrt(15))/60., (62. + 9.*sqrt(15))/60.},
            {(2.  - 3.*sqrt(15))/60.,         56./60.         , ( 2. + 3.*sqrt(15))/60.},
            {(62. - 9.*sqrt(15))/60., (-4. + 12.*sqrt(15))/60., ( 2. - 3.*sqrt(15))/60.}
            },
        };

 double GammaQuad[4][3] = 
        {{(1008. + 71.*sqrt(15))/5240.,403./655.,(1008. - 71.*sqrt(15))/5240.},
         {9./214.,196./214.,9./214.},
         {9./67.,49./67.,9./67.},
         {(1008. - 71.*sqrt(15))/5240.,403./655.,(1008. + 71.*sqrt(15))/5240.}
        };

  for(i=0;i<=3;i++)
    for(j=0;j<=3;j++)
       U_Quad_Stencil[i][j] = 0.;
 
  for(i=0;i<=3;i++)
    for(j=0;j<=2;j++)
      for(k=0;k<=2;k++)
         U_Quad_Stencil[i][j] += StencilQuadWeights[i][j][k] * WENO_Avg.x[0,j+k-2];
        
  Beta[0] = 13.*sq(WENO_Avg.x[0,-2] - 2.*WENO_Avg.x[0,-1] + WENO_Avg.x[0,0])/12. + sq(WENO_Avg.x[0,-2] - 4.*WENO_Avg.x[0,-1] + 3.*WENO_Avg.x[0,0])/4.;
  Beta[1] = 13.*sq(WENO_Avg.x[0,-1] - 2.*WENO_Avg.x[0,0]  + WENO_Avg.x[0,1])/12. + sq(WENO_Avg.x[0,-1]-WENO_Avg.x[0,1])/4.;
  Beta[2] = 13.*sq(WENO_Avg.x[0,0]  - 2.*WENO_Avg.x[0,1]  + WENO_Avg.x[0,2])/12. + sq(3.*WENO_Avg.x[0,0] - 4.*WENO_Avg.x[0,1] + WENO_Avg.x[0,2])/4.;

  for(i=0;i<=3;i++){
    for(j=0;j<=2;j++)
       Weights[i][j] = GammaQuad[i][j]/sq(epsilon + Beta[j]);          
    sum = Weights[i][0] + Weights[i][1] + Weights[i][2]; 
    for(j=0;j<=2;j++)
      Weights[i][j] /= sum;
   }
   
  for(i=0;i<=3;i++)
    for(j=0;j<=2;j++)
         U_Quad_Stencil[i][3] += Weights[i][j]*U_Quad_Stencil[i][j];

  static double QuadValues[3];
  QuadValues[0] = U_Quad_Stencil[0][3];
  QuadValues[1] = 214.*U_Quad_Stencil[1][3]/80. - 67.*U_Quad_Stencil[2][3]/40.;
  QuadValues[2] = U_Quad_Stencil[3][3];

  return(QuadValues);
}

double * QuadratureValues2DY(Point point, vector * WENO_AvgV ){

 face vector WENO_Avg = WENO_AvgV[0];
 int i,j,k;
 double U_Quad_Stencil[4][4],Weights[4][3];
 double Beta[3];
 double sum;
 
 double StencilQuadWeights[4][3][3]=   
        {
            {
            {(2.  - 3.*sqrt(15))/60., (-4. + 12.*sqrt(15))/60., (62. - 9.*sqrt(15))/60.},
            {(2.  + 3.*sqrt(15))/60.,         56./60.         , ( 2. - 3.*sqrt(15))/60.},
            {(62. + 9.*sqrt(15))/60., (-4. - 12.*sqrt(15))/60., ( 2. + 3.*sqrt(15))/60.}
            },
            {
            {-1./24. , 1./12.  , 23./24. },
            {-1./24. , 26./24. , -1./24. },
            {23./24. ,  2./24. , -1./24. }
            },
            {
            {-1./24. , 1./12.  , 23./24. },
            {-1./24. , 26./24. , -1./24. },
            {23./24. ,  2./24. , -1./24. }
            },
            {
            {(2.  + 3.*sqrt(15))/60., (-4. - 12.*sqrt(15))/60., (62. + 9.*sqrt(15))/60.},
            {(2.  - 3.*sqrt(15))/60.,         56./60.         , ( 2. + 3.*sqrt(15))/60.},
            {(62. - 9.*sqrt(15))/60., (-4. + 12.*sqrt(15))/60., ( 2. - 3.*sqrt(15))/60.}
            },
        };

 double GammaQuad[4][3] = 
        {{(1008. + 71.*sqrt(15))/5240.,403./655.,(1008. - 71.*sqrt(15))/5240.},
         {9./214.,196./214.,9./214.},
         {9./67.,49./67.,9./67.},
         {(1008. - 71.*sqrt(15))/5240.,403./655.,(1008. + 71.*sqrt(15))/5240.}
        };

  for(i=0;i<=3;i++)
    for(j=0;j<=3;j++)
       U_Quad_Stencil[i][j] = 0.;
 
  for(i=0;i<=3;i++)
    for(j=0;j<=2;j++)
      for(k=0;k<=2;k++)
         U_Quad_Stencil[i][j] += StencilQuadWeights[i][j][k] * WENO_Avg.y[j+k-2,0];
        
  Beta[0] = 13.*sq(WENO_Avg.y[-2,0] - 2.*WENO_Avg.y[-1,0] + WENO_Avg.y[0,0])/12. + sq(WENO_Avg.y[-2,0] - 4.*WENO_Avg.y[-1,0] + 3.*WENO_Avg.y[0,0])/4.;
  Beta[1] = 13.*sq(WENO_Avg.y[-1,0] - 2.*WENO_Avg.y[0,0]  + WENO_Avg.y[1,0])/12. + sq(WENO_Avg.y[-1,0]-WENO_Avg.y[1,0])/4.;
  Beta[2] = 13.*sq(WENO_Avg.y[0,0]  - 2.*WENO_Avg.y[1,0]  + WENO_Avg.y[2,0])/12. + sq(3.*WENO_Avg.y[0,0] - 4.*WENO_Avg.y[1,0] + WENO_Avg.y[2,0])/4.;

  for(i=0;i<=3;i++){
    for(j=0;j<=2;j++)
       Weights[i][j] = GammaQuad[i][j]/sq(epsilon + Beta[j]);          
    sum = Weights[i][0] + Weights[i][1] + Weights[i][2]; 
    for(j=0;j<=2;j++)
      Weights[i][j] /= sum;
   }
   
  for(i=0;i<=3;i++)
    for(j=0;j<=2;j++)
         U_Quad_Stencil[i][3] += Weights[i][j]*U_Quad_Stencil[i][j];

  static double QuadValues[3];
  QuadValues[0] = U_Quad_Stencil[0][3];
  QuadValues[1] = 214.*U_Quad_Stencil[1][3]/80. - 67.*U_Quad_Stencil[2][3]/40.;
  QuadValues[2] = U_Quad_Stencil[3][3];

  return(QuadValues);
}


double RiemannFlux (double LeftValue, double RightValue, double velocity){
 
   if(velocity >= 0)
      return (LeftValue);
   else
      return (RightValue);
}

void Weno_Flux_Order5(scalar *Xl, vector *ul, vector *Wenol){

  scalar X = Xl[0];
  vector u = ul[0];
  face vector ufavg[];
  face vector Weno = Wenol[0];

  vector grad_X[], grad_U[];
  foreach(){
    foreach_dimension(){
       grad_X.x[] = (X[1] - X[-1])/(2.*Delta);
       grad_U.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
      }
   }

  boundary((scalar *){grad_X, grad_U});

  #if dimension == 1

   foreach_face(){
     ufavg.x[] =   RiemannFlux( WENO5_Reconstruction_LeftX(point,{u.x},grad_U.x[-2]), WENO5_Reconstruction_RightX(point,{u.x},grad_U.x[1]), u.x[]  ); 
     Weno.x[] = ufavg.x[]*(RiemannFlux( WENO5_Reconstruction_LeftX(point,{X},grad_X.x[-2])  , WENO5_Reconstruction_RightX(point,{X},grad_X.x[1])  , ufavg.x[] ));
    }   

   #if TREE
     boundary_flux({Weno});
   #endif

  #elif dimension == 2

    face vector X_wenoavgL[], X_wenoavgR[];

    foreach_face(x){      
         X_wenoavgL.x[] = WENO5_Reconstruction_LeftX(point,{X},grad_X.x[-2,0]); 
         X_wenoavgR.x[] = WENO5_Reconstruction_RightX(point,{X},grad_X.x[1,0]);
         ufavg.x[] = RiemannFlux( WENO5_Reconstruction_LeftX(point,{u.x},grad_U.x[-2,0]), WENO5_Reconstruction_RightX(point,{u.x},grad_U.x[1,0]), u.x[]  ); 
    }

    foreach_face(y){      
         X_wenoavgL.y[] = WENO5_Reconstruction_LeftY(point,{X},grad_X.y[0,-2]); 
         X_wenoavgR.y[] = WENO5_Reconstruction_RightY(point,{X},grad_X.y[0,1]);
         ufavg.y[] = RiemannFlux( WENO5_Reconstruction_LeftY(point,{u.y},grad_U.y[0,-2]), WENO5_Reconstruction_RightY(point,{u.y},grad_U.y[0,1]), u.y[]  ); 
    }

    boundary((scalar *){X_wenoavgL,X_wenoavgR});

  
    double *XQL, *XQR;  // These 3 element arrays will hold the 3 quadtrature values obtained.

    foreach_face(x){
      XQL = QuadratureValues2DX(point,{X_wenoavgL});
      XQR = QuadratureValues2DX(point,{X_wenoavgR});
      Weno.x[] = ufavg.x[]*((5./18.)*RiemannFlux(*XQL,*XQR,ufavg.x[]) + (8./18.)*RiemannFlux(*(XQL+1),*(XQR+1),ufavg.x[]) + (5./18.)*RiemannFlux(*(XQL+2),*(XQR+2),ufavg.x[]));
    }

    foreach_face(y){
      XQL = QuadratureValues2DY(point,{X_wenoavgL});
      XQR = QuadratureValues2DY(point,{X_wenoavgR});
      Weno.y[] = ufavg.y[]*((5./18.)*RiemannFlux(*XQL,*XQR,ufavg.y[]) + (8./18.)*RiemannFlux(*(XQL+1),*(XQR+1),ufavg.y[]) + (5./18.)*RiemannFlux(*(XQL+2),*(XQR+2),ufavg.y[]));
    }

   #if TREE
     boundary_flux({Weno});
   #endif

  #endif
}