/**
# Adaptivity Header File - Has the 5th order Prolongation function modules
*/

double stencil_weightP[3][5], stencil_weightM[3][5];

void Prolongation_Weight_Initialization(){ 	 	

  int i,j;
  stencil_weightP[0][0] =  0.0103444235735617; 	 
  stencil_weightP[0][1] = -0.09792213521550909;
  stencil_weightP[0][2] =  1.107094656676454; 	
  stencil_weightP[0][3] = -0.01815143469025667; 	
  stencil_weightP[0][4] = -0.001365510344249945;

  stencil_weightP[1][0] =  0.02568359375; 	 
  stencil_weightP[1][1] = -0.188671875; 	
  stencil_weightP[1][2] =  1.026497395833333;
  stencil_weightP[1][3] =  0.1602864583333333;
  stencil_weightP[1][4] = -0.02379557291666666;
  
  stencil_weightP[2][0] =  0.0329368264264383;
  stencil_weightP[2][1] = -0.2189528647844909;
  stencil_weightP[2][2] =  0.8505095099902126;
  stencil_weightP[2][3] =  0.3804431013569233;
  stencil_weightP[2][4] = -0.04493657298908339;

  for(i=0; i<=2; i++)        
     for(j=0; j<=4; j++)          
        stencil_weightM[i][j] = stencil_weightP[2-i][4-j];
        
}

double * Prolongation_Quadrature_Order5 (double *T, int marker){

     int i,j;
     static double T_Quad[3];
    
     if(marker == 1){
        for(i=0; i<=2; i++){
            T_Quad[i] = 0.;
            for(j=0; j<=4; j++)
                  T_Quad[i] += stencil_weightP[i][j]*(*(T+j));
         }
     }
     
     else if(marker == -1){
         for(i=0; i<=2; i++){
            T_Quad[i] = 0.;
            for(j=0; j<=4; j++)
                  T_Quad[i] += stencil_weightM[i][j]*(*(T+j));
         }
     }
    
     return (T_Quad);  
}


#if dimension == 1

double Order5_Refine1D (Point point, scalar X){
     double T_hold[5];
     double * QuadV;
     for(int i = -2; i<=2; i++)
        T_hold[i+2] = coarse(X,i);
     QuadV = Prolongation_Quadrature_Order5 (T_hold, child.x);
     return(5.*(*QuadV)/18. + 8.*(*(QuadV+1))/18. + 5.*(*(QuadV+2))/18.);
}


#elif dimension == 2

static inline double Order5_Refine2D (Point point, scalar X){
     double T_hold[5];
     double * QuadPtr;
     double QuadArr1[5][3];
     double QuadArr2[3];
     int i,j;

     for(j = -2; j<=2; j++){
        for(i = -2; i<=2; i++)
            T_hold[i+2] = coarse(X,i,j);
        QuadPtr = Prolongation_Quadrature_Order5 (T_hold, child.x);
        
        for(i = 0; i<=2; i++)
          QuadArr1[j+2][i] = *(QuadPtr+i);         
     }

     for(j = 0; j<=2; j++){
        for(i = 0; i<=4; i++)
            T_hold[i] = QuadArr1[i][j];
        QuadPtr = Prolongation_Quadrature_Order5 (T_hold, child.y);
        QuadArr2[j] = 5.*(*(QuadPtr))/18. + 8.*(*(QuadPtr+1))/18. + 5.*(*(QuadPtr+2))/18.;
     }
  
    return(5.*QuadArr2[0]/18. + 8.*QuadArr2[1]/18. + 5.*QuadArr2[2]/18. );
}

#elif dimension == 3

static inline double Order5_Refine3D (Point point, scalar X){
     double T_hold[5];
     double * QuadPtr;
     double QuadArr1[5][5][3];
     double QuadArr2[5][3][3];
     double QuadArr3[3][3][3];

     int i,j,k;

     for(k= -2; k<=2; k++){
        for(j = -2; j<=2; j++){
           
           for(i = -2; i<=2; i++)
               T_hold[i+2] = coarse(X,i,j,k);

           QuadPtr = Prolongation_Quadrature_Order5 (T_hold, child.x);

           for(i = 0; i<=2; i++)
              QuadArr1[j+2][k+2][i] = *(QuadPtr+i);         
         }
     }

     for(i=0; i<=2; i++){
         
        for(j=0;j<=4;j++){
          
            for(k=0;k<=4;k++)
               T_hold[k] = QuadArr1[j][k][i];
            QuadPtr = Prolongation_Quadrature_Order5 (T_hold, child.z);
 
            for(k=0; k<=2; k++) 
               QuadArr2[j][k][i] = *(QuadPtr+k); 
         }

      }
             
     for(i=0; i<=2; i++){
         for(k=0; k<=2; k++){
             
            for(j=0;j<=4;j++)
               T_hold[j] = QuadArr2[j][k][i];
            QuadPtr = Prolongation_Quadrature_Order5 (T_hold, child.y);
          
            for(j=0;j<=2;j++)
               QuadArr3[j][k][i] = *(QuadPtr+j);

         }
      }

    double QSum[3] = {5./18., 8./18., 5./18. };
    double QuadSummation = 0.;
    for(k=0; k<=2; k++)
      for(j=0; j<=2; j++)
         for(i=0; i<=2; i++)
                QuadSummation += QSum[j]*QSum[k]*QSum[i]*QuadArr3[j][k][i];   
  
    return(QuadSummation);
}

#endif


static inline void refine_order5 (Point point, scalar s)
{
  #if dimension == 1
    foreach_child()
       s[] = Order5_Refine1D (point, s);
  #elif dimension == 2
    foreach_child()
       s[] = Order5_Refine2D (point, s);
  #elif dimension == 3
    foreach_child()
       s[] = Order5_Refine3D (point, s);
  #endif
}