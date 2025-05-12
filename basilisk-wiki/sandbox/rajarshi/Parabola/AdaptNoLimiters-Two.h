/**
# Prolongation function 5th order with 2 quadrature points.
*/

double stencil_weightP[2][5], stencil_weightM[2][5];

void Prolongation_Weight_Initialization(){ 	 	

  int i,j;
  stencil_weightP[0][0] =    863./36800. + 17.*sqrt(3)/3456.; 	 
  stencil_weightP[0][1] =  -6327./36800. - 5.*sqrt(3)/192.;			
  stencil_weightP[0][2] =  36803./36800. - 1.*sqrt(3)/18.; 	
  stencil_weightP[0][3] =   6323./36800. + 149.*sqrt(3)/1728.; 	
  stencil_weightP[0][4] =   -431./18400. - 11.*sqrt(3)/1152.;
 
  stencil_weightP[1][0] =    863./36800. - 17.*sqrt(3)/3456.; 	 
  stencil_weightP[1][1] =  -6327./36800. + 5.*sqrt(3)/192.;			
  stencil_weightP[1][2] =  36803./36800. + 1.*sqrt(3)/18.; 	
  stencil_weightP[1][3] =   6323./36800. - 149.*sqrt(3)/1728.; 	
  stencil_weightP[1][4] =   -431./18400. + 11.*sqrt(3)/1152.;


  for(i=0; i<=1; i++)        
     for(j=0; j<=4; j++)          
        stencil_weightM[i][j] = stencil_weightP[1-i][4-j];
        
}

double * Prolongation_Quadrature_Order5 (double *T, int marker){

     int i,j;
     static double T_Quad[2];
    
     if(marker == 1){
        for(i=0; i<=1; i++){
            T_Quad[i] = 0.;
            for(j=0; j<=4; j++)
                  T_Quad[i] += stencil_weightP[i][j]*(*(T+j));
         }
     }
     
     else if(marker == -1){
         for(i=0; i<=1; i++){
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
     return(1.*(*QuadV)/2. + 1.*(*(QuadV+1))/2.);
}


#elif dimension == 2

static inline double Order5_Refine2D (Point point, scalar X){
     double T_hold[5];
     double * QuadPtr;
     double QuadArr1[5][2];
     double QuadArr2[2];
     int i,j;

     for(j = -2; j<=2; j++){
        for(i = -2; i<=2; i++)
            T_hold[i+2] = coarse(X,i,j);
        QuadPtr = Prolongation_Quadrature_Order5 (T_hold, child.x);
        
        for(i = 0; i<=1; i++)
          QuadArr1[j+2][i] = *(QuadPtr+i);         
     }

     for(j = 0; j<=1; j++){
        for(i = 0; i<=4; i++)
            T_hold[i] = QuadArr1[i][j];
        QuadPtr = Prolongation_Quadrature_Order5 (T_hold, child.y);
        QuadArr2[j] = 1.*(*(QuadPtr))/2. + 1.*(*(QuadPtr+1))/2.;
     }
  
    return(1.*QuadArr2[0]/2. + 1.*QuadArr2[1]/2.);
}

#elif dimension == 3

static inline double Order5_Refine3D (Point point, scalar X){
     double T_hold[5];
     double * QuadPtr;
     double QuadArr1[5][5][2];
     double QuadArr2[5][2][2];
     double QuadArr3[2][2][2];

     int i,j,k;

     for(k= -2; k<=2; k++){
        for(j = -2; j<=2; j++){
           
           for(i = -2; i<=2; i++)
               T_hold[i+2] = coarse(X,i,j,k);

           QuadPtr = Prolongation_Quadrature_Order5 (T_hold, child.x);

           for(i = 0; i<=1; i++)
              QuadArr1[j+2][k+2][i] = *(QuadPtr+i);         
         }
     }

     for(i=0; i<=1; i++){
         
        for(j=0;j<=4;j++){
          
            for(k=0;k<=4;k++)
               T_hold[k] = QuadArr1[j][k][i];
            QuadPtr = Prolongation_Quadrature_Order5 (T_hold, child.z);
 
            for(k=0; k<=1; k++) 
               QuadArr2[j][k][i] = *(QuadPtr+k); 
         }

      }
             
     for(i=0; i<=1; i++){
         for(k=0; k<=1; k++){
             
            for(j=0;j<=4;j++)
               T_hold[j] = QuadArr2[j][k][i];
            QuadPtr = Prolongation_Quadrature_Order5 (T_hold, child.y);
          
            for(j=0;j<=1;j++)
               QuadArr3[j][k][i] = *(QuadPtr+j);

         }
      }

    double QSum[2] = {1./2., 1./2.};
    double QuadSummation = 0.;
    for(k=0; k<=1; k++)
      for(j=0; j<=1; j++)
         for(i=0; i<=1; i++)
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
