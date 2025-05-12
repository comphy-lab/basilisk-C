/**
## VOF to level-set mapping
The following functions are meant to obtain a level-set distance
field from the given volume fraction field. */

void vof_to_LS (scalar vof, scalar ls, vector gLSc)
{

  /**
  Redistancing parameters. Procedure based on Albadawi paper on 2013.*/
  int NCorr = 15;

  scalar lsSgn[],HS[];
  
  foreach(){
    HS[] = (2*vof[] - 1.)*Delta*0.75;
    lsSgn[] = (ls[] > 1e-6 ? 1. : ls[] < -1e-6 ? -1. : ls[]);
  }
  
    /**
  Smoothing the original step function:  */
  foreach()
    ls[] = (8.*HS[] + 
	    4.*(HS[0,0,1] + HS[0,0,-1] + HS[0,1,0] + HS[0,-1,0] + HS[1,0,0] + HS[-1,0,0] )+
            2.*(HS[0,1,1] + HS[0,1,-1] + HS[0,-1,1] + HS[0,-1,-1] +
            HS[1,0,1] + HS[-1,0,1] + HS[1,0,-1] + HS[-1,0,-1] +
            HS[1,1,0] + HS[-1,1,0] + HS[1,-1,0] + HS[-1,-1,0]) +
	    HS[-1,-1,1] + HS[1,-1,1] + HS[1,1,1] + HS[-1,1,1] + 
	    HS[-1,-1,-1] + HS[1,-1,-1] + HS[1,1,-1] + HS[-1,1,-1] )/64.;
  
  boundary ((scalar *){ls});
  /**
  Apply redistancing procedure *NCorr* iterations.*/

  int iters = 0;
  double minGrad = 0.;
  double maxGrad = 0.;
  face vector gLSf[];
     
  while(iters < NCorr){
    iters++;
    minGrad = HUGE;
    maxGrad = -HUGE;

  /**
  Gradient computation on faces.*/
  
    foreach_face()
      gLSf.x[] = (ls[] - ls[-1])/Delta;
    boundary_flux({gLSf});

    foreach(){
      double curGrad = 0;
      
  /**
  Gradient interpolation to cell centers.*/
        foreach_dimension(){
        curGrad += sq((gLSf.x[] + gLSf.x[1])/2.);
        gLSc.x[] = (gLSf.x[] + gLSf.x[1])/2.;
      }
      
  /**
  Update level-set function.*/
      ls[] = ls[] + lsSgn[]*(1. - sqrt(curGrad))*Delta*0.1;
      boundary ((scalar *){ls});
      
  /**
  Get maximum and minimum values for the gradient magnitude.*/
      if(sqrt(curGrad) < minGrad)
        minGrad = sqrt(curGrad);
      if(sqrt(curGrad) > maxGrad)
        maxGrad = sqrt(curGrad);
    }

  fprintf(stderr,"VOF_TO_LS Iters: %d \n",iters);
  fprintf(stderr,"maxGrad: %g \n",maxGrad);
  fprintf(stderr,"minGrad: %g \n\n",minGrad);

  }

}
