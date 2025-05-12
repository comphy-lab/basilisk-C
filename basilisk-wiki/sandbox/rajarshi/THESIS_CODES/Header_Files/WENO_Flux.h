/**
#WENO TRACER FLUX SCHEME

This header file implements the full WENO scheme for tracer-advection ( following the work of liu-et-al(1994) ). This header
file must be used when there are discontinuities in the advecting tracer field. For continuous tracer-fields we use the 
computationally cheaper version O5_Flux.h
*/

#define BGHOSTS 2
#define epsilon 1E-06

extern void face_velocity (face vector u, double t);

foreach_dimension(){

  double weno5_left_x (Point point, scalar X, double gradL){

    double T_Stencil[3], Beta[3], Weights[3];
    int i;
    double sum;

    double gammaL[3] = { 1./10., 3./5., 3./10. };

    T_Stencil[0]  =   1.*( X[-1] - 2.*Delta*gradL )/3.
                    - 7.*X[-2]/6. + 11.*X[-1]/6.;
    Beta[0]       =  13.*sq( X[-1] - 2.*Delta*gradL - 2.*X[-2] + X[-1] )/12.
                     + sq( X[-1] - 2.*Delta*gradL - 4.*X[-2] + 3.*X[-1] )/4.;

    T_Stencil[1]  =  -1.*X[-2]/6. + 5.*X[-1]/6. + 1.*X[]/3.;
    Beta[1]       =  13.*sq( X[-2] - 2.*X[-1] + X[] )/12. + sq ( X[-2]-X[] )/4.;

    T_Stencil[2]  =   1.*X[-1]/3. + 5.*X[]/6. - 1.*X[1]/6.;
    Beta[2]       =  13.*sq( X[-1] - 2.*X[] + X[1] )/12.
                     + sq( 3.*X[-1] - 4.*X[] + X[1] )/4.;
     
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


  double weno5_right_x (Point point, scalar X, double gradR){

    double T_Stencil[3], Beta[3], Weights[3];
    int i;
    double sum;

    double gammaR[3] = {3./10., 3./5., 1./10.};
    
    T_Stencil[0] = -1.*X[-2]/6. + 5.*X[-1]/6. + 1.*X[]/3.;
    Beta[0]      = 13.*sq( X[-2] - 2.*X[-1] + X[] )/12.
                   + sq( X[-2] - 4.*X[-1] + 3.*X[] )/4.;
    
    T_Stencil[1] =  1.*X[-1]/3. + 5.*X[]/6. - 1.*X[1]/6.;
    Beta[1]      = 13.*sq( X[-1] - 2.*X[] + X[1] )/12. + sq( X[-1]-X[1] )/4.;
    
    T_Stencil[2] = 11.*X[]/6. - 7.*X[1]/6. + ( X[]+2.*Delta*gradR )/3.;
    Beta[2]      = 13.*sq( 2.*X[] - 2.*X[1] + 2.*Delta*gradR )/12.
                   + sq( 4.*X[] - 4.*X[1] + 2.*Delta*gradR )/4.;

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

  #if dimension > 1
  static void quadrature3p_2D_y (Point point, face vector avg, double * q) {

    static const double StencilQuadWeights[4][3][3] = { 
      {
	{ -0.160315833977,  0.707930002575,  0.452385831402 }, 	
	{  0.226982500644,  0.933333333333, -0.160315833977 },	
	{  1.61428083526 , -0.841263335908,  0.226982500644 }	
      },

      {
	{ -0.0416666666667, 0.0833333333333,  0.958333333333  }, 	
	{ -0.0416666666667, 1.08333333333  , -0.0416666666667 }, 	
	{  0.958333333333 , 0.0833333333333, -0.0416666666667 }
      }, 	

      {
	{ -0.0416666666667, 0.0833333333333,  0.958333333333   },  	
	{ -0.0416666666667, 1.08333333333  , -0.0416666666667  },	
	{  0.958333333333 , 0.0833333333333, -0.0416666666667  }
      }, 	

      {
	{  0.226982500644, -0.841263335908,  1.61428083526    }, 	
	{ -0.160315833977,  0.933333333333,  0.226982500644   },	
	{  0.452385831402,  0.707930002575, -0.160315833977   }
      } 	
    }; 

    static const double GammaQuad[4][3] = {
      {0.244843858317 , 0.615267175573, 0.139888966111  },
      {0.0420560747664, 0.915887850467, 0.0420560747664 },
      {0.134328358209 , 0.731343283582, 0.134328358209  }
    };
 
    double U_Quad_Stencil[4][4];
    for (int i=0; i<=3; i++) {
      for (int j=0; j<=2; j++) {
	U_Quad_Stencil[i][j] = 0.; 
	for (int k=0;k<=2;k++)
	  U_Quad_Stencil[i][j] += StencilQuadWeights[i][j][k] * avg.y[j+k-2];
      }
    }
  
    double Beta[3];
    Beta[0] = 13.*sq ( avg.y[-2] - 2.*avg.y[-1] + avg.y[]  )/12.
                + sq (    avg.y[-2] - 4.*avg.y[-1] + 3.*avg.y[] )/4.;
		
    Beta[1] = 13.*sq ( avg.y[-1] - 2.*avg.y[]   + avg.y[1] )/12.
                + sq (    avg.y[-1] -    avg.y[1] )/4.;
		
    Beta[2] = 13.*sq ( avg.y[]   - 2.*avg.y[1]  + avg.y[2] )/12.
                + sq ( 3.*avg.y[]   - 4.*avg.y[1]  +    avg.y[2] )/4.;

 
    double Weights[4][3];
    double sum; 
    for ( int i = 0; i <= 3; i++){
      for ( int j = 0; j <= 2; j++)
	Weights[i][j] = GammaQuad[i][j]/sq(epsilon + Beta[j]);          
      sum = Weights[i][0] + Weights[i][1] + Weights[i][2]; 
      for(int j=0;j<=2;j++)
	Weights[i][j] /= sum;
    }
   
    for ( int i = 0; i <= 3; i++)
      for ( int j = 0; j <= 2; j++)
	U_Quad_Stencil[i][3] += Weights[i][j] * U_Quad_Stencil[i][j];

    q[0] = U_Quad_Stencil[0][3];
    q[1] = 214.*U_Quad_Stencil[1][3]/80. - 67.*U_Quad_Stencil[2][3]/40.;
    q[2] = U_Quad_Stencil[3][3];

  }
  #endif
}


static void fluxes (scalar f, face vector u, face vector flux)
{
  vector grad[];
  foreach()
    foreach_dimension()
    grad.x[] = (f[1] - f[-1])/(2.*Delta);
  boundary ((scalar *){grad});

  foreach_face() {
    if (u.x[] >= 0)
      flux.x[] = weno5_left_x (point, f, grad.x[-2]);
    else
      flux.x[] = weno5_right_x (point, f, grad.x[1]);
  }
}


void tracer_fluxes (scalar * fl, double t, scalar * gl)
{
  face_velocity (u, t);

  scalar f, g;  
  for (f,g in fl,gl) {

#if dimension == 1 
    face vector flux[];
    fluxes (f, u, flux);
    foreach_face()
      flux.x[] *= u.x[];

#elif dimension == 2
    face vector line_avg[];
    fluxes (f, u, line_avg);
    boundary ((scalar *) {line_avg});
    
    face vector flux[];
    foreach_face() {
      double fq[3];
      quadrature3p_2D_x (point, line_avg, fq);
      double uq[3];
      quadrature3p_2D_x (point, u, uq);
      flux.x[] = (5.*fq[0]*uq[0] + 8.*fq[1]*uq[1] + 5.*fq[2]*uq[2])/18.;
    }
#endif
    boundary_flux ({flux});

    foreach() {
      g[] = 0.;
      foreach_dimension()
	g[] += (flux.x[] - flux.x[1])/Delta;
    }
  }
}
