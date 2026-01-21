event AC_equation (i++) {
#ifdef two_phase
{ 
    vector nphi[],gradT[];
	foreach() {
		//Tm_eff[] = Tm;	
		foreach_dimension() {
			nphi.x[]  = center_gradient(phi1);
			gradT.x[] = center_gradient(T);
		}
	}

	face vector temp_la[];
	foreach_face(){
		temp_la.x[] = face_gradient_x(phi1,0);
    }
	
	scalar grad_magnitude[];
	foreach() {
		double temp = 0.0;
		foreach_dimension() {
			temp += sq(nphi.x[]);
		}
		grad_magnitude[] = sqrt(temp);
	}
	
	foreach_dimension() {
		foreach() {
		  double keep=1.e-6;
		  if (grad_magnitude[] >= keep){
			nphi.x[] = phi1[]*(1.0-phi1[])/sqrt(2.0)/Cn*nphi.x[]/grad_magnitude[];
		  }
		  else{
			nphi.x[] = 0.0;  
		  }  
		}
	}
	
	face vector nphiflux[];
	foreach_face() {
		nphiflux.x[] = (nphi.x[]+nphi.x[-1])/2.0;	
	}	
	
	scalar RHS_AC[];
	foreach() {
		double la_f=0.0, temp=0.0;
		foreach_dimension() {
			temp += (fm.x[1]*nphiflux.x[1] - fm.x[]*nphiflux.x[])/(Delta*max(cm[],1.e-15));
			la_f += (fm.x[1]*temp_la.x[1]  - fm.x[]*temp_la.x[]) /(Delta*max(cm[],1.e-15));
		}	
		RHS_AC[] = la_f-temp;		
	}

scalar phit[];
for(int kk = 0; kk< 1; kk++){
///%Prediction and correction step----------------------------------------------------------------------------------------------------------------------------------------
	foreach() {
		double temp;
		temp = RHS_AC[]
			 +1.0/Cn/Cn*A_AC*(T[]-Tm_eff[])*phi1[]*(1.0-phi1[]);

		phit[] = phi1[]+ dt*D_AC*temp;
		phit[] = clamp (phit[], 0., 1.);
	}

	boundary({phit});



//******************************-----------start Tm_eff_AC---------------------------------------------------------------------------------------------------------------------
	vector nphit[],gradphit[];
	foreach() {
		foreach_dimension() {
			gradphit.x[] = center_gradient(phit);
		}
	}
	
	scalar grad_magnitudet[];
	foreach() {
		double temp = 0.0;
		foreach_dimension() {
			temp += sq(gradphit.x[]);
		}
		grad_magnitudet[] = sqrt(temp);
	}
	
	foreach() {
		foreach_dimension() {
			double keep=1.e-6;
			if (grad_magnitudet[] >= keep){
				nphit.x[] = gradphit.x[]/grad_magnitudet[];
		    }
			else{
				nphit.x[] = 0.0;  
		    }  
		}
	}

	boundary({nphit,gradphit});	

    scalar cut_type[];
	foreach() {
		cut_type[] = 1.0;
	}
	boundary({cut_type});	


	foreach() {
		int di, dj, dk;
		double phit_center = phit[] - 0.5;
		bool condition_met = false;
		
	
		#if dimension == 2
		// 2D check
		for (dj = -1; dj <= 1; dj++) {
			for (di = -1; di <= 1; di++) {
				// Skip the center point itself
				if (di == 0 && dj == 0) continue;
				// Check the condition
				if (phit_center * (phit[di, dj] - 0.5) <= 0.0) {
					condition_met = true;
					break;
				}
			}
			if (condition_met) break;
		}
		#endif


		#if dimension == 3
		// 3D check
		for (dk = -1; dk <= 1; dk++) {
			for (dj = -1; dj <= 1; dj++) {
				for (di = -1; di <= 1; di++) {
					// Skip the center point itself
					if (di == 0 && dj == 0 && dk == 0) continue;
					// Check the condition
					if (phit_center * (phit[di, dj, dk] - 0.5) <= 0.0) {
						condition_met = true;
						break;
					}
				}
				if (condition_met) break;
			}
			if (condition_met) break;
		}
		#endif
		
		if (condition_met) {
			cut_type[] = -1.0;
		}
		
	}




	foreach() {
		double T_interf=0., delta_T=0.;
		if (cut_type[] < 0.0 ){  
			double T_n=0., temp=0.;  
			foreach_dimension() {
				T_n  += gradT.x[]*nphit.x[];
				temp += sq(gradphit.x[]);
			}
			temp = sqrt(temp);
			T_interf = T[]+(0.5- phit[])/(temp +1.e-15) *T_n;

			delta_T = T_interf -Tm;
			Tm_eff[] = Tm_eff[] -delta_T;	
		}
		
	}


	boundary({Tm_eff,cut_type}); 

//*****************************************************************************************************************************************
   //Average Tm_eff;
	scalar Tm_eff2[];
	#if dimension == 2
	foreach() {

		double L[] = {HUGE, HUGE, HUGE, HUGE, HUGE, HUGE, HUGE, HUGE, HUGE};
		double w[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
		static int coords[9][2] = {
			{0, 0}, {0, 1}, {0, -1}, {1, 0}, {-1, 0}, {1, 1}, {-1, -1}, {-1, 1}, {1, -1}
		};
		for (int i = 0; i < 9; i++) {
			int ii = coords[i][0], jj = coords[i][1];
			if (cut_type[ii, jj] < 0.0) {
				double norm = sqrt(sq(gradphit.x[ii, jj]) + sq(gradphit.y[ii, jj])) + DBL_EPSILON;
				L[i] = fabs((phit[ii, jj] - 0.5) / norm);
			}
		}	

		
		if (cut_type[0, 0] < 0.0) {
			
			double W = 0.0;
			for (int i = 0; i < 9; i++) {
					w[i] = 1.0 / (L[i] + DBL_EPSILON);
					W += w[i];
			}
			for (int i = 0; i < 9; i++) {
				w[i] /= W;
			}

			
			Tm_eff2[] = w[0] * Tm_eff[]      + w[1] * Tm_eff[ 0, 1] + w[2] * Tm_eff[0, -1] + 
						w[3] * Tm_eff[ 1, 0] + w[4] * Tm_eff[-1, 0] + 
						w[5] * Tm_eff[ 1, 1] + w[6] * Tm_eff[-1, -1] + 
						w[7] * Tm_eff[-1, 1] + w[8] * Tm_eff[1, -1];
		} else {
			Tm_eff2[] = Tm_eff[];
		}
		
	}
	#endif	



	#if dimension == 3
	foreach() {
		double L[3][3][3], w[3][3][3];
		double W = 0.0;

		// Initialize the L and w arrays
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					L[i][j][k] = HUGE;
					w[i][j][k] = 0.0;
				}
			}
		}

		// Compute L-array
		for (int i = -1; i <= 1; i++) {
			for (int j = -1; j <= 1; j++) {
				for (int k = -1; k <= 1; k++) {
					if (cut_type[i, j, k] < 0.0 ) {
						double norm =(sqrt(sq(gradphit.x[i, j, k]) + sq(gradphit.y[i, j, k]) + sq(gradphit.z[i, j, k])) + DBL_EPSILON);
						L[i + 1][j + 1][k + 1] = fabs((phit[i, j, k] - 0.5) /norm );
					}
				}
			}
		}

		// Calculate the weights w and W
		if (cut_type[] < 0.0) {
			W = 0.0;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 3; k++) {
						W += 1.0 / (L[i][j][k] + DBL_EPSILON);
					}
				}
			}

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					for (int k = 0; k < 3; k++) {
						w[i][j][k] = 1.0 / (L[i][j][k] + DBL_EPSILON) / W;
					}
				}
			}

			// Initialize Tm_eff2 to zero before accumulation
			Tm_eff2[] = 0.0;			
			// Accumulate the weighted sum	
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					for (int k = -1; k <= 1; k++) {
						Tm_eff2[] += w[i + 1][j + 1][k + 1] * Tm_eff[i, j, k];
					}
				}
			}

		} else {
			Tm_eff2[] = Tm_eff[];
		}
	}
	#endif	

	boundary({Tm_eff2});
	foreach() {
		Tm_eff[] =Tm_eff2[];
	}
	boundary({Tm_eff}); 	

//%%%%Expand Tm_eff along the normal direction

	
	scalar RHS[];
	double dtau = 0.8 * dxmin;

	for (int kkk = 0; kkk < 1; kkk++) {
		// Advection for phit > 0.5
		foreach() {
			if (cut_type[] > 0.0 && (phit[] - 0.5) > 0) {
				coord psp={0,0,0};
				psp.x = (nphit.x[] < 0.0) ? (Tm_eff[1, 0] - Tm_eff[]) / Delta : (Tm_eff[] - Tm_eff[-1, 0]) / Delta;
				psp.y = (nphit.y[] < 0.0) ? (Tm_eff[0, 1] - Tm_eff[]) / Delta : (Tm_eff[] - Tm_eff[0, -1]) / Delta;
				#if dimension == 3
				psp.z = (nphit.z[] < 0.0) ? (Tm_eff[0, 0, 1] - Tm_eff[]) / Delta : (Tm_eff[] - Tm_eff[0, 0, -1]) / Delta;
				#endif
				RHS[] = 0.0;
				foreach_dimension() {
					RHS[] += nphit.x[] * psp.x;
				}
			}
		}

		foreach() {
			if (cut_type[] > 0.0 && (phit[] - 0.5) > 0.0) {
				Tm_eff[] -= dtau * RHS[];
			}
		}

		// Advection for phit < 0.5
		foreach() {
			if (cut_type[] > 0.0 && (phit[] - 0.5) < 0) {
				coord psp={0,0,0};
				psp.x = (-nphit.x[] < 0.0) ? (Tm_eff[1, 0] - Tm_eff[]) / Delta : (Tm_eff[] - Tm_eff[-1, 0]) / Delta;
				psp.y = (-nphit.y[] < 0.0) ? (Tm_eff[0, 1] - Tm_eff[]) / Delta : (Tm_eff[] - Tm_eff[0, -1]) / Delta;
				#if dimension == 3
				psp.z = (-nphit.z[] < 0.0) ? (Tm_eff[0, 0, 1] - Tm_eff[]) / Delta : (Tm_eff[] - Tm_eff[0, 0, -1]) / Delta;
				#endif
				
				RHS[] = 0.0;
				foreach_dimension() {
					RHS[] -= nphit.x[] * psp.x;
				}
			}
		}

		foreach() {
			if (cut_type[] > 0.0 && (phit[] - 0.5) < 0.0) {
				Tm_eff[] -= dtau * RHS[];
			}
		}
	}


	foreach() {
		if ((phit[]<0.01 || phit[]>0.99)){
			Tm_eff[] = Tm;
		}
	}
	boundary({Tm_eff});
	
	
}
//******************************-----------end Tm_eff_AC-------------------------------------------------------------------------------------------
///%Update  phi1----------------------------------------------------------------------------------------------------------------------------------------

	foreach() {
		double temp;

		temp = RHS_AC[]
			 +1.0/Cn/Cn*A_AC*(T[]-Tm_eff[])*phi1[]*(1.0-phi1[]);		 

		phi1new[] = phi1[]+ dt*D_AC*temp;
		phi1new[] = clamp (phi1new[], 0., 1.);
	}

	boundary({phi1new});
	foreach()
	{
		phi1_ns[] = (phi1[]+phi1new[])/2.;
		phi1old[] = phi1[];
		phi1[] = phi1new[];
	}
	boundary({phi1_ns,phi1,phi1old});

}
#endif

  // update rho & alpha
  event ("properties");


} //event AC_equation