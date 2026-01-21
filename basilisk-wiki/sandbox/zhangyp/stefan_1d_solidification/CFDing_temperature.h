event Temperature (i++) {
#ifdef TmpE_on
int marcherTE = 1;
//1 means FT-CS
//2 means AB-CS
//3 means CNAB-CS
scalar Cp_TE[],k_TE[],rho_TE[];
scalar Tnew[];

foreach()
{
		double phi_1 = phi1[];
		double phi_2 = 1.0-phi1[];

		Cp_TE[] =  phi_1 + phi_2*rCp;
		rho_TE[]=  phi_1 + phi_2*rd;
		k_TE[]  =  phi_1 + phi_2*rk;
}



face vector temp_la[];
foreach_face(){
	temp_la.x[] = face_gradient_x(T,0)*(k_TE[]+k_TE[-1])/2.0;
}

scalar source_TE[];
foreach()
{   Tnew[]=T[];
	if(i == 0)
		source_TE[] =0.;
	else 
		source_TE[] = -0.5 * convprev_TE[];  
}


foreach(){
	double la_f=0.0, temp=0.0, temprhs=0.0;  	
	foreach_dimension()
		la_f +=  (fm.x[1]*temp_la.x[1]-fm.x[]*temp_la.x[])/(Delta*max(cm[],1.e-20));
	if(marcherTE == 1 || marcherTE == 2 ){  
		temprhs = 1.0/rho_TE[]/Cp_TE[]*
				  (la_f/Pe - St*rd*(phi1[]-phi1old[])/dt);        
	}
	else{
		temprhs = 1.0/rho_TE[]/Cp_TE[]*
				  (la_f/Pe/2.0 - St*rd*(phi1[]-phi1old[])/dt);
	}
	convprev_TE[] = temprhs;
	
	if( i == 0)
	   source_TE[] += 1. * temprhs;
    else source_TE[] += 1.5 * temprhs;
	
}

  if(marcherTE == 1){ 
    foreach()
        Tnew[] = T[] + dt*convprev_TE[];
	}
   else if (marcherTE == 2){
    foreach()
        Tnew[] = T[] + dt*source_TE[];
	}
  else if (marcherTE == 3){
   foreach()
	{ 
		Tnew[] = T[];
	}
	for(int ii = 0; ii< 30 ; ii++)
    { 

	   foreach_face(){
			temp_la.x[] = face_gradient_x(Tnew,0)*(k_TE[]+k_TE[-1])/2.0;
	   }
	   	boundary({temp_la});
		foreach() {
			double  a2 = 4. / (sq(Delta));
			double  diagonal = 1. + 0.5*dt*a2*k_TE[]/rho_TE[]/Cp_TE[]/Pe;
			double  la_f = 0.0;
			foreach_dimension() {
				la_f +=  (fm.x[1]*temp_la.x[1]-fm.x[]*temp_la.x[])/(Delta*max(cm[],1.e-20));
			}
			la_f = la_f/rho_TE[]/Cp_TE[]/Pe;
			double  some = dt * (source_TE[] + 0.5*la_f) + T[] -Tnew[];
			Tnew[] = Tnew[] + 0.8 * some /diagonal;	
	    }
	  boundary({Tnew});

    }

  }
  else{ 
    printf("Iteration error in advance temperature equation without Phase transformation.");
  }

	foreach()
    {
		T[] = Tnew[]; 
    }
	boundary({T});
	
#endif 
}// enf of 	event Temperature (i++) 
​