 /* 
 * functions.h
 *
 * Purpose:
 *   - Read the lookup table of the Equation of State (EOS) from input files
 *   - Compute EOS values (density) by interpolation of the lookup table
 *   - Read and interpolate the background (ambient) stratification
 *   - Define pipe geometry/markers
 *
 */
 
 
 #include <sys/time.h>                                                                
 #include <time.h>                                                                    
 #include <stdlib.h>                                                                  
 #include <stdio.h> 
 #include <math.h>
 #include <string.h>


/* Interpolate background temperature at a given vertical position */	
double get_background_temperature( double y_ref_f[], double b_ref_f[], int taille_data_f, double y_f )
{
  double indice, delta_yref, alpha, b_interpolation_f;  
  delta_yref = y_ref_f[1] - y_ref_f[0];
  indice = ( y_f - y_ref_f[0]) / delta_yref ;
  int i_inter = floor(indice);
  if (i_inter < taille_data_f-1 && i_inter >= 0 )
  {
      alpha = (y_ref_f[i_inter+1] - y_f) / (y_ref_f[i_inter+1] - y_ref_f[i_inter]);
      b_interpolation_f = alpha * b_ref_f[i_inter] + (1.0 - alpha) * b_ref_f[i_inter + 1];
  }
  if (i_inter == taille_data_f-1)
  {
      b_interpolation_f = b_ref_f[i_inter];
   }
   
  if (i_inter < 0 )
  {
  	b_interpolation_f = b_ref_f[0];
  }
  
  	
  return b_interpolation_f;
}

/* Read stratification data (y, b) from file */
int def_data( double *b_ref_f, double *y_ref_f, int L0_f, double origine, char * name )
{
  FILE *myfile;
  double b_fichier, y_fichier;
  int taille_data,width_data;
  int i_f;
  int result;

  myfile=fopen(name, "r");
    
  fscanf(myfile,"%d %d",&taille_data,&width_data);

 
  for(i_f = 0 ; i_f < taille_data; i_f++)
  {
      	result = fscanf(myfile,"%lf %lf",&y_fichier,&b_fichier);
        if(result != EOF)
        {
            b_ref_f[i_f] = b_fichier;
            y_ref_f[i_f] = y_fichier;      
        }
  }
  fclose(myfile);

#if 0
 for(i_f = 0 ; i_f < taille_data; i_f++)
  {
    y_ref_f[i_f] =  i_f * (double) L0_f/ (double)(taille_data-1) + origine; 
    

  }
#endif

  return 0;

}

/* Read number of entries from file header */
int read_size_file(char * name){
	FILE *file;
	int taille_data;
	file = fopen(name, "r");
	fscanf(file,"%d",&taille_data);
	fclose(file);
	return taille_data;
}

/* Read EOS state table (P,T,rho) from file */
int read_state_table( double *P_state_eq, double *T_state_eq, double *rho_state_eq, double *Pmin, double *Pmax, double *delta_P, double *Tmin, double *Tmax, double *delta_T, char * name )
{
  FILE *myfile;
  double P_file, T_file, rho_file;
  int taille_data,width_data;
  int i_f, result;
  char bin[200] = {0};

  myfile=fopen(name, "r");
    
  fscanf(myfile,"%d %d",&taille_data,&width_data); 
  fscanf(myfile,"%lf %lf %lf",Pmin, Pmax, delta_P); fgets(bin, 200, myfile);
  fscanf(myfile,"%lf %lf %lf",Tmin, Tmax, delta_T); fgets(bin, 200, myfile);
  fgets(bin, 200, myfile);
 
  for(i_f = 0 ; i_f < taille_data; i_f++)
  {
      	result = fscanf(myfile,"%lf %lf %lf",&P_file,&T_file,&rho_file);
        if(result != EOF)
        {
            P_state_eq[i_f] = P_file;
            T_state_eq[i_f] = T_file;
            rho_state_eq[i_f] = rho_file;   
        }
  }
  fclose(myfile);

#if 0
  FILE *fp1;
  fp1 = fopen("state_equation_table.txt", "w+");  
 for(i_f = 0 ; i_f < taille_data; i_f++)
  {
    fprintf(fp1, "%.2lf %.2lf %.2lf \n",P_state_eq[i_f],T_state_eq[i_f],rho_state_eq[i_f]); 
  }
  fclose(fp1);
#endif
  return 0;

}

/* Bilinear interpolation of density from EOS state table */
trace
double interpolation_eq_table( double pr, double temp, double *P_state_eq, double *T_state_eq, double *rho_state_eq, double *Pmin, double *Pmax, double *delta_P, double *Tmin, double *Tmax, double *delta_T )
{
    double rho_P1_T, alpha_P1_T, rho_P2_T, alpha_P2_T, alpha_P, rho_P;
    
    int size_Temp = (*Tmax - *Tmin) / *delta_T + 1;
    int size_pr = (*Pmax - *Pmin) / *delta_P + 1;

    int indice_pr = ( pr - *Pmin) / *delta_P ;
    int indice_pr_low = floor(indice_pr);
    int indice_temp = ( temp - *Tmin) / *delta_T ;
    int indice_temp_low = floor(indice_temp);

    int indice_P1_T1 = indice_pr_low * size_Temp + indice_temp_low;
    int indice_P1_T2 = indice_P1_T1 + 1;

    int indice_P2_T1 = ( indice_pr_low + 1) * size_Temp + indice_temp_low;
    int indice_P2_T2 = indice_P2_T1 + 1;
    int loop = 1;
    
    if ( (indice_pr_low > size_pr -1)){
    printf("Error1_P : outside state equation range\n");
    indice_pr_low = size_pr -1;
    loop = 0;
    rho_P = 1000;
    }
    if ( (indice_temp_low > size_Temp -1)){
    printf("Error1_T : outside state equation range\n");
    indice_temp_low = size_Temp -1;
        loop = 0;
    rho_P = 1000;
    }
    if ( (indice_temp_low < 0) ){
    printf("Error2_T : outside state equation range iP %d iT: %d\n", indice_pr_low, indice_temp_low);
    indice_temp_low  = 0;
        loop = 0;
    rho_P = 1000;
    }   
    if ( (indice_pr_low < 0)){
    printf("Error2_P : outside state equation range iP %d iT: %d\n", indice_pr_low, indice_temp_low);
    indice_pr_low = 0;
        loop = 0;
    rho_P = 1000;
    } 

    if (loop == 1 ) {
    if (indice_temp_low == size_Temp -1) { rho_P1_T = rho_state_eq[indice_P1_T1];}
    else {
    alpha_P1_T = (T_state_eq[indice_P1_T2] - temp) / (T_state_eq[indice_P1_T2] - T_state_eq[indice_P1_T1]);
    rho_P1_T = alpha_P1_T * rho_state_eq[indice_P1_T1] + (1.0 - alpha_P1_T) * rho_state_eq[indice_P1_T2];
    }

    if (indice_temp_low == size_Temp -1) { rho_P2_T = rho_state_eq[indice_P2_T1];}
    else {
    alpha_P2_T = (T_state_eq[indice_P2_T2] - temp) / (T_state_eq[indice_P2_T2] - T_state_eq[indice_P2_T1]);
    rho_P2_T = alpha_P2_T * rho_state_eq[indice_P2_T1] + (1.0 - alpha_P2_T) * rho_state_eq[indice_P2_T2];
    }

    if (indice_pr_low == size_pr -1) {rho_P = rho_P1_T;}
    else {  
    alpha_P = (P_state_eq[indice_P1_T1 + size_Temp] - pr) / (P_state_eq[indice_P1_T1 + size_Temp] - P_state_eq[indice_P1_T1]);
    rho_P = alpha_P * rho_P1_T + (1.0 - alpha_P) * rho_P2_T;
    }
    }

    return  rho_P;
}

/* Compute absolute pressure (bar) from depth z */
double pressure_from_z( double z )
{
  double a = -1.0114775;
  double b = -0.85843076;
  double p = (a * z + b +  10.1325) / 10. ; //absolute pressure in bars
  return p;
}



	
/* Define outlet geometry using scalar volume fractions */
void define_pipe( scalar f0, scalar f1, scalar cs, double xcenter, double ycenter, double radius_exit, double pipe_height, double hplan, double pipe_width, int maxlevel, int minlevel, double p1, double p2){
  scalar f2[];
	
	int iteration = 0;
  do{
	fraction (f0, z<=pipe_height+p1*Delta+hplan? (sq(radius_exit) - sq(y-ycenter) - sq(x-xcenter)):0);
	}
  while (adapt_wavelet({f0}, (double []){1e-15,1e-15}, maxlevel = maxlevel, maxlevel).nf &&
	 iteration++ <= 10);
  
  foreach(){
    if ( z <= (Z0 + hplan + p2*Delta))
      f2[] = 1.;
    else
      f2[] = 0.;
  }
  

	
	f0.refine = f0.prolongation = fraction_refine;
  f2.refine = f2.prolongation = fraction_refine;
  restriction ({f0,f2}); // for boundary conditions on levels	
  boundary ({f0,f2});

//somehow some pieces are missing
	int icheck,jcheck;
	foreach(){
		if( f0[] != f0[0,0,-1]){
			for(icheck = -1; icheck <= 1; icheck++){
				for(jcheck = -1 ; jcheck <= 1; jcheck++){
					if (f0[jcheck,icheck,0] > 1e-8){
						f0[] = 1.;
					}
				}
			}
		}
	}
	boundary ({f0});
			
	foreach(){
		f1[] = 0.;
		if(f0[]>1e-8){
			f0[] = 1.;
			f1[] = 1.;
		}
		else{ f0[] = 0.;}
	}
	 boundary ({f0,f1});
	int i,j;
	foreach(){
		for(i = -pipe_width ; i <= pipe_width; i++){
			for(j = -pipe_width ; j <= pipe_width; j++){
				if(f0[] < 1e-8){
					if (f0[i,j,0] > 1e-8 ){ //careful, f0 should be far from boundary (core dumped error)
						f1[] = 1.;
					}
				}
			}
		}
	}
	 boundary ({f1});
	foreach(){
	f1[] = (f1[]*(1-f0[])) ;
	}
	f1.refine = f1.prolongation = fraction_refine;
	restriction ({f1}); 
	boundary ({f1,f0});
  
#if EMBED	
	vertex scalar phi[];
	foreach_vertex()
		phi[] =  15; //rajour for reset phi in restore=> doesn't work with HUGE in mpi !!!!
	boundary ({phi});  
	fractions (phi, cs, fs);
	fractions_cleanup (cs, fs);
  boundary({cs,fs});				
	
	foreach_vertex(){
		phi[] =  HUGE;
		if (f1[] > 1e-8){
			phi[] = 0;
			phi[0,-1,0] = 0; //cs calculated based on vertex, this is needed to have the full pipe
		}
	}
  foreach_vertex(){
    if(f2[]>1e-8 && f0[] < 1e-8 ){
      phi[] = 0;

    }
  }
	boundary ({phi});  
	fractions (phi, cs, fs);
	fractions_cleanup (cs, fs);
	cs.refine=cs.prolongation=fraction_refine;
  boundary({cs,fs});
  restriction({cs,fs});
#else
  printf("error on f0, need embed");
#endif
  foreach(){
    f1[0,-1,0] = f1[];
    f0[0,-1,0] = f0[];
  }
	boundary ({f1,f0});
  foreach(){
    if(f0[] != f0[0,0,1]){
      f0[] = 0.;
    }
  }
  boundary ({f1,f0});
}	

/* Define marker fields at the outlet*/
void define_marker(scalar f0, scalar marker_analysis, scalar marker_injection, scalar marker_analysis_minus, scalar marker_injection_minus, scalar marker_ref){
	scalar * list_marker = {marker_analysis,marker_injection,marker_analysis_minus,marker_injection_minus,marker_ref};
	int i_marker,j_marker,k_marker;
	
	boundary({f0});	
	foreach(){
		for (scalar m in list_marker)
			m[] = 0.;
		
		for(k_marker = 0 ; k_marker >= -2; k_marker--){
				if(f0[] != f0[0,0,k_marker] &&  f0[0,0,k_marker] == 1. )
					marker_ref[] = 1.;							
				else
					marker_ref[] = 0.;
		}
	}
	for (scalar m in list_marker)
		boundary({m});
	
	foreach(){
		for(i_marker = -1; i_marker <= 1; i_marker++){
			for(j_marker = -1 ; j_marker <= 1; j_marker++){
				if( marker_analysis[] == 0. && marker_ref[i_marker,j_marker,0] == 1.){
					marker_analysis_minus[] = 1.;
					marker_analysis[0,0,1] = 1.;
				}				
			}
		}
		if( marker_injection[] != 1. && marker_ref[] == 1.){
			marker_injection_minus[] = 1.;
			marker_injection[0,0,1] = 1.;
		}
	}
	for (scalar m in list_marker){
		m.refine = m.prolongation = fraction_refine;
		restriction ({m}); 	
		boundary({m});
	}
}








