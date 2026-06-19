#include "grid/octree.h"               // Octree mesh
#include "embed.h"
#include "navier-stokes/centered.h"    // Navier-Stokes solver
#include "fractions.h"
#include "SGS.h"                    // Subgrid-scale turbulence model
#include "output_vtu_foreach.h"     // VTU output
#include "functions.h"        // Custom helper functions

//----------------- Simulation parameters ----------------//
int end_simu = 45000;
double lenght = 6;	
int maxlevel = 8;                     
int taille_grid = pow(2,8); 
int minlevel = 5;
double param_turb = 0.5;               // Turbulence forcing factor
double eps = 1e-8;
double timestep_ini = 10;

// Post-process
int paraview_record_step = 20;

// Input: jet/plume source
/*
To reproduce the different experiments presented in the manuscript,
users should change the input field, as detailed in Table 1 of the paper.
"Hydrothermal Plume Near-Field Dynamics from LES
and Observations"
*/
double T_input = 300.;
double u_exit = 0.7;
double radius_exit = 0.0076;
double pipe_height = 0., hplan;

// LES constants
double define_molvis = 1.58e-6; 
double molvis_paroi = 1.58e-6;
double define_csmag = 0.1;

// Stratification and EOS
double depth = -1690; 
double rho_offset = 1035.27; 
double rho_ref = 1035; 
double gravity = 9.803; 

// Equation of state table
#define taille_state_eq_cst 1403901
#define taille_data_T 3339
char * file_eq_state = "./input/equation_of_state";
char * file_T_data_stratification = "./input/stratification";

//----------------- Fields ----------------//
scalar b[], T_insitu[], * tracers = {T_insitu};
scalar T_strat[], rho_for_b[];
face vector av[];
const face vector muc[];

scalar absolute_pressure2[], z_rescale1[], f0[], f1[], rho_strat[];

scalar marker_analysis[], marker_injection[], marker_analysis_minus[], marker_injection_minus[], marker_ref[];

//----------------- State equation arrays ----------------//
double P_state_eq[taille_state_eq_cst], T_state_eq[taille_state_eq_cst], rho_state_eq[taille_state_eq_cst];
int taille_state_eq;
double T_momar[taille_data_T], z_momar[taille_data_T];
double Pmin, Pmax, delta_P, Tmin, Tmax, delta_T;
double z_rescale, absolute_pressure, temp_eq;


//----------------- Boundary conditions ----------------//
Evis[back]   = dirichlet (molvis_paroi);
Evis[embed]  = neumann (0.);
Evis[front]  = neumann (0.);
Evis[left]   = neumann (0.);
Evis[right]  = neumann (0.);
Evis[bottom] = neumann (0.);
Evis[top]    = neumann (0.);

// Injection BC
u.n[back]  = f0[] > eps ? dirichlet(u_exit) : dirichlet(0.);
T_insitu[back]  = f0[] > eps ? dirichlet(T_input) : dirichlet(4.6);

//Embed BC
u.n[embed] = f1[]>eps || f0[]>eps ? dirichlet(0.) : neumann(0.);
u.t[embed] = f1[]>eps || f0[]>eps ? dirichlet(0.) : neumann(0.);
u.r[embed]      = dirichlet (0.);
T_insitu[embed] = neumann(0.);

// Outlet BC
T_insitu[front] = neumann(0.);
uf.n[front]     = neumann(0.); 

// Horizontal BC
u.n[left]   = neumann(0.);
T_insitu[left]  = neumann(0.);
u.n[right]  = neumann(0.);
T_insitu[right] = neumann(0.);
u.n[bottom] = neumann(0.);
T_insitu[bottom]= neumann(0.);
u.n[top]    = neumann(0.);
T_insitu[top]   = neumann(0.);

// Pressure BC
p[front]  = dirichlet(0.);		
pf[front] = dirichlet(0.);
p[embed]  = neumann(0.);

//----------------- Main ----------------//
int main() {
  a = av;
	
  X0 = Y0 = Z0 = 0.;
	
  foreach_dimension()
    u.x.refine = refine_linear; 
		
  p.refine = p.prolongation = refine_linear;    
	
  T_insitu.gradient = minmod2; 
	
  for (scalar s in {T_insitu, u, g})
    s.prolongation = refine_embed_linear;	
	
  run();
}

//----------------- INIT ----------------//
event init (t = 0) {
  FILE * fp = fopen("count_cell", "w");
  fprintf (fp,"n_iteration\tn_cell\tdt\ttime\twall-time\tcheck-noise\tnb_parts\n");
  fflush (fp);
	fclose (fp);
	
	// Read EOS table
	/*
		The EOS used here is the parametric seawater equation of state 
		developed by Lemaréchal et al. (2024), designed for hydrothermal 
		vent conditions.
		Valid for: T = 0350 °C, P = 100550 bar, S_A = 35.2 g/kg.
		source :
			Lemaréchal, C., Roullet, G., & Gula, J. (2024).
			Parametric Equation of State of Seawater for
			Hydrothermal Vent Conditions.
			Zenodo. https://doi.org/10.5281/zenodo.14332032
	*/
  taille_state_eq = read_size_file(file_eq_state);
  read_state_table(P_state_eq, T_state_eq, rho_state_eq,
                   &Pmin, &Pmax, &delta_P, &Tmin, &Tmax, &delta_T,
                   file_eq_state);
  
  // Stratification profile
  def_data(T_momar, z_momar, L0, Z0, file_T_data_stratification);
  
  molvis = define_molvis;
  Csmag  = define_csmag;
  hplan = lenght/pow(2,4);
	 
	init_grid (taille_grid); 
	size (lenght + hplan );

	// Define pipe and markers
	define_pipe(f0,f1,cs,X0+L0/2.+L0/pow(2,maxlevel),Y0+L0/2.,
							radius_exit,pipe_height,hplan,2,maxlevel,minlevel,2,1);
	cs.refine=cs.prolongation=fraction_refine;
	restriction({cs,fs});		
	define_marker(f0,marker_analysis,marker_injection,
								marker_analysis_minus,marker_injection_minus,marker_ref);

	// Initial stratification
  foreach() {
		z_rescale = z + depth - pipe_height - hplan; 
		T_insitu[] = get_background_temperature(z_momar, T_momar, taille_data_T, z_rescale);	 
		temp_eq = T_insitu[];		
		absolute_pressure = pressure_from_z( z_rescale );
		rho_for_b[] = interpolation_eq_table(absolute_pressure,temp_eq,P_state_eq, T_state_eq, rho_state_eq, &Pmin, &Pmax, &delta_P, &Tmin, &Tmax, &delta_T);
		b[] = -gravity/rho_ref*(rho_for_b[] - rho_offset);
		rho_strat[] = rho_for_b[] ;
		T_strat[] = T_insitu[];		
		z_rescale1[] = z_rescale;
		absolute_pressure2[] = absolute_pressure;
	}
	
	foreach(){
		if (f0[] > eps ){
			T_insitu[] = T_input;
			u.z[] = u_exit;
		}
	}
	boundary({u,b,rho_strat,rho_for_b,T_strat,T_insitu, z_rescale1,absolute_pressure2});
						
	foreach() {
		if (cs[] < eps){
			T_insitu[] = 4.;
			foreach_dimension()
				u.x[] = 0.;
		}
	}
	boundary({T_insitu});
	boundary ((scalar*){u});
	
  DT = timestep_ini;
}

//----------------- ADD BUOYANCY FORCE & COMPUTE EOS ----------------//
event acceleration (i++) {
  boundary({T_insitu});	
  foreach() {
		z_rescale =  z + depth - pipe_height - hplan;
		absolute_pressure = pressure_from_z( z_rescale );
		temp_eq = T_insitu[];
		rho_for_b[] = interpolation_eq_table(absolute_pressure,temp_eq,P_state_eq, T_state_eq, rho_state_eq, &Pmin, &Pmax, &delta_P, &Tmin, &Tmax, &delta_T);
		b[] = -gravity/rho_ref*(rho_for_b[] - rho_offset);
 }
  boundary({b}); 	
  
  foreach_face(z){
    av.z[] = (b[] + b[0,0,-1])/2.;
  }

  // Add noise at injection to trigger turbulence
  boundary({u,marker_injection});
	
  foreach(){
    if(marker_injection[] > 0.5){
      u.z[] = u_exit * (1+noise()*param_turb);
    }
  }

  boundary({u,av,rho_for_b,T_insitu});
}

//----------------- CLIP ----------------//
event limiter (i++){
  foreach_boundary(front){
    if( T_insitu[] < T_strat[]){
      T_insitu[] = T_strat[]; //prevent spurious cooling due to numerical errors
    }
  }
}

//----------------- ADAPT ----------------//
event adapt (i++) {
  double uemax = 0.1;
  adapt_wavelet ({b,u,cs},(double[]){0.01,1e-3,0.01,uemax}, maxlevel, minlevel);
  unrefine ( f0[] < eps && (cs[] < eps || z <= hplan));
}

//------------------MONITOR THE RUN--------//
event count_cell(i++){
	timing count_cell_time = timer_timing (perf.gt, iter, perf.tnc, NULL);
	double delta_countcell = count_cell_time.cpu;
	timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
	if (pid() == 0){printf("starting time %g\n", delta_countcell);}
	int n = 0;
	foreach(reduction(+:n)){n++;}

	if (pid() == 0){
		FILE * fp = fopen("count_cell", "a");	
		fprintf (fp,"%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\n", i, n,dt,t,s.real, s.cpu, 0., 0);
		printf("cells : %d\n",n);
		fflush (fp);
		fclose (fp); 
	}
	boundary (all);

}

//----------------- EXPORT ----------------//
event export_to_paraview (i++){
  boundary (all);
  if (i%paraview_record_step == 0 || i==0){
    char name[350];
    sprintf(name, "wip%d",i);

    boundary ((scalar*){u});

		output_vtu((scalar *) {av, u,b,cs,T_insitu,T_strat,rho_strat,p,Evis,absolute_pressure2,z_rescale1,rho_for_b,mu,f0,f1}, (vector *) {0}, name);

	}
}

//----------------- END ----------------//
event end (i= end_simu) {
}