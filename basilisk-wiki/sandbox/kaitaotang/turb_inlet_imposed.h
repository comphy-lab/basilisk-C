/**
# Generation of artificial turbulent inlet conditions
See [Xie and Castro (2008)](https://doi.org/10.1007/s10494-008-9151-5) for details regarding the synthetic digital filtering method*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"

#include <sys/stat.h>
#include <sys/types.h>

int set_inlet_count = 0;
int n_yz, ngrid_yz, Nf_yz;
double stp; double t_pre = 0.;
#define N_max 1600

double Rtilde[3][N_max][N_max];
double psi_turb[3][N_max][N_max]; double psi_turb_old[3][N_max][N_max];
double ui[3][N_max][N_max];
scalar ux[]; scalar uy[]; scalar uz[];



// Take a 2D slice of the flow field at x = xp
void sliceYZ(double ** field, scalar s, double xp, int maxlevel){
  int nn = (1 << maxlevel);
  double stp = L0/(double)(nn);
  
  for (int i = 0; i < nn; i++) {
      double yp = stp*i + Y0 + stp/2.;
      for (int j = 0; j < nn; j++) {
	  double zp = stp*j + Z0 + stp/2.;
	  Point point = locate (xp, yp, zp);
	  field[i][j] = point.level >= 0 ? s[] : nodata;
      }
  }
  
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, field[0], sq(nn+1), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
}



/** 
Take a 2D slice of the flow field s at y = yp */
void sliceXZ(double ** field, scalar s, double yp, int maxlevel){
  int nn = (1 << maxlevel);
  double stp = L0/(double)(nn);
  
  for (int i = 0; i < nn; i++) {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) {
	  double zp = stp*j + Z0 + stp/2.;
	  Point point = locate (xp, yp, zp);
	  field[i][j] = point.level >= 0 ? s[] : nodata;
      }
  }
  
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, field[0], sq(nn+1), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
}



/** 
Take a 2D slice of the flow field s at z = zp */
void sliceXY(double ** field, scalar s, double zp, int maxlevel){
  int nn = (1 << maxlevel);
  double stp = L0/(double)(nn);
  
  for (int i = 0; i < nn; i++) {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) {
	  double yp = stp*j + Y0 + stp/2.;
	  Point point = locate (xp, yp, zp);
	  field[i][j] = point.level >= 0 ? s[] : nodata;
      }
  }
  
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, field[0], sq(nn+1), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
}


/** 
Auxiliary function: normalise a matrix to zero mean and unit variance */
void norm_matrix(double mat[3][N_max][N_max], int num_xy) {
  double r_ave[3]; double r_sq_ave[3]; double r_var[3];
  
  for (int i=0; i<3; i++) {
    r_ave[i] = 0.; r_sq_ave[i] = 0.;

    for (int j=0; j<num_xy; j++)
      for (int k=0; k<num_xy; k++) {
        r_ave[i] += mat[i][j][k];
        r_sq_ave[i] += sq(mat[i][j][k]);
      }
    
    r_ave[i] /= sq(num_xy); r_sq_ave[i] /= sq(num_xy);
    r_var[i] = r_sq_ave[i] - sq(r_ave[i]);    
  }      

  for (int i=0; i<3; i++)
    for (int j=0; j<num_xy; j++)
      for (int k=0; k<num_xy; k++) mat[i][j][k] = (mat[i][j][k] - r_ave[i]) / sqrt(r_var[i]);
}



void set_inlet(int lev, double L_turb, double u_ave, double u_rms) {
  double b_yz[N_max];
  
  ngrid_yz = 1 << lev; stp = L0/(double)(ngrid_yz);
  
  // Filter length setup
  n_yz = (int)(L_turb/stp); Nf_yz = 3 * n_yz; 
  
  
  
  if (set_inlet_count > 0) {
    for (int i=0; i<3; i++)
      for (int j=0; j<ngrid_yz; j++)
        for (int k=0; k<ngrid_yz; k++) psi_turb_old[i][j][k] = psi_turb[i][j][k];
  }
   
  // Set up random fields and the Gaussian filter 
  if (pid() == 0) {
    srand(time(NULL));
     
    for (int i=0; i<3; i++) 
      for (int j=0; j<2*Nf_yz+ngrid_yz; j++)
        for (int k=0; k<2*Nf_yz+ngrid_yz; k++) Rtilde[i][j][k] = 2*(rand()/(double)RAND_MAX)-1;
  }  
  
  MPI_Bcast (&(Rtilde[0][0][0]), 3*sq(N_max), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  norm_matrix(Rtilde, 2*Nf_yz+ngrid_yz);
    
    
    
  double norm_byz = 0.;    
  for (int i=0; i<2*Nf_yz; i++) {b_yz[i] = exp(-pi*fabs(i-Nf_yz)/2.0/n_yz); norm_byz += sq(b_yz[i]);}
  for (int i=0; i<2*Nf_yz; i++) b_yz[i] /= sqrt(norm_byz);
  
  // Apply filtering
  for (int i=0; i<3; i++) 
    for (int j=0; j<ngrid_yz; j++)
      for (int k=0; k<ngrid_yz; k++) {
        psi_turb[i][j][k] = 0.;
          
        if ((k + j*(ngrid_yz)) % npe() == pid()) {            
          for (int j_p=0; j_p<2*Nf_yz; j_p++)
            for (int k_p=0; k_p<2*Nf_yz; k_p++) psi_turb[i][j][k] += b_yz[j_p]*b_yz[k_p] * Rtilde[i][j+j_p][k+k_p];                
        }
      }
    
  MPI_Allreduce (MPI_IN_PLACE, &(psi_turb[0][0][0]), 3*sq(N_max), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);        
  norm_matrix(psi_turb, ngrid_yz);
    
      
  
  // Time advancing
  if (set_inlet_count > 0) {    
    double T_La = L_turb/u_ave;
  
    for (int i=0; i<3; i++)
      for (int j=0; j<ngrid_yz; j++)
        for (int k=0; k<ngrid_yz; k++)
          psi_turb[i][j][k] = psi_turb_old[i][j][k] * exp(-pi*(t-t_pre)/4.0/T_La) + psi_turb[i][j][k] * pow(1 - exp(-pi*(t-t_pre)/2.0/T_La), 0.5); 
      
    norm_matrix(psi_turb, ngrid_yz);
  }

  // Compute the full velocity profile
  for (int i=0; i<3; i++)
    for (int j=0; j<ngrid_yz; j++)
      for (int k=0; k<ngrid_yz; k++) ui[i][j][k] = psi_turb[i][j][k]*u_rms;
  


  // Inlet flux correction
  double u_in_ave[3] = {0., 0., 0.};
  for (int i=0; i<3; i++)
    for (int j=0; j<ngrid_yz; j++)
      for (int k=0; k<ngrid_yz; k++) u_in_ave[i] += ui[i][j][k] / sq(ngrid_yz);
  
  for (int i=0; i<3; i++)
    for (int j=0; j<ngrid_yz; j++)
      for (int k=0; k<ngrid_yz; k++) {
        ui[i][j][k] -= u_in_ave[i];
        if (i == 0) ui[i][j][k] += u_ave;
      }
       
         
         
  // Interpolate to the inlet
  foreach() ux[] = uy[] = uz[] = 0.; 

  foreach_boundary(left) {
    int i_tmp = (int)((y - stp/2. - Y0)/stp);
    int j_tmp = (int)((z - stp/2. - Z0)/stp);
    
    if (i_tmp < 0) i_tmp = 0; if (i_tmp > ngrid_yz-1) i_tmp = ngrid_yz-1;
    if (j_tmp < 0) j_tmp = 0; if (j_tmp > ngrid_yz-1) j_tmp = ngrid_yz-1;
    
    ux[] = ui[0][i_tmp][j_tmp]; 
    uy[] = ui[1][i_tmp][j_tmp]; 
    uz[] = ui[2][i_tmp][j_tmp];
  }
  
  u.n[left]  = dirichlet(ux[]);
  u.t[left]  = dirichlet(uy[]);
  u.r[left]  = dirichlet(uz[]);
  
  boundary ((scalar *){u});
  
  t_pre = t;
  
  set_inlet_count++;
  if (set_inlet_count == INT_MAX) set_inlet_count = 1; 
}



void match_inlet(int lev, double WLayer, double u_ave_in) {
  int nn = (1 << lev);
  double stp = L0/(double)(nn);
  
  double ** ux_samp = matrix_new (nn+1, nn+1, sizeof(double));
  double ** uy_samp = matrix_new (nn+1, nn+1, sizeof(double));
  double ** uz_samp = matrix_new (nn+1, nn+1, sizeof(double));
  
  sliceYZ(ux_samp, u.x, X0+WLayer, lev); sliceYZ(uy_samp, u.y, X0+WLayer, lev); sliceYZ(uz_samp, u.z, X0+WLayer, lev);
  
  
  
  double ux_samp_ave = 0.0; int num_samp_c = 0;
  
  foreach(reduction(+:ux_samp_ave) reduction(+:num_samp_c)) {
    if (X0 + WLayer - Delta/2.0 <= x && x <= X0 + WLayer + Delta/2.0) {
      ux_samp_ave += u.x[]; num_samp_c++;
    }
  }
  
  ux_samp_ave /= num_samp_c;
  
  
  
  double ux_intp_ave = 0.0; double s_in = 0;
  
  foreach_boundary(left, reduction(+:ux_intp_ave) reduction(+:s_in)) {
    int i_tmp = (int)((y - stp/2. - Y0)/stp);
    int j_tmp = (int)((z - stp/2. - Z0)/stp);
    
    if (i_tmp < 0) i_tmp = 0; if (i_tmp > nn) i_tmp = nn;
    if (j_tmp < 0) j_tmp = 0; if (j_tmp > nn) j_tmp = nn;
    
    ux[] = ux_samp[i_tmp][j_tmp]; ux_intp_ave += ux[]*sq(Delta); s_in += sq(Delta);
    uy[] = uy_samp[i_tmp][j_tmp];
    uz[] = uz_samp[i_tmp][j_tmp]; 
  } 
  
  ux_intp_ave = ux_intp_ave / s_in;
  
  foreach_boundary(left) ux[] = ux[] - ux_intp_ave + u_ave_in;
  
  u.n[left]  = dirichlet(ux[]);
  u.t[left]  = dirichlet(uy[]);
  u.r[left]  = dirichlet(uz[]);
  
  matrix_free (ux_samp); matrix_free (uy_samp); matrix_free (uz_samp);
}