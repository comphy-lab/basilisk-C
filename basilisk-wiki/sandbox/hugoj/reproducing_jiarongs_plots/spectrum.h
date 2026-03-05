/**
# PM spectrum read-in and wave field initialization functions

The initial condition is a random broad-banded wave field based on externally imported power spectrum. Some spectrum related variables: */

#define N_mode_ 32 // Corresponds to input of 32 modes in kx and 33 modes in ky. This number has to match the files being read in
double F_kxky_[N_mode_*(N_mode_+1)], omega[N_mode_*(N_mode_+1)], \
  phase[N_mode_*(N_mode_+1)];
double kx_[N_mode_], ky_[N_mode_+1];
double dkx_, dky_;
int RANDOM; // integer to seed random number generator

/**
   Random number generator. */
double randInRange(int min, int max)
{
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}

/** 
A MPI compatible function that reads in kx_, ky_, and F_kxky_. Next step is to generate F_kxky_ inside basilisk too. For now kx_, ky_ are 1D arrays while F_kxky_ is a 2D array. The root process reads in and then broadcast to all other processes. It looks for files named F_kxky, 
*/

void power_input() {
// Previously we were not reading in kx_, ky_.
//  for (int i=0; i<N_mode_; i++) { 
//     kx_[i] = 2.*pi/L0*(i+1); 
//     ky_[i] = 2.*pi/L0*(i-N_mode_/2); 
//  } 
//  ky_[N_mode_] = 2.*pi/L0*N_mode_/2; 

/* If using MPI. **/
#if _MPI 
  int length1D, length2D; // The length of the array to be read
  int rank, size;
  int root = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Read in on the zeroth thread
  if (rank == root) {
    // First read in F_kxky_
    length2D = N_mode_*(N_mode_+1);
    float * a = (float*) malloc (sizeof(float)*length2D);
    char filename[100];
    sprintf (filename, "F_kxky");
    FILE * fp = fopen (filename, "rb");
    fread (a, sizeof(float), length2D, fp);
    for (int i=0;i<length2D;i++) {
      F_kxky_[i] = (double)a[i];
    }
    fclose (fp);

    // Then read in kx_, ky_
    length1D = N_mode_;
    float * b1 = (float*) malloc (sizeof(float)*length1D);
    sprintf (filename, "kx");
    FILE *fp1 = fopen (filename, "rb");
    fread (b1, sizeof(float), length1D, fp1);
    for (int i=0;i<length1D;i++) {
      kx_[i] = (double)b1[i];
    }
    fclose (fp1);

    // One more mode in ky
    float * b2 = (float*) malloc (sizeof(float)*(length1D+1));
    sprintf (filename, "ky");
    FILE *fp2 = fopen (filename, "rb");
    fread (b2, sizeof(float), length1D+1, fp2);
    for (int i=0;i<length1D+1;i++) {
      ky_[i] = (double)b2[i];
    }
    fclose (fp2);

    // Wave frequency omega, and randomly generated phase
    double kmod = 0;
    int index = 0;
    srand(RANDOM); // We can seed it differently for different runs
    for (int i=0; i<N_mode_; i++) {
      for (int j=0; j<N_mode_+1; j++) {
	index = j*N_mode_ + i;
	kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
	omega[index] = sqrt(g_*kmod);
	phase[index] = randInRange (0, 2.*pi);
      }
    }
  }
  // Broadcast to other threads
  MPI_Bcast(&kx_, length1D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&ky_, length1D+1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&F_kxky_, length2D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&omega, length2D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&phase, length2D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  // Make sure that the inputs are correct by printing them out
  /* char checkout[100]; */
  /* sprintf (checkout, "F-%d", pid()); */
  /* FILE * fout = fopen (checkout, "w"); */
  /* for (int i=0; i<length2D; i++) */
  /*   fprintf (fout, "%g ", F_kxky_[i]); */
  /* fclose (fout); */
  /* sprintf (checkout, "ky-%d", pid()); */
  /* fout = fopen (checkout, "w"); */
  /* for (int i=0; i<length1D+1; i++) */
  /*   fprintf (fout, "%g ", ky_[i]); */
  /* fclose (fout); */

/**
   If not using MPI. */ 
#else
  int length1D, length2D;
  length2D = N_mode_*(N_mode_+1);
  float * a = (float*) malloc (sizeof(float)*length2D);
  char filename[100];
  sprintf (filename, "F_kxky");
  FILE * fp = fopen (filename, "rb");
  fread (a, sizeof(float), length2D, fp);
  for (int i=0;i<length2D;i++) {
    F_kxky_[i] = (double)a[i];
  }
  fclose (fp);

  // Then read in kx_, ky_
  length1D = N_mode_;
  float * b1 = (float*) malloc (sizeof(float)*length1D);
  sprintf (filename, "kx");
  FILE *fp1 = fopen (filename, "rb");
  fread (b1, sizeof(float), length1D, fp1);
  for (int i=0;i<length1D;i++) {
    kx_[i] = (double)b1[i];
  }
  fclose (fp1);

  // One more mode in ky
  float * b2 = (float*) malloc (sizeof(float)*(length1D+1));
  sprintf (filename, "ky");
  FILE *fp2 = fopen (filename, "rb");
  fread (b2, sizeof(float), length1D+1, fp2);
  for (int i=0;i<length1D+1;i++) {
    ky_[i] = (double)b2[i];
  }
  fclose (fp2);

  // Phase and omega, next focusing phase
  double kmod = 0;
  int index = 0;
  srand(0); 
  for (int i=0; i<N_mode_;i++) {
    for (int j=0;j<N_mode_+1;j++) {
      index = j*N_mode_ + i;
      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      omega[index] = sqrt(g_*kmod);
      phase[index] = randInRange (0, 2.*pi);
    }
  }
#endif
}

/**
   Functions that compute the orbital velocity (based on linear wave equations). */
double wave (double x, double y)
{
  double eta = 0;
  double ampl = 0, a = 0;
  int index = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      index = j*N_mode_ + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      eta += ampl*cos(a);
    }
  }
  return eta;
}
double u_x (double x, double y, double z) {
  int index = 0;
  double u_x = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0, theta = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      index = j*N_mode_ + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      z_actual = (z < ampl ? (z) : ampl);
      // fprintf(stderr, "z = %g, ampl = %g, z_actual = %g\n", z, ampl, z_actual);
      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      theta = atan(ky_[j]/kx_[i]);
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      u_x += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*cos(a)*cos(theta);
    }
  }
  return u_x;
}

double u_y (double x, double y, double z) {
  int index = 0;
  double u_y = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0, theta = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      index = j*N_mode_ + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      z_actual = (z < ampl ? (z) : ampl);
      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      theta = atan(ky_[j]/kx_[i]);
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      u_y += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*cos(a)*sin(theta);
    }
  }
  return u_y;
}

double u_z (double x, double y, double z) {
  int index = 0;
  double u_z = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0;
  for (int i=0; i<N_mode_; i++) {
    for (int j=0; j<N_mode_+1; j++) {
      index = j*N_mode_ + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      z_actual = (z < ampl ? (z) : ampl);
      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      u_z += sqrt(g_*kmod)*ampl*exp(kmod*z_actual)*sin(a);
    }
  }
  return u_z;
}
