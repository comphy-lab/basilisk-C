/**
# Setting dimension of array raise error at dimension check

T_profile should inherit the dimension of T. */

#define LAYERS 1

double * T_profile;
double T0 = 20. [0,0,1];
scalar T;

int main(){
  nl = 2;
  init_grid (1);
  T = new scalar[nl];
  //reset({T},0.);
  T_profile = (double *)calloc (nl, sizeof(double));
  for (int l = 0; l < nl; l++)
    T_profile[l] = 0.;
  foreach() {
    foreach_layer() {
      dimensional (T[] == T0);      // here we set T dimension
      T_profile[point.l] += T[];  // here T_profile should inherite T dimensions
    }
  }
  show_dimension (T);
  show_dimension (T_profile[0]);
  show_dimension (T_profile[nl-1]);
  
  free (T_profile);
}
