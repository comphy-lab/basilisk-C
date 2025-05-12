/**
# Many `bid`s

`n` to be exact 

![`n` = 4 (randomized)](bids/s.png) 
*/
#include "poisson.h"
#include "utils.h"

int n = 4;

scalar s[], g[];

s[left] = dirichlet (0.);
s[right] = dirichlet (0);

bid * id;
double bid_value = 1;
  
int main() {
  id = (bid*) malloc (n*sizeof(bid));
  for (int i = 0; i < n; i++)
    id[i] = new_bid();
  
  init_grid (N);
  
  srand(time(NULL));
  for (int i = 0; i < n; i++) {
    double xp = (noise() + 1.)/4. + 0.25 , yp = (noise() + 1.)/4 + 0.25;
    double R = (noise() + 1)/20. + 0.05;
    mask (sq(x - xp) + sq(y - yp) < sq(R) ? id[i] : none);
    bid temp = id[i]; //A place holder to set the identical conditions
    s[temp] = dirichlet (bid_value);
  }
  
  refine (level < 7);
  poisson (s, g);
  output_ppm (s, file = "s.png", n = 256, min = 0, max = 1);
  free (id);
}
  
