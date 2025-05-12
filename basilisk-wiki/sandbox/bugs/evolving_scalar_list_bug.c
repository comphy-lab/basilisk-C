/**
'Bug' which causes scrambling of scalar list 'evolving'
This was a mistake on my behalf not initialising evolving correctly so that it was only a local variable in defaults.  The commented out line 15 was the incorrect initialisation.  When correctly initialised using list_copy (which uses malloc() to globally allocate memory) this works properly.
*/
#define MAXLEVEL 2

scalar zb[], h[], eta[];
vector u[];

#include "predictor-corrector.h"
scalar * evolving = NULL;
int Nlist=3;

event defaults (i = 0){
  //evolving =(scalar *) {h, u};
  evolving = list_copy({h, u});
  foreach()
    for (scalar s in evolving)
      s[]=0.;
  boundary(evolving);
  for (int ii=0; ii<Nlist; ii++) {
    printf (" evolving[ %d ]: %s ", ii, _attribute[evolving[ii].i].name); 
  }
  printf ("\n");
  //  free(evolving);
}

event init (i = 0){
  for (int ii=0; ii<Nlist; ii++) {
    foreach()
        h[]=0.;
    printf (" evolving[ %d ]: %s ", ii, _attribute[evolving[ii].i].name); 
  }
  printf ("\n");
  //  free(evolving);
}

int main()
{
  N = 1 << MAXLEVEL;

  run();
}
