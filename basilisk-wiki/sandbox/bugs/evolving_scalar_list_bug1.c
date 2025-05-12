/**
Bug which causes scrambling of scalar list 'evolving'
compile with qcc -g evolving_scalar_list_bug.c -o evolving_scalar_list_bug.exe -lm

Flags -DMTRACE=3 and -catch can be added
*/
#include "two_layer_minimum.h"
#define MAXLEVEL 6

int main()
{
  N = 1 << MAXLEVEL;

  run();
}
/**
As I run through this I output the names of the scalars in the list evolving.  They should be in the order h, hs, u.x, u.y, us.x, us.y
They start in this order but somehow the process of going through a foreach loop scrambles them.  This can be seen by outputing the attribute names after each step.
*/
event init(i=0)
{
  printf("In init: ");
  for (int ii=0; ii<6; ii++) {
    printf (" evolving[ %d ]: %s ", ii, _attribute[evolving[ii].i].name); 
  }
  printf ("\n");
  int counti=0;
  do {
    printf("In do loop between foreach: ");
    for (int ii=0; ii<6; ii++) {
      printf (" evolving[ %d ]: %s ", ii, _attribute[evolving[ii].i].name); 
    }
    printf ("\n");
    foreach() {
      h[]=0;
    }
  } while (counti++<5);
}

