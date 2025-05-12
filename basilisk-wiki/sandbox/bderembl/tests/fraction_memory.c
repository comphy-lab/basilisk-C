/**

Illustration of the memory growth

~~~gnuplot mtrace
load "< tail -n2 mtrace"
~~~

*/


#include "fractions.h"
#include "run.h"

scalar topo[];

int main() {
  run();
}


#define MOUNTAIN (exp(-(sq((x)))) - z)

event adjust_topo (i++, i<100) {
  fraction (topo, MOUNTAIN);
}
