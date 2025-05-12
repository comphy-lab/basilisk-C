/**
# Test radix sort

~~~gnuplot
set key right outside
set size square
set xlabel 'Index'
set ylabel 'value'
plot 'out' t 'Original', 'log' t 'Sorted'
~~~
 */
#include "radix.h"

#define NUM 5000

int main() {
  int arr[NUM], arrb[NUM], ind[NUM];
  for (int i = 0; i < NUM; i++) {
    arrb[i] = arr[i] = 10000*(noise() + 1);
    printf ("%d %d\n", i, arr[i]);
  }
  radixsort (NUM, arr, ind);
  for (int i = 0; i < NUM; i++) 
    fprintf (stderr, "%d %d\n", i, arrb[ind[i]]);
  
}
