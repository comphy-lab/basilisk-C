/**
# `_OPENMP` not recognized

When compiling and running the script like so:

`$qcc example.c -lm -fopenmp` --> code should be run like '$qcc -fopenmp example.c -lm'

`$./a.out`

I obtain, 

`i = 0, j = 1`

This is not wat i expected.
*/

int main(){
  int i = 0, j = 0;
#if _OPENMP
  i++;
#endif
@if _OPENMP
  j++;
@endif
  printf("i = %d, j = %d\n", i, j);
}
/**
Noting that `#if _OPENMP` behaves as expected when using gcc `directly'.
*/