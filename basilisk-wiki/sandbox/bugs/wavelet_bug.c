/**
# `wavelet()` can lockup with MPI
*/
int main() {
  init_grid (N);
  refine (x + y > level/5.); // depth() varies per thread
  scalar s[], w[];
  wavelet (s, w);
}
/**
## A fix

I have `send` this patch

~~~literatec
[Fix a MPI lock bug with `wavelet()`
j.a.v.hooft@gmail.com**20200923124612
 Ignore-this: b1bb5176d1de39ceb4b9238271df49f1
] hunk ./src/grid/multigrid-common.h 89
-  for (int l = depth() - 1; l >= 0; l--) {
+  for (int l = grid->maxdepth - 1; l >= 0; l--) {
~~~
 */
