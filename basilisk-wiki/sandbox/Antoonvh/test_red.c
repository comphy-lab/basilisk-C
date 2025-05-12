int maxlevel = 8;

int main() {
  // setup grid
  init_grid (N);
  refine ((x + y) > 1.5 && level < maxlevel);
  unrefine (x < 0.5);
  output_cells (stderr);
  // Variables...
  int total = 0, depth = maxlevel + 2;
  int arr[depth], total_arr[1] = {0};
  double topleft = -HUGE, min_arr[5] = {1,1,1,1,1};
  for (int i = 0; i < depth; i++)
    arr[i] = 0;
  // Reductions
  foreach(reduction(+:arr[:depth]) reduction(+:total), // a comma and a comment
	  reduction(max:topleft) reduction(+:total_arr[:1]) 
	  reduction (min:min_arr[:5])) {
    total_arr[0]++;
    total++;
    arr[level]++;
    if (x + y > topleft)
      topleft = x + y;
    for (int i = 0; i < 5; i++) 
      if (i*x < min_arr[i])
	min_arr[i] = i*x;
  }
  // Output
  if (pid() == 0) {
    printf ("# Total: %d and %d\n"
	    "# Top left = %g and min_arr[3] = %g\n",
	    total, total_arr[0], topleft, min_arr[3]);
    for (int i = 0; i < depth; i++)
      printf ("%d %d\n", i, arr[i]); 
  }
}
