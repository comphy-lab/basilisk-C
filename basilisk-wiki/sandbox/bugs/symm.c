/**
# Bug with the tag "symmetric". */

symmetric tensor A[]; 

int main()
{
  init_grid (0);

  foreach() {
    A.x.x[] = 1.;
    A.y.y[] = 5.;
    A.x.y[] = 2.;
  }

  symmetric tensor B[];

  fprintf (stderr, "A: %d %d %d %d\n", A.x.x.i, A.y.y.i, A.x.y.i, A.y.x.i);
  fprintf (stderr, "B: %d %d %d %d\n", B.x.x.i, B.y.y.i, B.x.y.i, B.y.x.i);
  
  foreach() {
    B.x.y[] = 2;
    foreach_dimension() 
      B.x.x[] = A.x.x[];

    fprintf(stderr, " Axx %g Axy %g \n", A.x.x[], A.x.y[]); 
    fprintf(stderr, " Ayx %g Ayy %g \n", A.y.x[], A.y.y[]); 
    fprintf(stderr, " ******************\n"); 
    fprintf(stderr, " Bxx %g Bxy %g \n", B.x.x[], B.x.y[]); 
    fprintf(stderr, " Byx %g Byy %g \n", B.y.x[], B.y.y[]); 
    fprintf(stderr, " ******************\n"); 
  }  
}
