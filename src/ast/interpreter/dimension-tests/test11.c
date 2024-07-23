/**
Checks that the dimension of "out of bounds" array elements is the
same as that of initial elements. */

int main()
{

  /**
  The default maximum number of iterations of the interpreter is 32,
  so the last elements of the array will not be initialised. */
  
  double a[100];
  for (int i = 0; i < 100; i++)
    a[i] = 2. [1];
  
  /**
  The undefined condition below will return an error if a[0] and a[99]
  don't have the same dimensions. */
  
  double cond, val;
  if (cond)
    val = a[0];
  else
    val = a[99];
}
