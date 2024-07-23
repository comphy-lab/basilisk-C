/**
Checks for undefined pointers.

Should be run with valgrind */

int main()
{
  {
    interpreter_verbosity (2);

    int * a = NULL, b;
    if (b) {
      a = malloc (sizeof (int));
      a[0] = 1; // 'a[0]' is a local value
      free (a);
    }

    /**
    The expression below should in principle give "unallocated array
    access in 'a[0]'", however undefined pointers are not considered,
    so that "undefined" allocation of scalar field works (see
    test16.c). */
    
    a[0];

    coord * c = NULL;    
    if (b) {
      c = malloc (sizeof (coord));
      c->x = 1; // 'c->x' is a local value
      free (c);
    }

    /**
    The expression below should in principle give "undefined structure
    pointer in 'c->x'". Same explanation as above. */
    
    c->x;

    c = malloc (sizeof (coord));
    if (b)
      c->x = 1; // 'c->x' is not a local value
    
    if (c->x); // should give: undefined condition 'c->x'

    free (c);
    
    c = malloc (sizeof (coord));
    if (b) {
      coord * d = c;
      d->x = 1; // 'c->x' is not a local value
    }
    
    if (c->x); // should give: undefined condition 'c->x'

    free (c);

    a = malloc (sizeof (int));
    if (b)
      a[0] = 1; // 'a[0]' is not a local value
    
    if (a[0]); // should give: undefined condition 'a[0]'

    free (a);
    
    a = malloc (sizeof (int));
    if (b) {
      int * d = a;
      d[0] = 1; // 'd[0]' is not a local value
    }
    
    if (a[0]); // should give: undefined condition 'a[0]'

    free (a);
  }
}
