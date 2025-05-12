/**
# Default arguments do not work for functions returning a pointer

When applying default arguments to a function that returns a pointer, the C
compiler raises an error. This happens because qcc does not substitute the
default argument values during the function call (from the inspection of the
`-source` code).
*/

int aglobal = 1;

int foo1 (int a = 2) {
  printf ("a = %d\n", a); 
  return aglobal;
}

int * foo2 (int a = 2) {
  printf ("a = %d\n", a); 
  return &aglobal;
}

int main() {
  foo1();
  foo2();
}
