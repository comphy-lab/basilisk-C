/**
# Function pointer and optional arguments do not go together

Optional arguments do not work with function pointers.
*/
// this is fine
int not_optional (double a, int (*fun)(double)) {
  return fun(a);
}
// This is too
int part_optional (double a = 1, double b) {
  return 1;
}
// This isnt...
int optional (double a = 1, int (*fun)(double)) {
  return fun(a);
}
// ...but this is
int optional1 (double a = 1, int (*fun)(double)) {
  return (*fun)(a);
}
int my_fun (double a) {
  return 0;
}

int main() {
  not_optional (2.5, my_fun);
}
