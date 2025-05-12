/**
# Fields are not freed when using 'break'
*/

int main() {
  init_grid (1);
  
  while (true) {
    scalar s[];
    break;
  }

  for (scalar s in all) {
    fprintf (stderr, "%s\n", s.name);
    return 1; // error
  }
}
