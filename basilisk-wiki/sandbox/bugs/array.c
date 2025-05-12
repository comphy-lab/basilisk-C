/**
# Arrays and scalars cannot be mixed 

This is because the preprocessor does not realise that the local array takes precedence over the scalar. */

int main() {
  init_grid (64);
  scalar s[], p[];
  foreach()
    s[] = p[] = 1.;
  foreach() {
    double p[3]; // renaming p to something else fixes the problem
    for (int i = 0; i < 2; i++)
      p[i] = s[];
  }
}