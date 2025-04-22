/**
# Do not report unset values unless they are used in constraints */

int main()
{
  { // interpreter_verbosity (3);
    double s; // unset
    double val = s; // should not report this
    double s1; // unset
    val = s1; // should report this ...
    val == 1 [0]; // ... since it is used here
    init_grid(1);
    scalar s[];
    foreach() {
      if (x > 12)
	s[] = 33.;
      else
	s[] = 10;
      val = s[];
    }
    val == 2;
  }
}
