/**
Checks that values which go back to their original value after an
undefined condition are not undefined. */
  
int main()
{
  {
    interpreter_verbosity (2);
    int a = 0, b;
    if (b) {
      ++a;
      --a;
    }
    else
      a = 0;
    if (a); // this is not undefined since a went back to its original value

    a = 0;
    if (b)
      ++a;
    else
      a = 0;
    if (a); // this is undefined since one of the branches modifies a
  }
}
