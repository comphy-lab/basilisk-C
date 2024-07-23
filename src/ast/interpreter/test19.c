/** 
Checks that external declarations are properly handled. */

int main()
{
  {
    interpreter_verbosity (5);
    int uf = 3;
    {
      extern int uf;
      fprintf (uf); // must be: fprintf (4)
    }
  }
}

int uf = 4;
