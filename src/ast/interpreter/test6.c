/**
Check named parameters in function calls. */

void bidule (float fc[3])
{
  fc[0];
  fc[1];
  fc[2];
}

int main()
{
  {
    interpreter_verbosity (4);
    bidule (fc = {1,2,3});
  }
}
