/**
# Check that imbricated returns work properly */

int func()
{
  for (int i = 0; i < 3; i++)
    return i;
  return -1;
}

int main()
{
  display_value (func());
}
