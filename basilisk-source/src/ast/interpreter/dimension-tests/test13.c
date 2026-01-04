/**
Check that 'datasize' and other field properties are defined properly
even when allocating temporary fields within undefined
conditionals. */

int main()
{
  init_grid (1);
  scalar s[];
  int cond;
  if (cond) {
    scalar b[];
  }
  display_value (datasize); // must not be 'unset'
  display_value (_attribute[1].freed); // must not be 'unset'
}
