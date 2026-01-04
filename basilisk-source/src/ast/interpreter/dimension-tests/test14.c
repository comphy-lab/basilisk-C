/**
Checks that point variables are correctly defined. */

int main()
{
  init_grid (1);
  foreach_face()
    display_value (point);
}
