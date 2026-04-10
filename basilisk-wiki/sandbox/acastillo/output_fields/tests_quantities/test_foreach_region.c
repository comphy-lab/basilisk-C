/** 
This is just to verify that foreach_region uses the nearest neighbors.
*/

#define MAXLEVEL 4

int main()
{
  L0 = 1.0;
  X0 = Y0 = Z0 = -L0 / 2;
  N = 1 << MAXLEVEL;
  init_grid(N);

  /** We load the fields from */
  scalar f[], p[];
  vector u[];
  foreach(){
    f[] = x;
    p[] = y;
    u.x[] = y;
    u.y[] = x;
  }
  
  coord pp;
  coord box[2] = {{-0.5, -0.5}, {0.5, 0.5}};
  coord n = {16, 16};
  foreach_region (pp, box, n) {
    printf("%f %f %f %f %f %f\n", x, y, pp.x, pp.y, f[], p[]);
  }
}