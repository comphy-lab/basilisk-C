/**
# Global `vertex scalar`-field attributes are not set properly

The following code crashes due to `phi_global`, wheareas `phi_local`
works as intended.
*/
vertex scalar phi_global[];

int Maxlevel = 6;
#define FUNC (sin(x + y) + cos(x*y))

int main() {
  init_grid (1 << Maxlevel);
  vertex scalar phi_local[];
  
  for (scalar s in all)
    if(s.restriction != restriction_vertex)
      fprintf (stderr,
	       "A possible issue may arise with %s as the culprit.\n",
	       s.name);  
  
  foreach_vertex() 
    phi_local[] = phi_global[] = FUNC;
  unrefine (level == Maxlevel - 1);
  foreach_vertex() {
    assert (phi_local[]  == FUNC); // Well done restriction_vertex()!
    assert (phi_global[] == FUNC); // Triggers assertion
  }
}
