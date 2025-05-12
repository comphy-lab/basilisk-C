/**
# Work and torque

We compute work as the scalar product of force and
displacement; $W = \mathbf{F_1} \cdot \mathbf{s}$. Next, we compute the
torque as the vector product of a force and its point of application
vector $\mathbf{\tau} =\mathbf{r} \times \mathbf{F_2}$. */

int main() {
  double F = 1. [1,0]; // Force (not SI, but does not matter)
  double S = 1. [0,1]; // Length
  coord F1 = {F, 0, 0};
  coord s = {S, 0 ,0};
  double W = F1.x*s.x + F1.y*s.y + F1.z*s.z;

  coord F2 = {0, F, 0};
  coord r = {S, 0, 0};
  coord tau = {r.y*F2.z - r.z*F2.y,
	       r.z*F2.x - r.x*F2.z,
	       r.x*F2.y - r.y*F2.x};
  /**
     Finally, we add the z-component of $\mathbf{\tau}$ to the work ($W$)
   */
  double Work_and_torque_component = tau.z + W;
  NOT_UNUSED (Work_and_torque_component);
  if (Work_and_torque_component == 2)
    fprintf (stderr, "Test failed succesfully\n");
  return 1;
  
}