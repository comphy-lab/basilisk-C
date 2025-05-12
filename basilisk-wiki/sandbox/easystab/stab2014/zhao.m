 %{
# Contributions of ZHAO Dan
In this contribution page, I have worked on the equations of advection in 1D and 2D.  
There are two codes that I have mainly 
refered to : [vibrating_string.m](/sandbox/easystab/vibrating_string.m) and [poisson2D.m](/sandbox/easystab/poisson2D.m).
 %}
 
  %{
# Equation advection 1D
  My code : [advection_1D.m](/sandbox/easystab/advection_1D.m)   
  I do the marching in time of the advection equation $f_t+Uf_x=0$ with just one boundary condition at the left end (with U possitive)  
  This is an equation with one time derivative, to solve we need just one boundary condition. I have valide my curve with a theoretical solution。
  
  <center>
  ![visualisation of the advection_1D solution ](/advection_1d.gif)
  ![validation](/vibrating2.png)
 </center>
 %}
 
 %{
# Equation advection 2D
  My code: [advection_2D.m](/sandbox/easystab/advection_2D.m)   
  I do the marching in time of the advection equation $f_t+Uf_x+Vf_y=0$ in 2D (with U possitive and V=0)  
  In advection 2D, I have refered to the code [poisson2D.m](/sandbox/easystab/poisson2D.m) and use its methode to extract the boundary cell locations in order to set the boundary conditions more efficiently.
  <center>
 ![visualisation of the advection_2D solution ](/advection_2d.gif)
 ![validation](/validation.png)
 %}