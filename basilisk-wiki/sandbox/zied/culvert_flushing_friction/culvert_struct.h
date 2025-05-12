struct culvert {

/**
# Main variables

The structure “culvert” needs 13 parameters in order to compute the culvert dicharge “q_culvert” at each time-steps: the width, inlet and outlet bed elevations (CBE1,CBE2), the coordinates of the inlet and outlets, the culvert LENGTH, the contraction coefficient μ and inlet, outlet and linear friction coefficients.
*/

  double x_inlet,x_outlet,Dx_inlet,deltac,y_inlet,y_outlet,Dx_outlet;


  double C1,C2,C3,WIDTH,CBE1,CBE2,LENGTH,LHL,FRIC,C56;


  double W,hc,mu,g,q_culvert,h1,h2,hu,z1,z2,S1,S2;

  int FT;
};
/**
#    Culvert modelling
![Fig. 1: source-term coupling approach (green: $p$ inlet cells, purple: $p$ outlet cells, the culvert discharge is equally divided into each cells $\dfrac{Q_{c}}{p }$). ](method1.png)

We first need at each iteration to recover the upstream and downstream free surface water:

$S_{1}=h_{1}+z_{1}$, $S_{2}=h_{2}+z_{2}$, 


*/
double culvert_flux (struct culvert c){

  c.h1 = interpolate(h,c.x_inlet,c.y_inlet);
  c.h2 = interpolate(h,c.x_outlet,c.y_outlet);
  c.z1 = interpolate(zb,c.x_inlet,c.y_inlet)+c.CBE1;
  c.z2 = interpolate(zb,c.x_outlet,c.y_outlet)+c.CBE2;
  c.S1 = c.h1+c.z1;
  c.S2 = c.h2+c.z2;
/**

Using the Bodhaine (Bodhaine, 1968) equations for culvert discharge estimation, we can firstly segregate culvert flow into 6 types (see figure 1) depending on the upstream and downstream water depth :
	
*/

/**

## Flow type 1

-The first flow type is encountered when the critical depth is reached around the culvert entrance, the culvert barrel is not wholly full. It is linked to the slope of the culvert being greater than the critical slope. This case happens rarely in the case of real culverts, it will not be modeled in our simulation.

## Flow type 2

-Parallelly the second flow type is encountered when the critical depth is reached around the culvert exit, the culvert barrel is not wholly full. In this case the slope of the culvert is smaller than the critical slope. The boolean conditions linked to this case are listed below:
	
$$h_{1} < 1.5D\;\; \text{and} \;\; h_{2} \leq h_{c}$$
			
with $D$ the culvert diameter and $hc = \dfrac{2}{3}h_{1}$ the critical water depth in the culvert barrel. Using the simplifications from (carlier, 1972), we can estimate the culvert discharge for this flow type:
	
$$Q_{c} = m h_{c}D\sqrt{2g(S_{1}-z_{2}+hc)},$$
			
with $m$ the contraction coefficient, $S_{1}=h_{1}+z_{1}$ the water free surface at inlet and $z_{2}$ the culvert elevation at the outlet.

## Flow type 3

The third flow type is encountered when culvert barrel is partially full, with the outlet water elevation exceeding critical depth:
	
$$h_{1} < 1.5D\;\; \text{and} \;\; h_{2} > h_{c}$$
			
The culvert discharge for this flow type is then:
	
$$Q_{c} = m (S_{2}-z_{2})D\sqrt{2g(S_{1}-S_{2})},$$
			
with $S_{2}=h_{2}+z_{2}$ the water free surface at the outlet.	

	
*/
  
  
  if (c.S1>c.S2){
    if ((c.S1>c.z1) && (c.S1>c.z2)){
      if ((c.S1-c.z1<1.5*c.WIDTH) && (c.S2<(c.z2+c.WIDTH))){
	if(c.S2>2./3*(c.S1-c.z1)+c.z2){
	  c.FT = 3;}
	else{
	  c.FT=2;}}
      /**
## Flow type 4

Flow type 4 occurs when the culvert outlet is submerged:
	
$$S_{2} > z_{2}+D$$
			
The culvert discharge for this flow type is then:
	
$$Q_{c} = m D^{2}\sqrt{2g(S_{1}-S_{2})},$$	

## Flow type 5

-Flow type 5 occurs when only the inlet is submerged, and the culvert barrel is partially full:
	
$$S_{1} \geq z_{1}+1.5 D \;\; \text{and} \;\; S_{2} \leq z_{2}+D$$
			
The culvert discharge for this flow type is then:
	
$$Q_{c} = m D^{2}\sqrt{2gh_{1}},$$

## Flow type 6
-Flow type 6 is has identical conditions as 5 but describes the flow when the barrel is full, hence we can switch between using type 5 and 6 depending on the culvert geometry:
	
$$S_{1} \geq z_{1}+1.5 D \;\; \text{and} \;\; S_{2} \leq z_{2}+D$$
			
The culvert discharge for this flow type is:
	
$$Q_{c} = m D^{2}\sqrt{2g(S_{1}-(z_{2}+D))},$$

**/

      else if ((c.S1>c.z1+c.WIDTH) && (c.S2>c.z2+c.WIDTH)){
	c.FT = 4;}
      else if ((c.S1-c.z1>=1.5*c.WIDTH) && (c.S2<=c.z2+c.WIDTH)){
	if (c.LENGTH >=c.C56*c.WIDTH){
	  c.FT = 6;}
	else {
	  c.FT = 5;}
      }
      else{
	c.FT=0;}
    }
    else{
      c.FT =0;}
    c.hc = 2./3.*c.h1;
/**
We can now define the culvert discharge from the flow type $FT$:

*/
    if (c.S2>=c.S1){
      c.FT = 0;}
    if (c.FT ==0){
      c.q_culvert =0.;}
    if (c.FT ==2){
      c.q_culvert = c.mu*c.hc*c.W*sqrt(2*G*(c.S1-(c.z2+c.hc)));}
    if (c.FT ==3){
      c.q_culvert = c.mu*(c.S2-c.z2)*c.W*sqrt(2*G*(c.S1-c.S2));}
    if (c.FT ==4){
      c.q_culvert = c.mu*c.W*c.W*sqrt(2*G*(c.S1-c.S2));}
    if (c.FT ==5){
      c.q_culvert = c.mu*c.W*c.W*sqrt(2*G*c.h1);}
    if (c.FT ==6){
      c.q_culvert = c.mu*c.W*c.W*sqrt(2*G*(c.S1-(c.z2+c.W)));}
  }

return c.q_culvert;
}
/**
 
# Links
 
 [Ideal source/well on Basilisk](http://basilisk.fr/sandbox/M1EMN/Exemples/puits.c)
 
 see non viscous dam break with [standard C](http://basilisk.fr/sandbox/M1EMN/Exemples/svdb.c) 
and with [Basilisk](http://basilisk.fr/sandbox/M1EMN/Exemples/damb.c)

# Bibliographie

* [Cui,2019]
" Simulation of Hydraulic Structures in 2D High-Resolution Urban Flood Modeling"

* [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 UPMC


*/