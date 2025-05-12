/** 
#The eddy viscosity model by Vreman (2004)
We compute the eddy viscosity from a velocity field according to SGS-model proposed by Vreman in Physics of Fluids (13/10) from 2004 under the title 
*An eddy-viscosity subgrid-scale model for turbulent shear flow: Algebraic theory and applications*. An important aspect of this closure formulation is that the SGS-fluxes vanish in (some) cases of a laminar flow (unlike Standard Smagorinsky). 

Publication external link: [Vreman.pdf](http://www.vremanresearch.nl/Vreman-PF2004-subgridmodel.pdf).

I got my implementation inspiration from [this link @ www.vremanresearch.nl](http://www.vremanresearch.nl/Vreman_Subgridmodel_Fortran.txt)


The function requires a double `Cs` for the (classic) Smagorinsky constant, a (3D) vector field $\mathbf{u}$ and the molecular viscosity as input and it fills the Eddyviscosity field `Evis[]`. 
*/
void eddyviscosity(double Cs, vector u, double molv, scalar Evis){	
  double d1v1, d2v1, d3v1, d1v2, d2v2, d3v2, d1v3, d2v3, d3v3; 
  double b11, b12, b13, b22, b23, b33;
  double abeta, bbeta;	
  foreach(){
    d1v1 = (u.x[1,0,0] - u.x[-1,0,0])/2/Delta;
    d2v1 = (u.x[0,1,0] - u.x[0,-1,0])/2/Delta;
    d3v1 = (u.x[0,0,1] - u.x[0,0,-1])/2/Delta;
    d1v2 = (u.y[1,0,0] - u.y[-1,0,0])/2/Delta;
    d2v2= (u.y[0,1,0] - u.y[0,-1,0])/2/Delta;
    d3v2= (u.y[0,0,1] - u.y[0,0,-1])/2/Delta;
    d1v3= (u.z[1,0,0] - u.z[-1,0,0])/2/Delta;
    d2v3= (u.z[0,1,0] - u.z[0,-1,0])/2/Delta;
    d3v3= (u.z[0,0,1] - u.z[0,0,-1])/2/Delta;
    b11 = sq(Delta)*(d1v1*d1v1 + d2v1*d2v1 + d3v1*d3v1);
    b12 = sq(Delta)*(d1v1*d1v2 + d2v1*d2v2 + d3v1*d3v2);
    b13 = sq(Delta)*(d1v1*d1v3 + d2v1*d2v3 + d3v1*d3v3);
    b22 = sq(Delta)*(d1v2*d1v2 + d2v2*d2v2 + d3v2*d3v2);
    b23 = sq(Delta)*(d1v2*d1v3 + d2v2*d2v3 + d3v2*d3v3);
    b33 = sq(Delta)*(d1v3*d1v3 + d2v3*d2v3 + d3v3*d3v3);
    abeta = sq(d1v1) + sq(d2v1) + sq(d3v1)+
            sq(d1v2) + sq(d2v2) + sq(d3v2)+
            sq(d1v3) + sq(d2v3) + sq(d3v3);
    bbeta = b11*b22 - sq(b12) + b11*b33 - sq(b13) + b22*b33 - sq(b23);
    // Accordig to [1], Evis = 0 for abeta = 0.   
    Evis[] =  (abeta > 10E-5 && bbeta > (abeta/10E6))? 2.5*sq(Cs)*sqrt(bbeta/(abeta)) + molv: molv; 
  }	
}
/** 
## Usage
* [A implementation of this eddy viscosity formulation](SGS.h)

## Test
* [LES of isotropic turbulence](isotropicLES.c)
* There should be more tests

##Reference
[1] Vreman, A. W. "An eddy-viscosity subgrid-scale model for turbulent shear flow: Algebraic theory and applications." Physics of Fluids (1994-present) 16.10 (2004): 3670-3681.
*/ 