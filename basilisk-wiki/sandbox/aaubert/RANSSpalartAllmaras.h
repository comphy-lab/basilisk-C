/**
   Implementation of the Spalart-Allmaras model */

/**
The equation for the model is 
$$
\nu_t=\hat{\nu}f_{v1}
$$
with $f_{v1}=\frac{\chi^3}{\chi^3+c_{v1}^3}$ where $\chi=\frac{\hat{\nu}}{\nu}$ and $\hat{\nu}$ verify the equation (using Einstein notation)
$$
\frac{\partial\hat{\nu}}{\partial t}+u_j\frac{\partial\hat{\nu}}{\partial x_j}=
c_{b1}\left(1-f_{t2}\right)\hat{S}\hat{\nu}-\left(c_{w1}f_w-\frac{c_{b1}}{\kappa^2}f_{t2}\right) \left(\frac{\hat{\nu}}{d}\right)^2+
\frac{1}{\sigma}\left(\frac{\partial}{\partial x_j}\left(\left(\nu+\hat{\nu}\right)\frac{\partial\hat{\nu}}{\partial x_j}\right)+c_{b2}\frac{\partial\hat{\nu}}{\partial x_i}\frac{\partial\hat{\nu}}{\partial x_i}\right)
$$
where
$$
\hat{S}=\Omega+\frac{\hat{\nu}}{\kappa^2 d^2}f_{v2} ~~~~~   \Omega=\sqrt{2W_{ij}W_{ij}} ~~~~~~  \text{with} ~~~\bold{W}=\frac{1}{2}\left(\nabla \bold{u}-\nabla\bold{u}^t\right) ~~~~ \text{the rotation tensor}
$$
d is the distance to the nearest wall
$$
f_{v2}=1-\frac{\chi}{1+\chi f_{v1}}
$$
$$
f_{w}=g\left(\frac{1+c_{w3}^6}{g^6+c_{w3}^6}\right)^{\frac{1}{6}}  ~~~~ g=r+c_{w2}\left(r^6-r\right) ~~~  \text{with} ~~~~ r=min\left(\frac{\hat{\nu}}{\hat{S}\kappa^2 d^2},10\right)
$$
and
$$
f_{t2}=c_{t3}e^{-c_{t4}\chi^2}
$$
The boundary condition is $\hat{\nu}_{Wall}=0$

The value for the constants is
$$
c_{b1}=0.1355,~ \sigma=\frac{2}{3},~ c_{b2}=0.622,~ \kappa=0.41,~ c_{w2}=0.3,~ c_{w3}=2,~ c_{v1}=7.1, ~c_{t3}=1.2,~ c_{t4}=0.5 ~\text{and} ~c_{w1}=\frac{c_{b1}}{\kappa^2}+\frac{1+c_{b2}}{\sigma}
$$

For stability reason, we use
$$
\hat{S}=\Omega+\bar{S}~~~ \text{when}~~~\bar{S}\geq-c_2\Omega
$$ 
$$
\hat{S}=\Omega+\frac{\Omega\left(c_2^2\Omega+c_3\bar{S}\right)}{\left(c_3-2c_2\right)\Omega-\bar{S}} ~~~\text{when}~~~\bar{S}<-c_2\Omega
$$
with $c_2=0.7$ and $c_3=0.9$
*/

/**
We can use a wall function to prescribe the value of certain quantity where the grid is unresolved, mainly near a wall. The wall function that I use comes from Spalart and Allmaras and is derivated from a solution of the Spalart-Allmaras equation. Introducing the distance from the wall $y$ and the velocity tangential to the wall $u_t$, we have
$$
u_t=u_{\tau}f_{wall}(y^+) ~~\text{with}~~ f_{wall}(y^+)=\bar{B}+c_1log((y^++a_1)^2+b_1^2)-c_2log((y^++a_2)^2+b_2^2)-c_3arctan[y^++a1,b1]-c_4arctan[y^++a_2,b2]
$$
where $y^+=\frac{yu_{\tau}}{\nu}$ is the distance in wall unit and $arctan[x,y]$ is a function name $atan2(y,x)$ in C++ (more information on this function [here](https://en.wikipedia.org/wiki/Atan2)). The dimensionless coefficient have the value
$$
\bar{B}=5.033~~~a_1=8.148~~~a_2=-6.929~~~b_1=7.460~~~b_2=7.468~~~c_1=2.550~~~c_2=1.330~~~c_3=3.599~~~c_4=3.640
$$
Finally, the field transport by the Spalart-Allmaras equation have to be linear close to a wall,
$$
\hat{\nu}=\nu\kappa y^+
$$

The strategy to impose the wall function for embedded boundaries is

* Compute the distance from the wall
* If this distance is large, do nothing. If this distance is too small (typically, $2$ or $3$ times the minimum grid size), find the normal to the solid going through this point
* Find the point on this normal at a distance $d_{IP}$ from the wall. This point is called image point. $d_{IP}$ is chosen such that this point is surrounded by fluid cell (typically $2$ or $3$ times the minimum grid size)
* Interpolate the value of the tangential velocity $u_{t,IP} from the neighbour cell
* Find $u_{\tau}$ such that $u_{t,IP}=u_{\tau}f_{wall}(\frac{d_{IP}u_{tau}}{\nu})$, using Newton method
* Obtain the velocity and the viscosity at the initial point with the $u_{\tau}$ computed

The turbulent viscosity at the wall needs also to be modified. In fact, the wall shear stress will be computed using finite difference $\tau_w=\nu\frac{u_t}{d}$ where $u_t$ is the tangential velocity at the first point of the wall locate at a distance $d$. Due to the unresolved mesh near the wall, this value will be lower than the one excepted for a turbulent boundary layer $\tau_w=u_{\tau}^2$. As the turbulent viscosity is zero at the wall, we can use this value to obtain the correct wall shear stress by solving
$$
(\nu+\nu_t)\frac{u_t}{d}=u_{\tau}^2
$$

We thus need to solve
$$
\hat{\nu}f_{v1}(\frac{\hat{\nu}}{\nu})=\frac{du_{\tau}^2}{u_t}-\nu
$$
We define $C^+=\frac{du_{\tau}^2}{u_t}-\nu$, and we suppose that $C^+$ is positive.
The equation can be written
$$
\frac{\hat{\nu}^4}{\hat{\nu}^3+(c_{v1}\mu)^3}=C^+
$$
and finally, with $\sigma^+=(c_{v1}\mu)^3C^+$, we need to obtain the root of the polynomial
$$
X^4-C^+X^3-\sigma^+
$$
The discriminant of this polynomial is $\Delta=256(-\sigma^+)^3-27(-C^+)^4(-\sigma^+)^2=-256{\sigma^+}^3-27{C^+}^4 {\sigma^+}^2<0$ because $C^+$ and $\sigma^+$ are positive.
So this polynome got two real root and two complex conjugate root. We note $a$ and $b$ the real roots and $\lambda$ and $\bar{\lambda}$ the complex conjugate roots.
Using [Vieta's formulas](https://en.wikipedia.org/wiki/Vieta%27s_formulas), we obtain

$$
\begin{cases}
  |\lambda|^2 ab=-\sigma^+ \\
  2abRe(\lambda)+-|\lambda|^2 (a+b)=0 \\
  ab+|\lambda|^2+2(a+b)Re(\lambda)=0 \\
  a+b+2Re(\lambda)=C^+ 
\end{cases}
$$

With the first equation, we get that $a$ and $b$ do not have the same sign so there is a unique positive solution to the equation for $\hat{mu}$.
Combining the first, the second and the fourth equation, we obtain
$$
-2\sigma^+Re(\lambda)+|\lambda|^4\left(C^+-2Re(\lambda)\right)=0
$$
Thus $C^+-2Re(\lambda)$ and $Re(\lambda)$ need to have the same sign and the only possibility for this is to have $0<Re(\lambda)<\frac{C^+}{2}$.

By combining this four equation, one can obtain $|\lambda|^4$ in term of $Re(\lambda)$
$$
|\lambda|^4=\frac{2\sigma^+Re(\lambda)}{C^+-2Re(\lambda)}
$$
and then the product $ab$
$$
ab=-\frac{2Re(\lambda)\left(C^+-2Re(\lambda)\right)^2}{C^+-4Re(\lambda)}
$$
As $ab<0$, we thus have $Re(\lambda)<\frac{C^+}{4}$.
With the first equation, we obtain an equation for $Re(\lambda)$
$$
8Re(\lambda)^3\left(C^+-2Re(\lambda)\right)^3-\sigma^+\left(C^+-4Re(\lambda)\right)^2=0
$$
with $0<Re(\lambda)<\frac{C^+}{4}$
We define $f(x)=8x^3\left(C^+-2x\right)^3-\sigma^+\left(C^+-4x\right)^2$. We have $f(0)=-\sigma^+C^+<0$, $f(\frac{C^+}{4})>0$ and $f'(x)=24x^2\left(C^+-2x\right)^2\left(C^+-4x\right)+8\sigma^+(C^+-4x)>0$ for $0<x<\frac{C^+}{4}$. Thus the equation have a unique solution between $0$ and $\frac{C^+}{4}$ that we can obtain with a dichotomy.

Finally, we obtain the equation giving the two real root
$$
a^2-a\left(C^+-2Re(\lambda)\right)-\frac{2Re(\lambda)\left(C^+-2Re(\lambda)\right)^2}{C^+-4Re(\lambda)}=0
$$
The discriminant is $\Delta=\left(C^+-2Re(\lambda)\right)^2\left(1+\frac{8Re(\lambda)}{C^+-4Re(\lambda)}\right)>0$.

The final solution is then
$$
\hat{\nu}=\frac{1}{2}\left(C^+-2Re(\lambda)+\sqrt{\Delta}\right)
$$

*/

#include "bcg.h"   //needed for the advection part


#if !defined(WALL_FUNCTION)      //use non linearize version for the wall function if not state otherwise
#define WALL_FUNCTION 0
#endif

#include "functionSA.h"

scalar muhat[];   //field that the Spalart-Allmaras equation propagated
scalar muhatSA[];
scalar muhatSAdiss[];
face vector mut[];    //turbulent viscosity

scalar omega[];

scalar source[];   //source term

scalar molvis[];  //molecular viscosity

double Deltamin, d_IP;   //distance for the image point

double leading; //position of the leading edge
double trailing;  //position of the trailing edge


event defaults(i=0) {
  mu=mut;
  leading=X0;
  trailing=X0+L0;
  
#if TREE
#if EMBED
  foreach_dimension();
  for (scalar s in {muhat,muhatSA,source}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
    s.depends = list_add (s.depends, cs);
  }
#endif // EMBED
#endif // TREE
  
  foreach() {
    dimensional (muhat[] == Delta*Delta/t);
    dimensional (muhatSA[]== Delta*Delta/dt);
    dimensional (muhatSAdiss[]==Delta*Delta/dt);
    source[]=0.;
  }
}

#if EMBED

#undef center_gradient_x
#define center_gradient_x(a) (fs.x[] && fs.x[1] ? (a[1] - a[-1])/(2.*Delta) : \
			     fs.x[1] ? (a[1] - a[])/Delta :		    \
			     fs.x[]  ? (a[] - a[-1])/Delta : 0.)

#undef center_gradient_y
#define center_gradient_y(a) (fs.y[] && fs.y[0,1] ? (a[0,1] - a[0,-1])/(2.*Delta) : \
			      fs.y[0,1] ? (a[0,1] - a[])/Delta :        \
			      fs.y[]  ? (a[] - a[0,-1])/Delta : 0.)

#undef center_gradient_z
#define center_gradient_z(a) (fs.z[] && fs.z[0,0,1] ? (a[0,0,1] - a[0,0,-1])/(2.*Delta) : \
			      fs.z[0,0,1] ? (a[0,0,1] - a[])/Delta :	\
			      fs.z[]  ? (a[] - a[0,0,-1])/Delta : 0.)

#else

#undef center_gradient_x
#define center_gradient_x(a) ((a[1] - a[-1])/(2.*Delta))
			     

#undef center_gradient_y
#define center_gradient_y(a) ((a[0,1] - a[0,-1])/(2.*Delta))

#undef center_gradient_z
#define center_gradient_z(a) ((a[0,0,1] - a[0,0,-1])/(2.*Delta))

#endif


/**
   We first advect the viscosity */

void advection_source() {
  foreach_face() {
    uf.x[]=fm.x[]*face_value(u.x,0);
  }

  advection((scalar *) {muhatSA},uf,dt,(scalar *) {source});

}


/**
   We need to add a correction to the viscosity computed to account for the term in the right hand side */

void correction_nu(double dt) {
  foreach() {
    muhatSA[]+=dt*source[];
  }
}


event reaction_diffusion(i++) {


  advection_source();
  correction_nu(dt);
  
    
  face vector muhat2[];
  face vector amuhat[];
  double chi2;
  foreach_face() {
    chi2=face_value(muhatSAdiss,0)/face_value(molvis,0);//(muhat[]+muhat[-1])/2./molvis;
    muhat2.x[]=fm.x[]*face_value(molvis,0)/sigma*(1.+chi2);
    amuhat.x[]=fm.x[]*face_value(molvis,0)/sigma*(1.+chi2)*face_gradient_x(muhatSA,0);
  }
  
  scalar sourcediff[];
    foreach() {
      if (cm[]>0.) {
	double value=0.;
	double boundary_flux=0.;
#if EMBED
	double coef=embed_flux(point,muhatSA,muhat2,&boundary_flux);
#else
	double coef=0.;
#endif
	foreach_dimension() {
	  value+=(amuhat.x[1]-amuhat.x[])/Delta;
	}
	sourcediff[]=(value-(boundary_flux+coef*muhatSA[]))/(cm[]+SEPS);
      }
      else {
	sourcediff[]=0.;
      }
    }
    foreach() {
      muhatSA[]+=dt*sourcediff[];
    }
    

    correction_nu(-dt);

    foreach() {
      double dstar2=distance_to_wall(x,y);
      if ((dstar2>d_IP)||(x<leading)||(x>trailing)) {
        muhatSA[]=max(muhatSA[],0.);
      }
    }

    trash({source});
    
    double chi;

    double sum;
    double Omega2;                 //norm of the rotation tensor
    double S;
    double d;              //distance to the nearest wall
  
    foreach() {
      sum=0.;
      foreach_dimension() {
        sum+=center_gradient_y(u.x)*(center_gradient_y(u.x)-center_gradient_x(u.y));
#if dimension==3
	sum+=(u.x[0,0,1]-u.x[0,0,-1])/(2.*Delta)*((u.x[0,0,1]-u.x[0,0,-1])/(2.*Delta)-(u.z[1]-u.z[-1])/(2.*Delta));
#endif      
      }
      Omega2=max(sum,0.);
      chi=muhatSA[]/molvis[];
#if WALL
      d=distance_to_wall(x,y);
#else
      d=HUGE;
#endif
      if ((i>100)&&(d<d_IP)) {
	Omega2=sq(omega[]);
      }
      //S=max(sqrt(Omega2)+muhat[]/(pow(kappa*d,2.)+1e-10)*fv2(chi),0.3*sqrt(Omega2));
      S=muhatSA[]/(pow(kappa*d,2.)+1e-10)*fv2(chi);
      if (S>=-c2*sqrt(Omega2)) {
	S+=sqrt(Omega2);
      }
      else {
	S=sqrt(Omega2)+sqrt(Omega2)*(pow(c2,2.)*sqrt(Omega2)+c3*S)/((c3-2*c2)*sqrt(Omega2)-S);
      }
      source[]=cb1*(1-ft2(chi))*S*muhatSA[];

      double r=min(muhatSA[]/(S*pow(kappa*d,2.)+1e-20),10);
      source[]+=-(cw1*fw(r)-cb1/(pow(kappa,2.))*ft2(chi))*pow(muhatSA[]/(d+1e-20),2.);

      foreach_dimension() {
	source[]+=cb2/sigma*sq((muhatSA[1]-muhatSA[-1])/(2.*Delta));//center_gradient(muhat));
      }

      if ((d<d_IP+Deltamin)&&(x>leading)&&(x<trailing)) {
	source[]=0.;
      }

    }

    correction_nu(dt);  
    
    foreach() {
      double dstar3=distance_to_wall(x,y);
      if ((dstar3>d_IP)||(x<leading)||(x>trailing)) {
        muhatSA[]=max(muhatSA[],0.);
      }
    }

    foreach() {
      double dstar=distance_to_wall(x,y);
      if (i>100) {
        if ((dstar>d_IP)||(x<leading)||(x>trailing)) {
          muhat[]=muhatSA[];
          muhatSAdiss[]=muhatSA[];
        }
      }
      else {
        muhat[]=muhatSA[];
        muhatSAdiss[]=muhatSA[];
      }
    }

    foreach_face() {
      chi=face_value(muhat,0)/face_value(molvis,0);
      mut.x[]=fm.x[]*face_value(molvis,0)*(1.+chi*fv1(chi));
    }
}

#if WALL_FUNCTION
#include "wallmodelSA.h"
#endif

    
