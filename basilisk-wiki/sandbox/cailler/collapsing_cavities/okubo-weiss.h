/**
Following [Shivamoggi *et al.*, (2022)](#shivamoggi2022), 
we compute the *Okubo--Weiss* criterion in 
cylindrical coordinates $(z, r)$, adapted from eq. (40c) p.11 given for 
3D--axisymmetric flows in spherical polar coordinates $(R,\theta)$:

$$
\lambda^2
=
\left(\dfrac{\partial v_R}{\partial R} \right)^2
+ \cot \theta \dfrac{v_\theta}{R} \dfrac{\partial v_R}{\partial R}
+ \dfrac{1}{R}\dfrac{\partial v_R}{\partial \theta} 
  \dfrac{\partial v_\theta}{\partial R}
- 2 \dfrac{v_\theta}{R} \dfrac{\partial v_\theta}{\partial R}
\equiv 
Q
$$

where $z = R \cos(\theta)$ and $r = R \sin(\theta)$ (that is: 
$R = \sqrt{z^2 + r^2}$ and $\theta = \arctan(r/z)$).
Also, we have that the velocity components in spherical coordinates 
are related to cylindrical ones by:

$$
\left\{\begin{matrix}
v_R & = & v_z \cos \theta + v_r \sin \theta\\ 
v_\theta & = & - v_z \sin \theta + v_r \cos \theta 
\end{matrix}\right.
$$
*/

void compute_lambda2 (scalar lambda2, vector u){
  foreach (){
    double R = sqrt( sq(x) + sq(y) ); 
    double th = atan2(y, x);
    
    double vth = - u.x[]*sin(th) + u.y[]*cos(th);

    // cylindrical derivatives:
    double dvz_dz = ( u.x[1] - u.x[-1] )/( 2.*Delta + SEPS);
    double dvz_dr = ( u.x[0,1] - u.x[0,-1] )/( 2.*Delta + SEPS);
    double dvr_dz = ( u.y[1] - u.y[-1] )/( 2.*Delta + SEPS);
    double dvr_dr = ( u.y[0,1] - u.y[0,-1] )/( 2.*Delta + SEPS);

    // spherical derivatives:
    double dvR_dR = sq(cos(th))*dvz_dz + sq(sin(th))*dvr_dr 
      + cos(th)*sin(th)*( dvr_dz + dvz_dr );

    double dvR_Rdth = - cos(th)*sin(th)*(dvz_dz - dvr_dr)
      - sq(sin(th))*dvr_dz + sq(cos(th))*dvz_dr
      + pow(R, -4.)*(
        - u.x[]*pow(y, 3.) + u.y[]*x*sq(y) - u.x[]*sq(x)*y + u.y[]*pow(y, 3.)
      );

    double dvth_dR = - cos(th)*sin(th)*(dvz_dz - dvr_dr) 
      + sq(cos(th))*dvr_dz - sq(sin(th))*dvz_dr;

    // Okubo-Weiss parameter:
    lambda2[] = ( sq(dvR_dR) + ( vth/( tan(th)*R + SEPS ) )*dvR_dR 
      + dvR_Rdth*dvth_dR - 2.*( vth/(R + SEPS) )*dvth_dR );
  }
}

/** 
We rescale the *Okubo--Weiss* scalar field in the range [-1 ; 1] for better 
visualization. To do so, we use the formula:

$$
x' = 2 \dfrac{x - \min x}{\max x - \min x} - 1
$$
*/

void rescaled_lambda2 (scalar lambda2_rs, scalar lambda2){
  stats s = statsf(lambda2);
  foreach ()
    lambda2_rs[] = 2.*( lambda2[] - s.min )/( s.max - s.min ) - 1.;
}


void compute_lambda2_v2 (scalar lambda2, vector u){
  foreach (){
    if (sqrt(sq(x) + sq(y)) > 3.*Delta){
      double R = sqrt( sq(x) + sq(y) ); 
      double cth = x/( R + SEPS );
      double sth = y/( R + SEPS );
      double tth = y/( x + SEPS );
      
      double vth = - u.x[]*sth + u.y[]*cth;

      // cylindrical derivatives:
      double dvz_dz = ( u.x[1] - u.x[-1] )/( 2.*Delta + SEPS);
      double dvz_dr = ( u.x[0,1] - u.x[0,-1] )/( 2.*Delta + SEPS);
      double dvr_dz = ( u.y[1] - u.y[-1] )/( 2.*Delta + SEPS);
      double dvr_dr = ( u.y[0,1] - u.y[0,-1] )/( 2.*Delta + SEPS);

      // spherical derivatives:
      double dvR_dR = sq(cth)*dvz_dz + sq(sth)*dvr_dr 
        + cth*sth*( dvr_dz + dvz_dr );

      double dvR_Rdth = - cth*sth*(dvz_dz - dvr_dr)
        - sq(sth)*dvr_dz + sq(cth)*dvz_dr
        + pow(R, -4.)*(
          - u.x[]*pow(y, 3.) + u.y[]*x*sq(y) - u.x[]*sq(x)*y + u.y[]*pow(y, 3.)
        );

      double dvth_dR = - cth*sth*(dvz_dz - dvr_dr) 
        + sq(cth)*dvr_dz - sq(sth)*dvz_dr;

      // Okubo-Weiss parameter:
      lambda2[] = ( sq(dvR_dR) + ( vth/( tth*R + SEPS ) )*dvR_dR 
        + dvR_Rdth*dvth_dR - 2.*( vth/(R + SEPS) )*dvth_dR );
    }
  }
}


/** 
# Reference

~~~bib
@article{shivamoggi2022,
  title={The Okubo--Weiss criterion in hydrodynamic flows: geometric aspects and further extension},
  author={Shivamoggi, BK and Van Heijst, GJF and Kamp, LPJ},
  journal={Fluid Dynamics Research},
  volume={54},
  number={1},
  pages={015505},
  year={2022},
  url={https://iopscience.iop.org/article/10.1088/1873-7005/ac495d/pdf},
  doi={10.1088/1873-7005/ac495d},
  publisher={IOP Publishing}
}
~~~
*/