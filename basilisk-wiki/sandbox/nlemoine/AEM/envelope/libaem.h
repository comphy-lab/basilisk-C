/**
# Prototypes of FORTRAN functions in the static library *libaem.a*
## FORTRAN subroutine *conformsystem_( )*

The purpose of this function/subroutine is to compute the Jacobian matrices of real potential $\Phi = \mathrm{Re}\!\left\{\Omega\right\}$, stream function $\Psi = \mathrm{Im}\!\left\{\Omega\right\}$, and velocity components $v_x = \mathrm{Re}\!\left\{-\overline{\frac{d\Omega}{dz}}\right\}$ and $v_y = \mathrm{Re}\!\left\{-\overline{\frac{d\Omega}{dz}}\right\}$ at control points on the slit elements, with respect to the degrees of freedom of the elements:
*/
void conformsystem_ ( 	double *, int *,
	       		double *, double *, double *, double *, int *,
			double *, double *, int *, double *, double *, double *, double *,
			int *, double *, double *,
			double *, double *, double *, double *, double *, double *, double *, double *,
			double *, double *, double *, double *,
			double *, int *);
/*
	SUBROUTINE CONFORMSYSTEM ( THETA,ntheta,
     &                             ReZmin,ImZmin,ReZmax,ImZmax,nslit,   
     &                             ReCn,ImCn,ncoef,Q,ReV0,ImV0,Phi0,
     &                             varout,ReZ,ImZ,
     &                             PHI,PSI,Vx,Vy,J_PHI,J_PSI,J_Vx,J_Vy,
! work arrays:
     &                             J_PHI1,J_PSI1,J_Vx1,J_Vy1,
     &                             Xdof,ndof)
*/

/**
Following the approach of [Steward (2015)](https://doi.org/10.1002/2015WR017526), we consider a collection of slit elements, each element $j$ creating a complex potential $\Omega_j(z)$ at location $z=x+iy$ in the $\mathbb{C}$-plane, with the form:
$$\Omega_j(z)= Q_j\times\left(\frac{\ln Z_j(z)}{2\pi}-1\right) + \sum_{k=1}^{n_\mathrm{coef}} c_{j,k}Z_j^{-k}(z)$$
$Q_j$ is a sink term, and the $c_{j,k}$ are the complex coefficients of a Laurent series that can be adjusted to match various boundary conditions. 
The mapping function $Z_j(z)$ for each slit element is given by:

$$\begin{cases}
\zeta_j(z;a_j,b_j)= \displaystyle\frac{z-\frac{1}{2}(a_j+b_j)}{b_j-a_j} & \\
& \\
Z_j(z) = \zeta_j(z)+\sqrt{\zeta_j(z)+1}\sqrt{\zeta_j(z)-1} & \\
\end{cases}$$
 
where $a_j$ is the complex affix of the start point of slit element $j$, and $b_j$ the affix of its end point. The function $\zeta\mapsto \zeta + \sqrt{\zeta+1}\sqrt{\zeta-1}$ is the Joukowsky transformation (Joukowsky, 1910).

By superposition, the total (complex) potential at location $z$ in the complex plane is simply the sum of the contributions of all slit elements in the domain:
$$
\Omega_\mathrm{tot}(z) = \Phi_\mathrm{tot}(z) + i\,\Psi_\mathrm{tot}(z) = \Phi_0 - \overline{v_0}\,z + \sum_{j=1}^{n_\textrm{slit}} \Omega_j(z)
$$
where $v_0 = v_{0,x}+i\,v_{0,y}$ is a complex number giving the 'background' uniform flow, $\overline{v_0}$ its complex conjugate, and $\Phi_0$ an offset value for the real potential.

Denoting $n_\textrm{slit}$ the number of slit elements in the system and $n_\textrm{coef}$ the number of complex-valued coefficients in the Laurent series, there are $(2 n_\textrm{coef} +1)n_\textrm{slit}$ real degrees of freedom. We define $n_\theta$ control points along each slit element, such that:

$$\begin{cases}
\theta_m = \pi\frac{m-\frac{1}{2}}{n_\theta}\qquad 1\leq m\leq n_\theta & \\
& \\
z(\theta_m) = \frac{1}{2}(b_j-a_j)\cos(\theta_m) + \frac{1}{2}(a_j+b_j)& \\
\end{cases}$$

Then the matrices $J_\Phi$, $J_\Psi$, $J_{v_x}$, and $J_{v_y}$ are of the form:


$$\small\begin{array}{r|c|}
& \overbrace{
\mathrm{Re}\!\left\{c_1\right\}\ 
\mathrm{Im}\!\left\{c_1\right\}
\ \cdots\ 
\mathrm{Re}\!\left\{c_{n_\textrm{coef}}\right\}\ 
\mathrm{Im}\!\left\{c_{n_\textrm{coef}}\right\}\ 
Q
}^{\small\textrm{slit } 1}
\quad\overbrace{
\mathrm{Re}\!\left\{c_1\right\}\ 
\mathrm{Im}\!\left\{c_1\right\}
\ \cdots\ 
\mathrm{Re}\!\left\{c_{n_\textrm{coef}}\right\}\ 
\mathrm{Im}\!\left\{c_{n_\textrm{coef}}\right\}\ 
Q
}^{\small\textrm{slit } 2}
\qquad\cdots\qquad
\overbrace{
\mathrm{Re}\!\left\{c_1\right\}\ 
\mathrm{Im}\!\left\{c_1\right\}
\ \cdots\ 
\mathrm{Re}\!\left\{c_{n_\textrm{coef}}\right\}\ 
\mathrm{Im}\!\left\{c_{n_\textrm{coef}}\right\}\ 
Q
}^{\small\textrm{slit } n_\textrm{slit}}
\\
& \\
\hline
\textrm{slit }1 \left\{\begin{array}{c}
\theta_1 \\
\theta_2 \\
\vdots \\
\theta_{n_\theta}
\end{array}\right.
& \\
\textrm{slit }2 \left\{\begin{array}{c}
\theta_1 \\
\theta_2 \\
\vdots \\
\theta_{n_\theta}
\end{array}\right. & \Large\texttt{J(nslit\,*\,ntheta,(2\,*\,ncoef+1)*\,nslit)}\\
\vdots\quad\left\{
\begin{array}{c}
~\\
~\\
~\\
~\\
\end{array}
\right.\quad & \\
\textrm{slit }n_\textrm{slit} \left\{\begin{array}{c}
\theta_1 \\
\theta_2 \\
\vdots \\
\theta_{n_\theta}
\end{array}\right. & \\
\hline
\end{array}$$
*/

/**
![](https://comptes-rendus.academie-sciences.fr/geoscience/article/CRGEOS_2023__355_S1_79_0/jpg/src/tex/figures/fig05.jpg){width="1200px"}
*/

/**
General form of the potential created by a slit element. In each figure, the gray levels show the value of $\Phi = \mathrm{Re}\left\{\Omega(z)\right\}$ while the dashed white lines are contours of the stream function
$\Psi = \mathrm{Im}\left\{\Omega(z)\right\}$. Black arrows show the flux vector $q=-\overline{\frac{d\Omega}{dz}}$. (Left) Potential created by a pure sink. (Center) Potential created by a single, real-valued coefficient 
in the Laurent series; note that the real potential $\Phi$ is continuous between the 'top'  and 'bottom' of the slit, while the stream function $\Psi$ is discontinuous. (Right) Potential created by a single, pure imaginary coefficient 
in the Laurent series; there is now a jump in the real potential $\Phi$ between the two sides of the slit while the stream function, in turn, is continuous.

## FORTRAN subroutine *slit_( )*
*/

void slit_ ( double *, double *, int *,
	     double *, double *, double *, double *, int *,
             double *, double *, int *,double *,double *,double *,double *,
	     int *,
             double *, double *, double *, double *,double *,double *,double *,double *);
/*
	SUBROUTINE SLIT (ReZ,ImZ,nz,
     &                   ReZmin,ImZmin,ReZmax,ImZmax,nslit,
     &                   ReCn,ImCn,ncoef,Q,ReV0,ImV0,Phi0,
     &                   varout,
     &                   PHI,PSI,Vx,Vy,J_PHI,J_PSI,J_Vx,J_Vy)
*/

/**
## FORTRAN subroutine *arealsink_( )* 
*/

void arealsink_ ( double *, double *, int *,   
                  double *, double *, int *, double *,
                  double *, double *, double *, double *);
/*
	SUBROUTINE AREALSINK ( Xvert,Yvert,nvert,   
     &                         Xquery,Yquery,nquery,gamma0,
! output arrays:
     &                         PHI,Vx,Vy,A)
*/
/**
## FORTRAN subroutine *inpolygon_( )* 
*/

void inpolygon_ ( double *, double *, int *,
                  double *, double *, int *,
                  double *);
/*
        SUBROUTINE INPOLYGON( Xvert,Yvert,nvert,   
     &                        Xquery,Yquery,nquery,
! Output array
     &                        RES)
*/



/**
## LAPACK FORTRAN subroutine *dgels_( )*
*/

void dgels_ ( char *, int *, int *, int *, double *, int *, double *, int *, double *, int *, int *);

/*
      SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
     $                  INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*/