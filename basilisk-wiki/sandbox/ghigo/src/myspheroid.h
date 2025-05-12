/**
# Utility functions for spheroids

## Geometry of a spheroid

The equation for a spheroid with an $x$-symmetry is given by:
$$
\frac{x^2}{a^2} + \frac{y^2 + z^2}{r^2} - 1 = 0,
$$
where $a$ is the polar radius along the symmetry axis ($x$-direction
here) and $r$ is the equatorial radius. We also denote $d=2r$ the
equatorial diameter of the spheroid.

The aspect ratio:
$$
\gamma = \frac{a}{r},
$$
is the ratio of the polar to equatorial radii. It is a usefull
quantity to characterize the shape of a spheroid. We distinguish three
cases:

* $\gamma > 1$: prolate spheroid;
* $\gamma < 1$: oblate spheroid;
* $\gamma = 1$: sphere.

Other usefull quantities are the flattening (or oblateness):
$$
f = \frac{r - a}{r} = 1 - \gamma,
$$
and the (first) eccentricity:
$$
e = \left\{
\begin{aligned}
&\mathrm{if}\: \gamma < 1: \quad \sqrt{1 - \gamma^2} \\
&\mathrm{if}\: \gamma > 1: \quad \sqrt{1 - \frac{1}{\gamma^2}},
\end{aligned}
\right.
$$
of the spheroid.

In the following, we choose to describe the spheroid's geometry in the
$\left[r, \gamma\right]$ phase space.
*/

double spheroid_flattening (double g)
{
  return 1. - g;
}

double spheroid_eccentricity (double g)
{
  double e = 0.;
  if (g <= 1.)
    e = sqrt (1. - sq (g));
  else
    e = sqrt (1. - 1./sq (g));
  return e;
}

/**
## Volume of a spheroid

The volume of a spheroid is given by:
$$
V = \frac{4}{3}\pi r^3 \gamma.
$$
*/

double spheroid_volume (double r, double g)
{
  return 4./3.*(M_PI)*cube (r)*g;
}

/**
Another usefull quantity linked to the volume of a spheroid is the
radisu $r_s$, or diameter $d_s = 2r_s$, of the sphere of equivalent
volume:
$$
\begin{aligned}
& V_s = V \\
\Leftrightarrow & \frac{4}{3}\pi r_s^3 = \frac{4}{3}\pi r^3 \gamma \\
\Leftrightarrow & r_s = r \gamma^{\frac{1}{3}}.
\end{aligned}
$$
*/

double spheroid_diameter_sphere (double r, double g)
{
  return 2.*r*pow (g, 1./3.);
}

/**
## Surface area of a spheroid

The formula for the surface area of a spheroid varies depending on
whether the spheroid is oblate ($\gamma < 1$) or prolate ($\gamma >
1$):
$$
S = 4\pi r^2 S_{\mathrm{corr}},
$$
with
$$
S_{\mathrm{corr}} =
\left\{
\begin{aligned}
&\mathrm{if}\: \gamma < 1: \quad \frac{1 + \frac{1 - e^2}{e}\mathrm{arctanh}\, e}{2} \\
&\mathrm{if}\: \gamma > 1: \quad \frac{1 + \frac{\gamma}{e}\arcsin e}{2}.
\end{aligned}
\right.
$$
Note here that the expression for the eccentricity $e$ changes
depending on whether the spheroid is oblate and prolate. */

double spheroid_surface_correction (double g)
{
  double e = spheroid_eccentricity (g);
  double s = 1.;
  if (g < 1.)
    s *= 0.5*(1. + (1. - sq (e))/e*atanh (e));
  else if (g > 1.)
    s *= 0.5*(1. + g/e*asin (e));
  return s;
}

double spheroid_surface (double r, double g)
{
  return 4.*(M_PI)*sq (r)*spheroid_surface_correction (g);
}

/**
Other usefull quantities linked to the surface area of a spheroid are
its sphericity $\varphi$ and crosswise sphericity $\varphi_{\perp}$:

* the sphericity $\varphi$ is the ratio of the surface area of the
  sphere of equivalent volume over the surface area of the spheroid:
  $$
  \varphi = \frac{4 \pi r_s^2}{S} = \frac{4 \pi r^2 \gamma^{\frac{2}{3}}}{4 \pi r^2 S_{\mathrm{corr}}} = \frac{\gamma^{\frac{2}{3}}}{S_{\mathrm{corr}}};
  $$
* the crosswise sphericity $\varphi_{\perp}$ is the ratio of the
  cross-sectional area of the sphere of equivalent volume over the
  projected cross-sectional area of the spheroid perpendicular to the
  flow direction:
  $$
  \varphi_{\perp} = \frac{\pi r_s^{2}}{\pi r^2 \cos^2\phi} = \frac{\gamma^{\frac{2}{3}}}{\cos^2\phi},
  $$
  where $\phi$ is the incidence angle characterizing the rotation
  along the $z$-axis. */

double spheroid_sphericity (double g)
{
  return pow (g, 2./3.)/spheroid_surface_correction (g);
}

double spheroid_crosswise_sphericity (double g, double phi)
{
  return pow (g, 2./3.)/sq (cos (phi));
}

/**
## Correlations for the drag, lift and pitching torque coefficients of a spheroid

We define the drag coefficient of a spheroid in a uniform flow in the
$x$-direction as:
$$
C_D = \frac{F_x}{\frac{1}{2}\rho u^2 S_{s,\perp}},
$$
where $S_{s,\perp} = \frac{\pi}{4}d_s^2$ is the cross-sectional area
perpendicular to the flow direction of the sphere of equivalent
volume. Equivalently, we also define the drag correction function:
$$
f_D = \frac{C_D}{C_{D,\mathrm{Stokes}}},
$$
where $C_{D,\mathrm{Stokes}}$ is the analytical drag for creeping flow
conditions ($Re=0$).

Similarly, we define the lift coefficient as:
$$
C_L = \frac{F_{y\text{ or }z}}{\frac{1}{2}\rho u^2 S_{s,\perp}}.
$$
Finally, we define the pitching torque coefficient (along the
$z$-axis) as:
$$
C_T = \frac{T_z}{\frac{1}{2}\rho u^2 S_{s,\perp} \frac{d_s}{2}}.
$$

Many correlations exist in the litterature to predict the drag, lift
and pitching torque coefficients of a spheroid in a uniform flow as a
function typically of the Reynolds number $Re$, the aspect ration
$\gamma$ and the incidence angle $\phi$. We provide expressions below
for a few of them, and focus in particular on those that provide a
dependance with the incidence angle $\phi$. */

/**
#### Correlations for a sphere

A sphere is a spheroid with an aspect ratio $\gamma = 1$. As a
reference, we therefore provide in the following well-know correlation
for the drag coefficient of a sphere as a function of the Reynolds
number $Re$.

In creeping flow conditions, the drag coefficient is given by the
Stokes formula:
$$
C_{D,\mathrm{Stokes}} = \frac{24}{Re}.
$$
*/

double drag_sphere_Stokes (double Re)
{
  return 24./Re;
}

/**
[Schiller and Naumann, 1933](#Schiller1933) then proposed the
following correlation for the large range of Reynolds number $0 \leq
Re \leq 1000$:
$$
f_{D,SN} =
\left\{
\begin{aligned}
&\mathrm{if}\:  Re \leq 1000:& 1 + 0.15 Re^{0.687}\\
&\mathrm{else}:& 0.44.
\end{aligned}
\right.
$$
*/

double drag_sphere_Schiller_Naumann (double Re)
{
  double f = 0.44;
  if (Re <= 1000)
    f = 1. + 0.15*pow (Re, 0.687);
  return drag_sphere_Stokes (Re)*f;
}

/**
[Clift et al., 2005](#Clift2005) finally proposed the following piecewise defined
standard drag curve after carefully reviewing experimental and
numerical results for the range of Reynolds number $0 \leq Re \leq
260$:
$$
f_{D,\mathrm{Clift}} = 
\left\{
\begin{aligned}
&\mathrm{if}\: Re < 0.01:& 1 + \frac{Re}{128} \\
&\mathrm{if}\: 0.01 \leq Re \leq 20:& 1 + 0.1315 Re^{0.82 - 0.05\log_{10} Re} \\
&\mathrm{if}\: 20 < Re \leq 260:& 1 + 0.935 Re^{0.6305}.
\end{aligned}
\right.
$$
We do not show here the additional correlations proposed by [Clift et
al., 2005](#Clift2005) for larger Reynolds numbers.
*/

double drag_sphere_Clift (double Re)
{
  double f = 1.;
  if (Re < 0.01)
    f = 1. + Re/128.;
  else if (Re >= 0.01 && Re <= 20.) {
    double w = log10 (Re);
    f = 1. + 0.1315*pow (Re, 0.82 - 0.05*w);
  }
  else if (Re > 20. && Re <= 260)
    f = 1. + 0.935*pow (Re, 0.6305);
  return drag_sphere_Stokes (Re)*f;
}

/**
#### Low Reynolds numbers correlations for a spheroid

Several correlations are summarized in [Masliyah and Epstein,
1970](#Masliyah1970) for the drag coefficient of a spheroid at low
Reynolds numbers:

* Stokes' solution:
  $$
  C_{D,\mathrm{Stokes}} = \frac{24}{Re} K,
  $$
* Oseen's solution:
  $$
  C_D = \frac{24}{Re}K\left(1 + \frac{3}{16} K Re \right),
  $$
* Breach's solution:
  $$
  C_D = \frac{24}{Re}K\left(1 + \frac{3}{16} K Re + \frac{9}{160} \left(K Re\right)^2 \log\left(\frac{Re}{2}\right)\right),
  $$

where $K$ is a correction factor derived in [Happel and Brenner,
1967?](#Happel2012) for both a prolate and an oblate spheroid aligned
with the flow direction:
$$
K = \left\{
\begin{aligned}
&\mathrm{if}\: \gamma < 1: \frac{1}{\frac{3}{4}\sqrt{\lambda_0^2 + 1}\left(\lambda_0 - \left(\lambda_0^2 - 1\right) \cot^{-1} \left(\lambda_0\right)\right)} \quad \mathrm{with} \quad \lambda_0 = \frac{\gamma}{e} \\
&\mathrm{if}\: \gamma > 1: -\frac{1}{\frac{3}{4}\sqrt{\lambda_0^2 - 1}\left(\lambda_0 - \left(\lambda_0^2 + 1\right) \coth^{-1} \left(\lambda_0\right)\right)} \quad \mathrm{with} \quad \lambda_0 = \frac{1}{e}.
\end{aligned}
\right.
$$
Note here that the expression of the eccentricity $e$ changes between
oblate and prolate spheroids.
*/

double drag_correction_Happel (double g)
{
  double e = spheroid_eccentricity (g);
  double K = 1.;
  if (g < 1) {
    double l0 = g/e;
    K = 4./3./sqrt (sq (l0) + 1.)/
      (l0 - (sq (l0) - 1.)*((M_PI)/2. - atan (l0)));
  }
  else if (g > 1) {
    double l0 = 1./e;
    K = -4./3./sqrt (sq (l0) - 1.)/
      (l0 - (sq (l0) + 1.)*(0.5*(log (1. + 1./l0) - log (1. - 1/l0))));
  }
  return K;
}

double drag_Stokes (double g, double Re)
{
  double K = drag_correction_Happel (g);
  return 24./Re*K;
}

double drag_Oseen (double g, double Re)
{
  double K = drag_correction_Happel (g);
  return 24./Re*K*(1. + 3./16.*K*Re);
}

double drag_Breach (double g, double Re)
{
  double K = drag_correction_Happel (g);
  return 24./Re*K*(1. + 3./16.*K*Re + 9./160.*sq (K*Re)*log (Re/2.));
}

/**
#### Rosendahl, 2000

The correlation provided in [Rosendahl, 2000](#Rosendahl2000) uses a
blending between two drag coefficients, obtained when the spheroid is
aligned with the flow direction ($\phi = 0$) and when the spheroid is
perpendicular to the flow direction ($\phi = \frac{\pi}{2}$):
$$
C_D = C_{D,\phi = 0} + \left(C_{D,\phi=90} - C_{D,\phi=0}\right)\sin^{3}\phi.
$$
where the values of $C_{D,\phi = 0}$ and $C_{D,\phi = 90}$ are taken
from previous experimental or numerical results. Typically, we use
here the correlation from [Holzer and Sommerfeld, 2008](#Holzer2008).
*/

/**
#### Holzer and Sommerfeld, 2008

The correlation provided in [Holzer and Sommerfeld, 2008](#Holzer2008)
expresses the drag as a function of the Reynolds number $Re$ and the
sphericity $\varphi$. It also accounts for the orientation of the
spheroid (i.e. incidence angle $\phi$) using the crosswise sphericity
$\varphi_{\perp}$. It is the sum of a correlation representing the
drag in the Stokes regime and a correlation representing the drag in
the Newton regime:
$$
C_D = \frac{1}{Re_{p}} \left(\frac{8}{\sqrt{\varphi_{\perp}}} +
\frac{16}{\sqrt{\varphi}} +
\frac{3\sqrt{Re_{p}}}{{\varphi}^{\frac{3}{4}}} \right) +
\frac{0.42\times 10^{0.4\left(-\log
\varphi\right)^{0.2}}}{\varphi_{\perp}}.
$$
*/

double drag_Holzer2008 (double g, double angle, double Re)
{
  double phi = spheroid_sphericity (g);
  double phi_perp = spheroid_crosswise_sphericity (g, angle);
  return 1./Re*(
		8./sqrt (phi_perp) +
		16./sqrt (phi) +
		3.*sqrt (Re)/pow (phi, 3./4.)
		) +
    0.42*pow (10., 0.4*pow (-log (phi), 0.2))/phi_perp;
}

/**
#### Zastawny et al., 2012 (immersed boundary method)

The correlation for the drag coefficient provided in [Zastawny et al.,
2012](#Zastawny2012) uses a blending between two drag coefficients,
obtained when the spheroid is aligned with the flow direction ($\phi =
0$) and when the spheroid is perpendicular to the flow direction
($\phi = \frac{\pi}{2}$):
$$
C_D = C_{D,\phi = 0} + \left(C_{D,\phi=90} - C_{D,\phi=0}\right)\sin^{a_0}\phi,
$$
with:
$$
C_{D,\phi=0} = \frac{a_1}{Re^{a_2}} + \frac{a_3}{Re^{a_4}} \quad
\mathrm{and} \quad C_{D,\phi=90} = \frac{a_5}{Re^{a_6}} +
\frac{a_7}{Re^{a_8}}.
$$

The correlation for the lift coefficient provided in [Zastawny et al.,
2012](#Zastawny2012) writes:
$$
C_L = \left(\frac{b_1}{Re^{b_2}} +
\frac{b_3}{Re^{b_4}}\right)\left(\sin \phi \right)^{b_5 + b_6 Re^{b_7}}
\left(\cos \phi\right)^{b_8 + b_9 Re^{b_{10}}}.
$$

Finally, the correlation for the pitching torque coefficient provided
in [Zastawny et al., 2012](#Zastawny2012) writes as the one for the
lift coefficient:
$$
C_T = \left(\frac{c_1}{Re^{c_2}} +
\frac{c_3}{Re^{c_4}}\right)\left(\sin \phi \right)^{c_5 + c_6 Re^{c_7}}
\left(\cos \phi\right)^{c_8 + c_9 Re^{c_{10}}}.
$$

Note that the coefficients $a_i$, $b_i$ and $c_i$ are provided in
[Zastawny et al., 2012](#Zastawny2012) only for two spheroids,
respectively defined by $\gamma = 5/2$ and $\gamma = 5/4$, and are
valid for $0.1 \leq Re \leq 300$ and $0 \leq \phi \leq 90$. */

double drag_Zastawny2012 (double g, double angle, double Re)
{
  double a0 = 0., a1 = 0., a2 = 0., a3 = 0., a4 = 0.;
  double a5 = 0., a6 = 0., a7 = 0., a8 = 0.;
  if (g == 5./2.) {
    a0 = 2., a1 = 5.1, a2 = 0.48, a3 = 15.52, a4 = 1.05;
    a5 = 24.68, a6 = 0.98, a7 = 3.19, a8 = 0.21;
  }
  else if (g == 5./4.) {
    a0 = 1.95, a1 = 18.12, a2 = 1.023, a3 = 4.26, a4 = 0.384;
    a5 = 21.52, a6 = 0.99, a7 = 2.86, a8 = 0.26;
  }
  double CD0  = a1/pow (Re, a2) + a3/pow (Re, a4);
  double CD90 = a5/pow (Re, a6) + a7/pow (Re, a8);
  return CD0 + (CD90 - CD0)*pow (sin (angle), a0);
}

double lift_Zastawny2012 (double g, double angle, double Re)
{
  double b1 = 0., b2 = 0., b3 = 0., b4 = 0., b5 = 0.;
  double b6 = 0., b7 = 0., b8 = 0., b9 = 0., b10 = 0.;
  if (g == 5./2.) {
    b1 = 6.079, b2 = 0.898, b3 = 0.704, b4 = -0.028, b5 = 1.067;
    b6 = 0.0025, b7 = 0.818, b8 = 1.049, b9 = 0., b10 = 0.;
  }
  else if (g == 5./4.) {
    b1 = 0.083, b2 = -0.21, b3 = 1.582, b4 = 0.851, b5 = 1.842;
    b6 = -0.802, b7 = -0.006, b8 = 0.874, b9 = 0.009, b10 = 0.57;
  }
  return (b1/pow (Re, b2) + b3/pow (Re, b4))*
    pow (sin (angle), b5 + b6*pow(Re, b7))*
    pow (cos (angle), b8 + b9*pow(Re, b10));
}

double torque_Zastawny2012 (double g, double angle, double Re)
{
  double c1 = 0., c2 = 0., c3 = 0., c4 = 0., c5 = 0.;
  double c6 = 0., c7 = 0., c8 = 0., c9 = 0., c10 = 0.;
  if (g == 5./2.) {
    c1 = 2.078, c2 = 0.279, c3 = 0.372, c4 = 0.018, c5 = 0.98;
    c6 = 0., c7 = 0., c8 = 1., c9 = 0., c10 = 0.;
  }
  else if (g == 5./4.) {
    c1 = 0.935, c2 = 0.146, c3 = -0.469, c4 = 0.145, c5 = 0.116;
    c6 = 0.748, c7 = 0.041, c8 = 0.221, c9 = 0.657, c10 = 0.044;
  }
  return (c1/pow (Re, c2) + c3/pow (Re, c4))*
    pow (sin (angle), c5 + c6*pow(Re, c7))*
    pow (cos (angle), c8 + c9*pow(Re, c10));
}

/**
#### Ouchenne et al., 2016

The correlation provided in [Ouchenne et al., 2016](#Ouchenne2016)
follows the same approach as in [Rosendahl, 2000](#Rosendahl2000) and
[Holzer and Sommerfeld, 2008](#Holzer2008). 

*/

/**
## Plots

~~~gnuplot Low Reynolds number correlations for $\gamma = 6/1$
reset
set terminal svg font ",16"
set key top right spacing 1.1
set grid ytics
set xtics 0,1.e-1,10
set ytics 0,10,10000
set xlabel 'Re'
set ylabel 'C_D'
set xrange [1.e-2:10]
set yrange [1:10000]
set logscale

# Spheroid eccentricity (x is aspect ratio)
ecc(x) = (x <= 1.) ? sqrt(1. - x*x) : \
	 sqrt(1. - 1./(x*x));

# Spheroid surface correction (x is aspect ratio)
corr(x) = (x == 1) ? 1. :						\
	  ((x < 1) ? 0.5*(1. + (1. - ecc(x)*ecc(x))/ecc(x)*atan(ecc(x))) : \
	  0.5*(1. + x/ecc(x)*asin(ecc(x))));

# Sphericity (x is apsect ratio)
phi(x) = x**(2./3.)/corr(x);

# Crosswise sphericity for phi=0 (x is apsect ratio)
phiperpzero(x) = x**(2./3.);

# Crosswise sphericity for phi=90 (x is apsect ratio)
phiperpninety(x) = x**(-1./3.);

# Happel and Brenner's drag correction (x is aspect ratio)
K(x) = (x == 1.) ? 1. : \
       ((x < 1) ? 4./3./sqrt(x*x/(ecc(x)*ecc(x)) + 1.)/(x/ecc(x) - (x*x/(ecc(x)*ecc(x)) - 1.)*(pi/2. - atan(x/ecc(x)))) : \
       -4./3./sqrt(1./(ecc(x)*ecc(x)) - 1.)/(1./ecc(x) - (1./(ecc(x)*ecc(x)) + 1.)*(0.5*(log(1. + 1./ecc(x)) - log(1./ecc(x) - 1.)))));

# Stokes drag (x is aspect ratio, y is Reynolds number)
fStokes(x,y) = 24./y*K(x);

# Oseen's drag (x is aspect ratio, y is Reynolds number)
fOseen(x,y) = fStokes(x,y)*(1. + 3./16.*x*y);

# Breach's drag (x is aspect ratio, y is Reynolds number)
fBreach(x,y) = fStokes(x,y)*(1. + 3./16.*x*y + 9./160.*(x*y)*(x*y)*log(y/2.));

# Holzer and Sommerfeld, 2008 correlation for phi=0 (x is aspect ratio, y is Reynolds number)
fHSzero(x,y) = 1/y*(8./sqrt(phiperpzero(x)) + 16./sqrt(phi(x)) + 3.*sqrt(y)/(phi(x))**(3./4.)) + 0.42/phiperpzero(x)*(10.**(0.4*(-log(phi(x)))**(0.2)));

# Holzer and Sommerfeld, 2008 correlation for phi=90 (x is aspect ratio, y is Reynolds number)
fHSninety(x,y) = 1/y*(8./sqrt(phiperpninety(x)) + 16./sqrt(phi(x)) + 3.*sqrt(y)/(phi(x))**(3./4.)) + 0.42/phiperpninety(x)*(10.**(0.4*(-log(phi(x)))**(0.2)));

# Rosendahl, 2000 correlation (x is aspect ratio, y is Reynolds number, z is incidence angle)
fRosendahl(x,y,z) = fHSzero(x,y) + (fHSninety(x,y) - fHSzero(x,y))*(sin(z)**3.);

# Zastawny et al., 2012 correlation for x=5/2 (y is Reynolds number, z is incidence angle)
a0 = 2.; a1 = 5.1; a2 = 0.48; a3 = 15.52; a4 = 1.05;
a5 = 24.68; a6 = 0.98; a7 = 3.19; a8 = 0.21;
CDzero(y)  = a1/(y**a2) + a3/(y**a4);
CDninety(y) = a5/(y**a6) + a7/(y**a8);
fZastawny(y,z) = CDzero(y) + (CDninety(y) - CDzero(y))*(sin(z)**a0);

plot fStokes(6,x)      w l lw 2 lc rgb "blue"      t "Stokes",			\
     fOseen(6,x)       w l lw 2 lc rgb "red"       t "Oseen",			\
     fBreach(6,x)      w l lw 2 lc rgb "sea-green" t "Breach", \
     fRosendahl(6,x,0) w l lw 2 lc rgb "black"     t "Rosendahl, 2000",		\
     fHSzero(6,x)      w l lw 2 lc rgb "orange"    t "Holzer and Sommerfeld, 2008, phi=0", \
     fHSninety(6,x)    w l lw 2 lc rgb "brown"     t "Holzer and Sommerfeld, 2008, phi=90", \
     fZastawny(x,0)    w l lw 2 lc rgb "purple"    t "Zastawny et al., 2012"
~~~
*/

/**
## References

~~~bib
@article{Clift2005,
  title={Bubbles, drops, and particles},
  author={Clift, Roland and Grace, John R and Weber, Martin E},
  year={2005},
  publisher={Courier Corporation}
}

@book{Happel2012,
  title={Low {R}eynolds number hydrodynamics: with special applications to particulate media},
  author={Happel, John and Brenner, Howard},
  volume={1},
  year={2012},
  publisher={Springer Science \& Business Media}
}

@article{Masliyah1970,
  title={Numerical study of steady flow past spheroids},
  author={Masliyah, Jacob H and Epstein, Norman},
  journal={Journal of Fluid Mechanics},
  volume={44},
  number={3},
  pages={493--512},
  year={1970},
  publisher={Cambridge University Press}
}

@article{Holzer2008,
  title={New simple correlation formula for the drag coefficient of non-spherical particles},
  author={H{\"o}lzer, Andreas and Sommerfeld, Martin},
  journal={Powder Technology},
  volume={184},
  number={3},
  pages={361--365},
  year={2008},
  publisher={Elsevier}
}

@article{Rosendahl2000,
  title={Using a multi-parameter particle shape description to predict the motion of non-spherical particle shapes in swirling flow},
  author={Rosendahl, Lasse},
  journal={Applied Mathematical Modelling},
  volume={24},
  number={1},
  pages={11--25},
  year={2000},
  publisher={Elsevier}
}

@article{Schiller1933,
  title={Fundamental calculations in gravitational processing},
  author={Schiller, L and Naumann, A},
  journal={Z. Ver. Dtsch. Ing},
  volume={77},
  pages={318--320},
  year={1933}
}

@article{Zastawny2012,
  title={Derivation of drag and lift force and torque coefficients for non-spherical particles in flows},
  author={Zastawny, Marian and Mallouppas, George and Zhao, Fan and Van Wachem, Berend},
  journal={International Journal of Multiphase Flow},
  volume={39},
  pages={227--239},
  year={2012},
  publisher={Elsevier}
}
~~~
*/


