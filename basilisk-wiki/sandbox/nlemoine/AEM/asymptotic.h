/**
Translation of $+n\Delta x$:
$$\zeta = \left(z-\frac{z_\mathrm{max}-z_\mathrm{min}}{2}-n\Delta x\right)\frac{2}{z_\mathrm{max}-z_\mathrm{min}}$$
Translation of $+n\Delta x + \frac{\Delta x}{2}$ of conjugate element:
$$\zeta' = \left(z-\frac{\overline{z}_\mathrm{max}-\overline{z}_\mathrm{min}}{2}-n\Delta x-\frac{\Delta x}{2}\right)\frac{2}{\overline{z}_\mathrm{max}-\overline{z}_\mathrm{min}}$$
Translation of $-n\Delta x$:
$$\zeta'' = \left(z-\frac{z_\mathrm{max}-z_\mathrm{min}}{2}+n\Delta x\right)\frac{2}{z_\mathrm{max}-z_\mathrm{min}}$$
Translation of $-n\Delta x + \frac{\Delta x}{2}$ of conjugate element:
$$\zeta''' = \left(z-\frac{\overline{z}_\mathrm{max}-\overline{z}_\mathrm{min}}{2}+n\Delta x-\frac{\Delta x}{2}\right)\frac{2}{\overline{z}_\mathrm{max}-\overline{z}_\mathrm{min}}$$

The residual (complex) potential which must be added to the one created by the truncated series is the sum of the potentials created by all missing quadruplets: 

$$\Omega_\mathrm{res} = 
\sum_{n>n_\mathrm{trunc}}\left(\textstyle\frac{Q}{2\pi}\log\zeta
-\frac{Q}{2\pi}\log\zeta'+\frac{Q}{2\pi}\log\zeta''-\frac{Q}{2\pi}\log\zeta'''\right)
=
\frac{Q}{2\pi}\sum_{n>n_\mathrm{trunc}}\log\left(\frac{\zeta\zeta''}{\zeta'\zeta'''}\right)$$

Denoting $z_\mathrm{mid}=\frac{z_\mathrm{max}+z_\mathrm{min}}{2}$ and $(z_\mathrm{max}-z_\mathrm{min}) = R e^{i\alpha}$ we have:

$$\Omega_\mathrm{res}(z)\quad=\quad \frac{Q}{2\pi}\sum_{n>n_\mathrm{trunc}}\log\left(
\frac{(z-z_\mathrm{mid}-n\Delta x)(z-z_\mathrm{mid}+n\Delta x)}
{\left(z-\overline{z}_\mathrm{mid}-\frac{\Delta x}{2}-n\Delta x)(z-\overline{z}_\mathrm{mid}-\frac{\Delta x}{2}+n\Delta x\right)}
\right)
\quad+\quad
\frac{Q}{\pi}\sum_{n>n_\mathrm{trunc}}\log\left(\frac{\overline{z}_\mathrm{max}-\overline{z}_\mathrm{min}}{z_\mathrm{max}-z_\mathrm{min}}\right)
$$

$$\Omega_\mathrm{res}(z)\quad=\quad \frac{Q}{2\pi}\sum_{n>n_\mathrm{trunc}}\log
\left(
\frac{1-\frac{(z-z_\mathrm{mid})^2}{n^2\Delta x^2}}
{1-\frac{(z-\overline{z}_\mathrm{mid}-\frac{\Delta x}{2})^2}{n^2\Delta x^2}}
\right)
\quad+\quad
\frac{Q}{\pi}\sum_{n>n_\mathrm{trunc}}-2i\alpha$$

Let denote

$$u=\frac{z-z_\mathrm{mid}}{\Delta x}, \qquad u'=\frac{z-\overline{z}_\mathrm{mid}}{\Delta x}-\frac{1}{2}$$

then

$$\Omega_\mathrm{res}(z)\quad=\quad \frac{Q}{2\pi}\sum_{n>n_\mathrm{trunc}}\left[\log\left(1-\left(\textstyle\frac{u}{n}\right)^2\right)
-\log\left(1-\left(\textstyle\frac{u'}{n}\right)^2\right)
\right]
\quad+\quad
\frac{Q}{\pi}\sum_{n>n_\mathrm{trunc}}-2i\alpha$$

If the truncation index $n_\mathrm{trunc}$ is large enough, we can use a Taylor expansion for the logarithms which are in the form $\log(1-\varepsilon)$ with $\left|\varepsilon\right| \ll 1$, hence:

$$\log\left(1-\left(\textstyle\frac{u}{n}\right)^2\right) = -\left[
\left(\frac{u}{n}\right)^2
+ \frac{1}{2}\left(\frac{u}{n}\right)^4
+ \frac{1}{3}\left(\frac{u}{n}\right)^6
+ \cdots 
\right] = -\sum_{k=1}^\infty \frac{1}{k}\left(\frac{u}{n}\right)^{2k}$$

The residual potential is then:

$$\Omega_\mathrm{res}(z) \underset{n_\mathrm{trunc}\mathrm{\ large}}{\qquad\sim\qquad} \frac{Q}{2\pi}\sum_{n>n_\mathrm{trunc}}\left[
\sum_{k=1}^\infty \frac{1}{k}\left(\frac{u'}{n}\right)^{2k}
-\sum_{k=1}^\infty \frac{1}{k}\left(\frac{u}{n}\right)^{2k}
\right]
\quad+\quad
\frac{Q}{\pi}\sum_{n>n_\mathrm{trunc}}-2i\alpha$$

$$\Omega_\mathrm{res}(z) \underset{n_\mathrm{trunc}\mathrm{\ large}}{\qquad\sim\qquad} \frac{Q}{2\pi}\sum_{n>n_\mathrm{trunc}}\left[
\sum_{k=1}^\infty \frac{1}{n^{2k}}\frac{\left(u'^{\,2k}-u^{2k}\right)}{k}
\right]
\quad+\quad
\frac{Q}{\pi}\sum_{n>n_\mathrm{trunc}}-2i\alpha$$

Switching the order of the two sums yields:

$$\Omega_\mathrm{res}(z) \underset{n_\mathrm{trunc}\mathrm{\ large}}{\qquad\sim\qquad} \frac{Q}{2\pi}\sum_{k=1}^\infty 
\left[\frac{\left(u'^{\,2k}-u^{2k}\right)}{k}
\sum_{n>n_\mathrm{trunc}}
\frac{1}{n^{2k}}
\right]
\quad+\quad
\frac{Q}{\pi}\sum_{n>n_\mathrm{trunc}}-2i\alpha$$

Using the Riemann $\zeta$ function we define the correcting function:
$$F_\mathrm{corr}(s\,;n_\mathrm{trunc}) = \sum_{n>n_\mathrm{trunc}}\frac{1}{n^s} = 
\left[\sum_{n=1}^\infty\frac{1}{n^s} - \sum_{n=1}^{n_\mathrm{trunc}}\frac{1}{n^s}\right]
= \zeta(s) - \sum_{n=1}^{n_\mathrm{trunc}}\frac{1}{n^s}
$$
then
$$\Omega_\mathrm{res}(z) \underset{n_\mathrm{trunc}\mathrm{\ large}}{\qquad\sim\qquad} \frac{Q}{2\pi}\sum_{k=1}^\infty 
\left[F_\mathrm{corr}(2k\,;n_\mathrm{trunc})\frac{\left(u'^{\,2k}-u^{2k}\right)}{k}
\right]
\quad+\quad
\frac{Q}{\pi}\sum_{n>n_\mathrm{trunc}}-2i\alpha \qquad(1)$$

The last term does not have a limit, but since it is a pure imaginary number that is independent of $z$ we can evaluate the real part of the residual potential:

$$\Phi_\mathrm{res} = \mathcal{Re}\left\{\Omega_\mathrm{res}\right\}$$

as well as the residual velocity component:

$$v_\mathrm{res} = -\overline{\left(\frac{d\Omega_\mathrm{res}}{dz}\right)}$$

The advantage of expression (1) is that instead of using a large number of repetitions of the periodic slit pattern, we can use a moderate truncation order $n_\mathrm{trunc}$ (e.g. about 10) and then use a finite sum on $k$ (a second or third order expansion is just enough).
