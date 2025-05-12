%{

# Page de Alexandre Nazarians et Tony Spriet

## Projets développés

* We studied a wave propagation on a string which lenght is L = 2 pi.
However we have defined different boundaries conditions from the original code : [../vibrating_string.m]()

One sinusoidal excitation on x = 0, and one straight point on the other end x= L
This can be translated mathematically by these expressions : 
 
$$f(x=0,t)= sin(\frac{t}{2\pi})$$
$$f(x=L,t)= 0$$

Our code is available here : [../vibrating string forced.m]()

* We studied St-Venant equations which describe wave dynamics when the wavelenght is larger than the depth of the sea.
The code is based on [../wave_like.m]().
$$
\begin{array}{l}
u_t=-g\eta_x-bu\\
\eta_t=-(Hu)_x\\
\end{array}
$$
In our study, we have only considered 1D, infinite environment and a null fluid viscosity. 

* We created a Gif animation on our code [../free_surface_gravity_particles_gif.m]() from the original code [../free_surface_gravity_particles.m]()

![Free Gravity Surface particles](../freegravitysurfaceparticles.gif)


# Notes

De la part de Jérôme

~~~matlab
domaine	        valeur	note
connectivité 	2	1
recyclage	2	2
graphiques	2	1
théories	4	0
Originalité	4	1
note /14	14	5
~~~

Le code de corde vibrante avec forçage est interessant, mais c'est votre seule contribution originale. Pour StVenant, le lien n'est même pas sur votre page de contribution, il a fallu que j'aille le chercher dans les suggestions de contributions de "wave like". Ca a l'air de marcher, mais vous n'expliquez rien et il n'y a même pas de figure sur le wiki. Dommage...
}%