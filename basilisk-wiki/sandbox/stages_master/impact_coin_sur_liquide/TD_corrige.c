# Introduction

L'impact hydrodynamique, entre un objet solide et un fluide est très souvent modélisé et étudié dans le domaine naval, mais également dans le domaine aéronautique en terme de sécurité notamment pour l'amerrissage des avions ou hydravion. Nous nous intéressons à l'impact de carènes de bateaux au contact d'un fluide. La carène d'un navire constitue la partie immergée de celui-ci. De part la simulation numérique ainsi que l'expérience pratique, nous allons déterminer les champs aérodynamiques de notre système et visualiser l'écoulement.

## Modélisation du problème

Notre problème peut être, dans un premier temps modélisé de la façon suivante:

Avec :


$\rho$ la masse volumique de l'eau ($kg/m^3$) 

$\mu$, la viscosité dynamique de l'eau ($Pa.s$)

$\mathbf{u}$, la vitesse  du fluide ($m.s^{-1}$)
  
$U\mathbf{z}$, la vitesse  d'impact de la coque avec le fluide ($m.s^{-1}$)

$p_0$, la pression atmosphérique à la surface de l'eau ($Pa$)

$\alpha$, angle dièdre que forme la coque ($rad$)

$\mathbf{n}$, vecteur normale sortant à la coque


## Hypothèses

L'écoulement est très fortement instationnaire, de plus on considère que l'écoulement est à haut nombre de Reynolds donc les zones rotationnelles sont confinés dans des couches limites d'épaisseur très fines le longs des parois. En première approximation, on peut alors considérer l'écoulement comme irrotationnel. On a la relation suivante: $$ \nabla \wedge\mathbf{u}=\mathbf{0}$$
Or, on constate que quelque soit la fonction scalaire ${\phi}$, la relation $\nabla \wedge( \nabla{\phi})= \mathbf{0}$ est vérifiée.
Il existe alors une fonction ${\phi}$, appelée potentiel des vitesses, tel que $\mathbf{u}=\nabla{\phi}$.

Les efforts compressibles commencent à se manifester dès lors que les vitesses sont comparables à la vitesse du son (dans l'eau $c \sim 1500 m.s^{-1}$). On considère l'écoulement isovolume ce qui implique que le Mach est faible:

$Ma=\frac{U}{c}\ll1$

Nous faisons également l'hypothèse de fluide parfait.

## Equations régissant l'écoulement

Notre système est régi par les équations de Naviers-Stokes suivantes:

* Équation de continuité (équation de bilan de la masse)

$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0$.

Or, on considère l'écoulement comme étant incompressible, donc $div (\mathbf{u}) = 0\Leftrightarrow div(\nabla{\phi})=0$, soit: $\Delta{\phi}=0$.


* Équation de bilan de la quantité de mouvement

$\rho \frac{D\mathbf{u}}{D t} =-\frac{\partial p}{\partial x} + \rho \mathbf{g}  \Leftrightarrow \rho\frac{\partial u}{\partial t} + (\mathbf{u}. \nabla)\mathbf{u}=-\nabla p+ \rho \mathbf{g} $

En projetant selon les axes:

Selon $\mathbf{x}$ : $\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x}+w \frac{\partial u}{\partial z}=-\frac{1}{\rho}\frac{\partial p}{\partial x} $ (1)

Selon $\mathbf{z} $: $\frac{\partial w}{\partial t} + u \frac{\partial w}{\partial x}+w \frac{\partial w}{\partial z}=-\frac{1}{\rho}\frac{\partial p}{\partial z}+g $ (2), avec $g$  accélération de la pesanteur ($m/s^{-2}$)

En dérivant partiellement l'équation (1) par rapport à z et l'équation (2) par rapport à x et en soustrayant l'équation (1) à l'équation (2), on a:

$\frac{\partial \Omega_y}{\partial t} + u \frac{\partial  \Omega_y}{\partial x}+w \frac{\partial  \Omega_y}{\partial z}+\Omega_y div(\mathbf{u})=0 $ (3)

Avec $\Omega_y $ la composante selon $\mathbf{y} $ du vecteur vorticité $\vec{\Omega}$ défini tel que: $\vec{\Omega}={{\mathbf{rot}}}\ \mathbf{u}  $.
D'où $\Omega_y= \frac{\partial u}{\partial z}-\frac{\partial w}{\partial x}$

Or l'écoulement étant incompressible $ div(\mathbf{u})=0$, alors la relation (3) devient: 

$\frac{\partial \Omega_y}{\partial t} + u \frac{\partial  \Omega_y}{\partial x}+w \frac{\partial  \Omega_y}{\partial z}=0 $

On peut appliquer le théorème de bernoulli, en effet, l'écoulement est incompressible, irrationnel et le fluide parfait.

On a la relation vectorielle suivante:

$(\mathbf{u}.\nabla)\mathbf{u}=\vec{\Omega}\wedge\mathbf{u}+\frac{1}{2}\nabla(\mathbf{u}^2)$

Or: $\vec{\Omega}=\mathbf{0}$ et $\mathbf{u}=\nabla \phi$

Donc: $$(\frac{\partial \nabla \phi}{\partial t}+ \nabla\frac{1}{2}(\nabla \phi)^2 +\frac{1}{\rho}\nabla p+g\mathbf{z})=0$$

$$\nabla(\frac{\partial  \phi}{\partial t}+ \frac{1}{2}(\nabla \phi)^2 +\frac{p}{\rho}+gz)=0 \leftrightarrow \frac{\partial  \phi}{\partial t}+ \frac{1}{2}(\nabla \phi)^2 +\frac{p}{\rho}+gz =F(t)$$

## Conditions aux limites 

Les conditions aux limites associées à cet écoulement sont: 

* à la surface libre
$p=p_0$

* L'hypothèse de fluide parfait impliqué que le fluide peut glisser le long de la paroi mais l'imperméabilité de la coque implique que la composante normale de la vitesse du fluide à la paroi est égale à la composante de la vitesse de la coque soit:
$\mathbf{u}.\mathbf{n}=\mathbf{U}.\mathbf{n}$, le long de la coque, soit pour: $z=-tan(\frac{\pi-\alpha}{2})\mid x \mid +\, Ut $


* Condition de champ lointain
$\mathbf{u}\rightarrow\mathbf{0}$



# Ordre de grandeur et recherche de solution autosimilaire
## Orde de grandeur et loi d'échelle

L'ordre des termes intervenant dans l'équation de Bernoulli est:

$$\rho\frac{\partial \phi}{\partial t}\sim\rho\frac{UL}{T}  $$

$$\frac{1}{2}\rho(\nabla{\phi})^2\sim\rho U^2 $$

$$\rho g z \sim \rho g L$$

$$ p \sim P  $$

Dans ce problème les termes dominants sont le terme instationnnaire et la pression.

Lorsque $\rho\frac{UL}{T}\sim\rho g L  $ , donc que $T \sim\frac{U}{g}  $, le terme de gravité devient lui aussi important: c'est la fin du régime inertiel. Pendant cette durée, on peut négliger le terme inertiel.

A l'aide du théorème $\Pi $, nous allons proposer une loi d'échelle pour la force par unité de longueur $f  $ subie par la coque.

La force $f $ par unité de longueur dépend des paramètres:
$f=h(\rho,U,t,\alpha) $, avec $h$ fonction quelconque.

Chaque paramètre a pour dimension:

$$\left[\rho\right] =L^2T^{-1} \quad \left[x\right] =L \quad \left[y\right] =L  \quad \left[\rho\right] =ML^{-3}  \quad \left[U\right] =LT^{-1} \quad \left[\alpha\right] =1$$

On a donc une relation entre 4 paramètres plus 1 observable et 3 dimensions différentes pour caractériser le système. On peut donc exprimer $f $ en fonction de 2 grandeurs sans dimension.

Par analyse dimensionnelle, on trouve:

$$f=h(\rho U^3t,F(\alpha)) $$

A la fin du régime purement inertiel la force est maximale:

 $$ f_{max}=f(t\simeq\frac{U}{g}) \sim \frac{\rho U^4}{g} F(\alpha) $$
 
 Nous allons maintenant chercher une loi d'échelle pour le potentiel des vitesses <math> \phi </math>, faisant intervenir des variables d'espaces adimensionnées.

 
 Le potentiel  $\phi$ dépend de: $\phi=\Phi(x,y,\rho,U,t,\alpha)$, avec $\Phi$ fonction quelconque.
Chaque paramètre a pour dimension:

$$\left[f\right] =MT^{-2} \quad \left[\rho\right] =ML^{-3} \quad \left[t\right] =T  \quad \left[\alpha\right] =1  \quad \left[U\right] =LT^{-1}$$

On a donc une relation entre 6 paramètres plus 1 observable et 3 dimensions différentes pour caractériser le système. On peut donc exprimer $phi  $ en fonction de 4 grandeurs sans dimension.

Par analyse dimensionnelle, on trouve:

$$\pi=\Phi(\pi_1,\pi_2,\alpha),\quad  avec \quad  \pi=\frac{\phi}{U^2 t},\quad  \pi_1=\frac{x}{U t},\quad  \pi_2=\frac{y}{U^2 t}$$

Soit: $$\phi(x,y,t)=U^2 t \Phi  (\frac{x}{U t},\frac{y}{U t},\alpha),\quad \pi=\frac{\phi}{U^2 t}$$

## Adimensionnement

Nous allons maintenant adimensionner l'équation de Laplace $\Delta \phi$.

On note:  $\tilde{\phi}=\frac{\phi}{U^2 t}, \quad \tilde{x}=\frac{x}{U t},  \quad \tilde{y}=\frac{y}{U t} $

d'où:

 $$\frac{\partial \phi}{\partial x}=\frac{U^2 t}{U t} \frac{\partial \tilde{\phi}}{\partial \tilde{x}}=U\frac{\partial \tilde{\phi}}{\partial \tilde{x}}$$
<br />
 <center> et  </center>
 <br />

$$\frac{\partial ^2 \phi}{\partial x^2}=\frac{U^2 t}{U^2 t^2} \frac{\partial ^2 \tilde{\phi}}{\partial \tilde{x}^2}=\frac{1}{t} \frac{\partial ^2 \tilde{\phi}}{\partial \tilde{x}^2}$$

<br />
alors: $$\Delta \phi=\frac{1}{t}\left( \frac{\partial ^2 \tilde{\phi}}{\partial \tilde{x}^2}+ \frac{\partial ^2 \tilde{\phi}}{\partial \tilde{y}^2} \right)=0 \quad \leftrightarrow \quad \tilde{\Delta}\tilde{\phi}=0 $$

## Solution autosimilaire

Dans le domaine réel, la côte de la coque est défini par: $z=-tan(\frac{\pi-\alpha}{2})\mid x \mid +\, Ut $
<br />
Soit dans l'espace autosimilaire: 

$$\frac{z}{U t}=-tan(\frac{\pi-\alpha}{2})\frac{\mid x \mid}{ U t} +\frac{U t}{U t}\quad \leftrightarrow \quad \tilde{z}=-tan(\frac{\pi-\alpha}{2}) \mid \tilde{x} \mid +1 $$


<span style="color: #FF0000;">Le problème est devenu stationnaire dans l'espace autosimilaire.</span>

