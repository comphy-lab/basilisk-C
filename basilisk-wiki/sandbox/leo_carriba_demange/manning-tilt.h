/**
## Principe
L'équation scalaire de continuité et l'équation vectorielle de conservation de la quantité de mouvement peuvent se rassembler en une seule équation vectorielle décrivant l'évolution du vecteur de variables
conservatives :
$$
\displaystyle\frac{\partial}{\partial t} \left[\begin{array}{c} h \\ h\,u_x \\ h\,u_y \end{array}\right] = -\boldsymbol{\nabla}\cdot
\left[\begin{array}{cc} h\,u_x & h\,u_y \\ h\,u_x^2+ \frac{1}{2} g h^2 & h\,u_x\,u_y \\ h\,u_x\,u_y & h\,u_y^2 + \frac{1}{2} g h^2\end{array}\right]
\ -\ g\,h\left[\begin{array}{c} 0 \\ \partial_x z_b \\ \partial_y z_b \end{array}\right]
\ -\ \frac{1}{\rho}\left[\begin{array}{c} 0 \\ \tau_x \\ \tau_y \end{array}\right]$$

$$
\partial_t\mathbf{W}\qquad\quad = \quad\qquad-\boldsymbol{\nabla}\cdot\mathbf{F}(\mathbf{W})
\qquad\qquad +\qquad\qquad \mathbf{S_0}(\mathbf{W})
\qquad +\qquad \mathbf{S_f}(\mathbf{W})
$$
Par défaut, le solveur Saint-Venant de Basilisk prend en compte le terme source $\mathbf{S_0}$ dû à la bathymétrie, mais il faut ajouter le terme source dû aux frottements. On peut réarranger cette décomposition en 2 termes.

Supposons que l'on puisse décomposer la cote du fond $z_b$ en une composante de moyenne nulle (ou constante) spatialement, et une tendance linéaire (plan) :

$$z_b(x,y) = z^\ast_b(x,y) - I_x\,x - I_y\,y$$

où $z^\ast_b(x,y)$ est la cote du fond débarrassée de sa tendance régionale (c'est-à-dire "plate en moyenne"), et $I_x$ et $I_y$ les inclinaisons ("*tilt *") dans les directions $x$ et $y$ respectivement (pour le dire autrement, attention donc, si $z_b$ est un plan leurs définitions sont respectivement $I_x=-\partial_x z_b$ et $I_y=-\partial_y z_b$). Dans ce cas, les équations précédentes deviennent :
$$\partial_t\mathbf{W}\quad = \quad-\boldsymbol{\nabla}\cdot\mathbf{F}(\mathbf{W})
\ -\ g\,h\left[\begin{array}{c} 0 \\ \partial_x z^\ast_b \\ \partial_y z^\ast_b \end{array}\right]
\qquad+\underbrace{g\,h\left[\begin{array}{c} 0 \\ I_x \\ I_y \end{array}\right]
\ -\ \frac{1}{\rho}\left[\begin{array}{c} 0 \\ \tau_x \\ \tau_y \end{array}\right]}_{}$$
$$
\partial_t\mathbf{W}\quad = \quad-\boldsymbol{\nabla}\cdot\mathbf{F}(\mathbf{W})
\quad+\qquad\mathbf{S^\ast_0}(\mathbf{W})
\quad+\qquad\qquad \mathbf{S'_f}(\mathbf{W})\qquad\qquad
$$

On peut donc fournir en entrée à Basilisk la topographie $z^\ast_b$ sans tendance, et rajouter l'effet de l'inclinaison régionale dans le second terme source $\mathbf{S'_f}$ qui doit de toutes façons être
fourni pour prendre en compte le frottement.

##Implémentation avec le frottement Manning-Strickler

Dans le modèle de Manning-Strickler, le cisaillement au fond est donné par :
$$\boldsymbol{\tau}\quad=\quad\rho\,g\,n^2 h^{-\frac{1}{3}}\left|\mathbf{u} \right|\mathbf{u}\quad=\quad\rho\,g\,n^2 h^{-\frac{7}{3}}\left|\mathbf{q} \right|\mathbf{q}$$
Le schéma numérique résout d'abord l'équation d'évolution sans terme de frottement, puis calcule l'effet de ce terme (correction). Cela revient à résoudre l'équation de correction suivante :

$$
\displaystyle\frac{\partial}{\partial t} \left[\begin{array}{c}h \\ h\,u_x \\ h\,u_y \end{array}\right] =
\left[\begin{array}{c} 0 \\ +g\,h\,I_x - g\,n^2\,h^{-\frac{1}{3}}\left|\mathbf{u} \right|u_x \\
 +g\,h\,I_y - g\,n^2\,h^{-\frac{1}{3}}\left|\mathbf{u}\right|u_y
\end{array}\right]
$$

Sur cette étape de correction, la hauteur est supposée constante ($\partial_t h = 0$). Il reste donc :

$$
\displaystyle\frac{\partial}{\partial t} \left[\begin{array}{c}u_x \\ u_y \end{array}\right] =
\left[\begin{array}{c} +g\,I_x - g\,n^2\,h^{-\frac{4}{3}}\left|\mathbf{u} \right|u_x \\
 +g\,I_y - g\,n^2\,h^{-\frac{4}{3}}\left|\mathbf{u}\right|u_y
\end{array}\right]
$$

Cette équation peut être assez raide pour de fortes vitesses, on ne peut donc pas la traiter avec un schéma totalement explicite. Considérerons l'évolution de la composante $u_x$. Pour calculer son évolution du temps $t_k$ au temps $t_{k+1}=t_k+\Delta t$, on va linéariser le terme quadratique de la vitesse et l'exprimer implicitement :

$$\left|\mathbf{u}\right|\,u_x \simeq |\mathbf{u}^{(k)}|\,u_x^{(k+1)}$$

On obtient donc :

$$
u_x^{(k+1)}-u_x^{(k)} \simeq \Delta t \left[g\,I_x - g\,n^2\,|\mathbf{u}^{(k)}|h^{-\frac{4}{3}}u_x^{(k+1)}\right]$$
$$u_x^{(k+1)}\Big[1 + \underbrace{g\,n^2\,|\mathbf{u}^{(k)}|h^{-\frac{4}{3}}\Delta t}_{s}\Big] \simeq u_x^{(k)} + g\,I_x\,\Delta t
$$
On définit la quantité adimensionnelle :
$$
s = g\,n^2\,|\mathbf{u}^{(k)}|h^{-\frac{4}{3}}\Delta t $$
D'où finalement :
$$u_x^{(k+1)}=\frac{u_x^{(k)} + g\,I_x\,\Delta t}{1+s}$$
*/

// Manning coefficient
scalar nmanning[];
// Tilt
coord tilt = {.x = 0., .y = 0.};

event manningsourceterm(i++){
  
  foreach(){
    if(h[] > dry){      
      double s = dt*G*sq(nmanning[])*norm(u)/pow(h[],4./3.);
      foreach_dimension()
        u.x[] = (u.x[]+G*dt*tilt.x)/(1.+s);
    }
  }
  boundary ((scalar *){u});
}