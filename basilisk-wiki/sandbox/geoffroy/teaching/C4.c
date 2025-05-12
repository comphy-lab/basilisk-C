/**
#  Les coefficients de friction dans l'équation de Saint-Venant et leur implémentation

## Les coefficients de friction

Commencez par coder les trois coefficients de friction
vu dans le cours précédant le TP : les coefficients de Manning, de Darcy-Weisbach ainsi que
celui de Poiseuille. Chacun de ces coefficients sera codé dans un
fichier séparé qui portera son nom. Par exemple, le code permettant de
prendre en compte le coefficient de friction de Poiseuille de manière
explicite sera codé dans le fichier "poiseuille-expl.h".

## Pluie sur un plan incliné

Nous allons tester les différents coefficients de friction ainsi que
leur implémentation (explicite ou implicite) en simulant un épisode de
pluie sur un plan incliné. Les résultats de nos simulations seront
comparés aux résultats expérimentaux de Frédéric Darboux : il a
réalisé des expériences de pluie sur plan incliné à l'INRA d'Orléans.

### Les Includes

Vous allez utiliser le solveur des équations de Saint-Venant sur une
grille 1D. Ajoutez également un des termes de friction que vous avez
codé. Il n'y aura qu'à changer le nom pour changer de terme de
friction (ou pour changer son implémentation).

### Fonction main()

A ce stade, vous devez être en mesure de savoir ce que vous devez
mettre dans la fonction main. Nous prendrons comme viscosité
cinématique de l'eau : $\nu = 10^{-6}$, comme coefficient de
Darcy-Weisbach : $f = 0.5$ et comme coefficient de Manning : $n =
0.025 s.m^{−1/3}$. La longueur de la pente est de 6.04 m avec X0 = -1.
L'expérience dure 1000 secondes.

### Terme source de pluie

Ajoutez un event qui ajoute de la pluie grâce au champs P[]. Vous
mettrez dans ce champ la valeur de l'intensité de la pluie. L'épisode pluvieux
dure 600 secondes.

### Conditions initiales

La pente est de 2% et l'intensité de la pluie est de 25 mm / h (attention aux unités
!). La pente est initiallement sèche, son point le plus haut se situera à gauche. 
Il ne pleut qu'entre $X = 0 m$ et $X = 4.1 m$. Ajoutez une chute d'eau en $X = 4.1 m$
d'environ 1 m de profondeur, elle représente le bac qui récupère l'eau dans le cas réel.

### Conditions aux bords

Fixez les bonnes conditions aux bords. Rappelez-vous qu'à droite nous avons un bac retenant l'eau.

### Output 

Ajoutez un event qui permet de calculer le débit en $g.s^{-1}$ à l'éxutoire (X = 4.04 m) pour une
largeur de canal de 11,5 cm. Imprimez le débit et le temps dans un
fichier qui sera lu par gnuplot. 

### Animation 

Ajoutez l'event suivant qui va vous permettre de voir le résultat de
votre simulation à l'aide de gnuplot :*/

event initplot(t = 0) {
   printf("set xlabel 'X'\n");
  }
event plot ( t <= tmax; t += 0.5) {
  printf("set title 'Ideal rain ----- t= %.3g '\n"
	 "p[%g:%g][0:2]  '-' u 1:($2*1e3) t'free surface (mm)' w linespoints lt 3,"
	 "'' u 1:($3*10) t'velocity (10 cm/s)' w l lt 4\n",t,X0+1,X0+L0-1);
  foreach(){
    printf ("%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
  } 
  printf ("e\n"
	  "pause %.4g \n\n",pause);
}
/**
### Résultats 

Sur le même graphique, comparez le débit à l'éxutoire pour chaque
terme de friction traité implicitement avec les données expérimentales
de F. Darboux. Quel terme de friction permet de reproduire le mieux l'hydrographe du cas réel ?
Conclure sur le type d'écoulement.

Faire une animation du cas.
