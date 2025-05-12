/**
#  Mini-projet : Pluie idéale sur un plan incliné

## Télécharger b-flood

b-flood est une partie de code basilisk spécialement développée pour faire de la simulation d'hydraulique (inondation, crue, etc...)

Télécharger le code b-flood et copier-le dans le dossier /basilisk/src/ (si possible) : [http://basilisk.fr/sandbox/b-flood/Readme](http://basilisk.fr/sandbox/b-flood/Readme)

## Pluie sur un plan incliné

Ce cas test sera utilisé dans le cadre du projet menant à une conférence scientifique. 

Nous nous concentrons ici sur une cas idéal d'une pluie homogène tombant sur une surface plane inclinée et imperméable. Le même montage expérimental a été utilisé avant d'évaluer la validité des schémas numériques. 

La topographie plate est inclinée d'un angle a et une intensité de pluie constante égale à I (mm.h 1 ) est imposée. Le canal a une longueur L = 4,04 m (direction x) et largeur l = 11,5 cm (direction y), et est initialement sèche. La pluie entraîne un écoulement qui se caractérise par la profondeur de l'eau et le profil de vitesse, et enfin S 0 = tan(a) est la valeur absolue de la pente du canal. 

### Création du script

Nommer votre script, "ideal_rain_man_R25_P2.c", le mettre dans un dossier portant exactement le même nom sans le .c : "ideal_rain_man_R25_P2" . 
Tous les fichiers qui seront créés par l'éxécution de votre code seront également dans ce dossier.


### Les Includes

Vous allez utiliser le solveur des équations de Saint-Venant sur une
grille 1D. Faire un include de la grille 1D (voir TP1)


Nous allons utiliser un code Saint-Venant spécial hydraulique, 
nous allons sprendre en compte un coefficient de frottement de type Manning et nous allons faire pleuvoir sur notre topo. Pour cela, copier/coller les includes suivants : 

COpier/coller depuis le dossier b-flood les fichiers saint-venant-topo.h darcy.h manning.h et rain.h dans le dossier où vous travaillez*/

#include "./saint-venant-topo.h"
#include "./manning.h"
#include "./rain.h"
/**

Pour comprendre comment imposer une pluie : ouvrir [b-flood/rain.h](http://basilisk.fr/sandbox/b-flood/rain.h)

Pour comprendre comment imposer une friction de manning : ouvrir [b-flood/manning.h](http://basilisk.fr/sandbox/b-flood/manning.h)

### Fonction main()

A ce stade, vous devez être en mesure de savoir ce que vous devez
mettre dans la fonction main. Nous prendrons comme coefficient de Manning : $n =
0.025 s.m^{−1/3}$. La longueur de la pente est de 6.04 m avec X0 = -1.
L'expérience dure 1000 secondes.

### Terme source de pluie

Ajoutez un event qui ajoute de la pluie grâce au champs scalar rain[]. Faire bien attention à ce qu'il ne pleuve UNIQUEMENT sur la zone choisie ! (entre X = 0 et X = 4.04 m ).
Vous mettrez dans ce champ la valeur de l'intensité de la pluie. L'épisode pluvieux
dure 600 secondes. 
S'arranger pour que la pluie s'arrête à cette date.

### Conditions initiales

La pente est de 2% et l'intensité de la pluie est de 25 mm / h (attention aux unités !). La pente est initiallement sèche, son point le plus haut se situera à gauche.  Il ne pleut qu'entre $X = 0 m$ et $X = 4.1 m$. La pente continuera pendant 1 m derrière pour ne pas perturber l'écoulement

### Conditions aux bords

A droite, l'eau  doit sortir librement : ajouter les conditions aux bords suivants :
*/
u.n[right] = neumann(0.);
/**
A gauche, nous avons ajouté 1 m de pente qui ne sert à rien pour ne pas avoir à nous soucier des conditions aux bords de ce coté.


### Débit

(Warning : maths ici !)

Ajoutez un event qui permet de calculer le débit en $g.s^{-1}$ à l'éxutoire (X = 4.04 m) pour une largeur de canal de 11,5 cm.  Appeler ce débit "Qm".
On appelle débit le débit massique d'eau passant au travers d'une surface S : [voir la page wiki](https://fr.wikipedia.org/wiki/D%C3%A9bit_(physique)#D%C3%A9bit_volumique)
Pour passer du débit volumique au débit massique, il faut multiplier par la masse volumique de l'eau : $1000 kg / m^3$. 


Dans cet event, copier/coller ce code : */
Point point = locate ( 4.04 , 0 );
/**
Cette ligne "place" le code au point de coordonnée x = 4,04 et y= 0 (nous sommes en 1D...). Vous pouvez utiliser "h[]" et "u.x[]" dans tout l'event comme étant la valeur des champs scalaires uniquement en ce point.

Le calcul sera fait à chaque pas de temps.

Vous imprimerez ce débit dans un fichier de sortie à l'aide du code suivant : 
*/
 static FILE * fexu = fopen("debit_exut.dat","w");
 fprintf(fexu,"%g \t %g \n",t, Qm);
/**



### Animation de la surface libre

Ajoutez l'event suivant qui va vous permettre de voir le résultat de
votre simulation à l'aide de gnuplot :*/
double pause1 = 0.05; // Valeur de la pause entre chaque image en seconde
event initplot(t = 0) {
   fprintf(stderr,"set xlabel 'X'\n");
  }
event plot ( t <= tmax; t += 0.5) {
 fprintf(stderr,"set title 'Ideal rain ----- t= %.2g '\n"
	 "p[%g:%g][0:2]  '-' u 1:($2*1e3) t'free surface (mm)' w linespoints lt 3,"
	 "'' u 1:($3*10) t'velocity (10 cm/s)' w l lt 4\n",t,X0+1,X0+L0-1);
  foreach(){
    fprintf (stderr,"%g %g %g %g %g\n", x, h[], u.x[], zb[], t);
  } 
  fprintf (stderr,"e\n"
                  "pause %.4g \n\n",pause1);
}
/**
### Profile de hauteur d'eau et de vitesse

Créer un event qui enregistre dans un fichier les profiles de hauteur d'eau et de vitesse en régime permanent lorsqu'il pleut. Pour cela, regarder l'animation et choisir un temps inférieur à 600 secondes où le régime permanent est établi.

Vous écrirez dans un fichier texte nommé "profil.dat" : la variable x, la variable h[], la variable u.x[] et la variable zb[].

Inspirez-vous des events précédents.


## Plusieurs pentes et pluies                                
  
Nous souhaitons étudier les trois cas suivants, avec à chaque fois la vitesse de pluie et la pente :

* 25 mm.h 1 and 2%,
* 25 mm.h 1 and 5%,
* 50 mm.h 1 and 2%.

Nous venons de faire le premier cas : faire les deux suivants.

_ Pour 25 mm.h 1 and 5% : Nommer votre script, "ideal_rain_man_R25_P5.c", le mettre dans un dossier portant exactement le même nom sans le .c : "ideal_rain_man_R25_P5" 

_ Pour 50 mm.h 1 and 2% : Nommer votre script, "ideal_rain_man_R50_P2.c", le mettre dans un dossier portant exactement le même nom sans le .c : "ideal_rain_man_R50_P2" 

## Différents coefficients de friction 

Refaire les trois cas précédents mais avec le coef de Darcy, vous prendrez comme valeur pour le coefficient de Darcy : f=0.5.
Pour cela, vous remplacerez l'include de manning.h par : */
#include "b-flood/darcy.h"
/**

et vous fixerez la valeur de f dans la fonction main.

## Les livrables 

Pour chaque cas, vous devrez livrer le fichier du débit à l'éxutoire en fonction du temps, le fichier "profile.dat" contenant les profiles de hauteur d'eau et de vitesse en régime permanent ainsi que le fichier permettant de voir l'animation gnuplot. 

## Friction de Poiseuille

Créer le fichier : "ideal_rain_pois_R25_P2.c" et le dossier dans lequel il va.
Dans le cas d'une friction parfaite dans un écoulement de Poiseuille, nous sommes capables de calculer un terme de friction analytique. Ce terme modifie la vitesse de l'écoulement de la manière suivante

$$ 
u.x = \frac{u.x}{1+dt*3.*nu/(h*h);}, 
$$ 
avec $nu$ la viscosité dynamique de l'eau, qui est égale à 1e-6. (dtmax est déjà connu dans le code).


Coder un event qui réalise cette modification du champ u.x. Cet event se déclenchera à chaque pas de temps.

Refaire les calculs pour les trois cas à étudier (les mêmes que pour les autres frictions).

## Résultats

Pour chaque couple pente/pluie, vous ferez un graph du profile de hauteur avec les trois coefficients de friction. Les trois profiles devront apparaitre dans le même graph : utiliser la commande "replot" de gnuplot.


*/