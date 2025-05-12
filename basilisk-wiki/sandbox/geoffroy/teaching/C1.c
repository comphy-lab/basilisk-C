/**
# Prise en main du langage Basilisk 

## Basilisk, qu'est ce que c'est ?

  Basilisk est un code créé par
  [S. Popinet](http://gerris.dalembert.upmc.fr/user.html:Popinet)
  à l'[Institut Jean le Rond
  d'Alembert](http://www.dalembert.upmc.fr/ijlrda/) (Paris VI) et est
  développé par [différents
  chercheurs](http://basilisk.fr/src/AUTHORS) à travers le monde. On
  peut considérer Basilisk comme étant le petit frère de
  [Gerris](http://gerris.dalembert.upmc.fr/main_page.html) du
  même auteur. En effet, tout comme Gerris, Basilisk est un [logiciel
  libre](https://fr.wikipedia.org/wiki/Logiciel_libre) qui permet de
  résoudre des équations différentielles partielles sur un maillage
  cartesien adaptatif, et plus particulièrement les équations du
  mouvement des fluides (Navier-stokes, Saint-Venant, etc). 

## Installer Basilisk

Pour commencer, installez Basilisk en suivant [ce
lien](http://basilisk.fr/src/INSTALL). Veillez à bien faire
l'installation à l'aide de darcs et n'oubliez pas d'installer tous les
paquages additionels précisées en pied de page. Une fois que vous avez
compilé le code à l'aide de l'instructions "make", verifiez que tout
fonctionne à l'aide de la commande : 

~~~bash
qcc --version
~~~

Appelez votre professeur pour qu'il valide votre installation puis
passez à la suite.

## Inclure des fichiers sources

 Nous allons programmer le cas d'une [Vague qui déferle sur un plage
1D](http://basilisk.fr/sandbox/M1EMN/Exemples/slope.c). (Ce cas est 
reproduit d'après les cours de [Pierre-Yves Lagrée](http://www.lmm.jussieu.fr/~lagree/), avec
son aimable autorisation)

Pour cela, nous utiliserons comme modèle physique les équations de
Saint-Venant, également appelées "équations en eaux peu profondes". Ce
système d'équation s'écrit en 1D de la manière suivante :

$$
 \left\{\begin{array}{l}
\partial_t h+\partial_x Q=0\\ \partial_t Q+ \partial_x
\left[\dfrac{Q^2}{h}+g\dfrac{h^2}{2}\right] = - gh \partial_x Z
\end{array}\right.  
$$ 

Nous reconaissons l'équation de conservation de la masse ainsi que
celle de conservation de la quantité de mouvement. Nous remarquons
également que, dans cette modélisation, nous ne prenons pas en compte
les frottements. Ces équations, ainsi que la façon de les résoudre est
déjà programmée dans Basilisk et se trouve dans le fichier
"saint-venant.h". De plus, nous allons résoudre ce cas sur un maillage
purement cartesien (non-adaptatif). Commencez donc votre programme en
incluant ces deux fichiers :
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"

/**
Vous pouvez ouvrir ces fichiers dans votre éditeur afin de voir ce
qu'ils contiennent. Je vous déconseille pour l'instant d'allez voir ce
qu'il se trouve dans le fichier "grid/cartesian1D.h". Ouvrez plutôt le
fichier "saint-venant.h". Il se trouve, comme tout les autres fichiers
source de Basilisk, dans le répertoire :

~~~literatec
 ~/basilisk/src/
~~~
Vous êtes libre d'écrire dans ces fichiers donc faites attention ou il
vous faudra réinstaller Basilisk. Vous pouvez également lire ce fichier 
sur le net en suivant ce lien : [http://basilisk.fr/src/saint-venant.h](http://basilisk.fr/src/saint-venant.h).
Vous pouvez remarquer que le code est commenté et mis en page. Vous pouvez voir que le fichier
"saint-venant.h" est composé d'une première partie où sont déclarées les
variables, suivie par la déclaration des fonctions
"advance_saint_venant" et "update_saint_venant, et enfin une dernière
partie composée de différentes déclarations d'"event".

Exercice : Comment s'appelle la fonction qui appelle le solveur de Riemann
utilisée dans le fonction update_saint_venant ?  (répondez en
commentaire à la fin de votre programme à la suite de la question que
vous aurez copié/collé).

## Les champs dans Basilisk : scalar, vector et tensor

Dans Basilisk, il est aisé de déclarer des champs scalaires,
vectoriels ou tensoriels. Notons qu'un champs vectoriel en 1 dimension
est un champ scalaire, ce qui ne sera pas le cas en 2D. Nous déclarons
un nouveau champ scalaire dh ainsi que quelques variables :*/
double a, tmax;
scalar dh[];
/**

Exercice : Quels sont les champs scalaires déclarés dans "saint-venant.h" ? Même
question pour les champs vectoriels.

## La fonction "main"

Pour qu'un code Basilisk soit compilable, il faut obligatoirement la
présence d'un fonction main(), à l'instar d'un programme écrit en
C. Cette fonction renvoit un entier et n'a aucun argument.  A
l'intérieur de cette routine, nous définissons les différents
paramêtres de la simulation et appellons la fonction run() de la façon
suivante :*/
 
int main(){
  X0 = -15.;
  L0 = 60.;
  N = 1024; 
  G = 9.81;
  tmax = 50; // temps final
  a = 0.1; // amplitude de la vague
  run();
}
/** 
Les paramêtres L0 et X0 définissent la taille du domaine (L0) ainsi
que son origine (X0 en 1D, {X0,Y0} en 2D), ces paramêtres sont
déclarés de base dans tous les programmes Basilisk. N défini le nombre
de volumes finis (ou éléments ou dalles) utilisés pour découper le
domaine dans une dimension, nous reviendrons sur ce point lorsque nous
passerons à 2 dimensions. N est également déclaré de base dans
Basilisk et doit être obligatoirement une puissance de 2. A contrario,
le paramêtre G, qui fixe la valeur de l'accélération de la gravité,
est déclaré dans "saint-venant.h". La routine run() est déclaré dans
le fichier "predictor-corrector.h", elle sert à définir le schéma
d'avancé en temps (ici predictor corrector qui est un schéma d'ordre
2).

Exercice : Quelle est la valeur de G avant que nous la fixions à 9.81 ?

## Les conditions aux bords

Les conditions aux bords s'écrivent "à nue" derrière la fonction
main(), comme une déclaration globale de variable. Nous allons fixer
des conditions sur le bord gauche du domaine, la plage se trouvera à
droite. Nous fixons une condition de neumann(0) sur la vitesse u.n et
la hauteur d'eau h : */

u.n[left] = neumann(0);
h[left] = neumann(0);

/**
Nous remarquons la présence de la composante 'n' du champ vectoriel
"u", cette composante signifie que nous fixons la composante normale
au bord, soit u.x[]. La composante 't' signifiera que l'on fixe la
composante tangentielle (qui n'existe pas en 1D). Ces composantes (
'n' et 't') n'existent que lors de la déclaration des conditions aux
bords.

## Les events 

Les events permettent d'executer du code à des instants précis lorsque
le solveur tourne. Ils se placent après la fonction main et leur syntaxe
est la suivante :

~~~literatec
event nom_event ( cond_event ) { code_event } 
~~~ 

"nom_event" est le nom de l'event choisi par l'utilisateur. Attention
tout de même aux noms "init" et "default" qui sont des noms d'events
spéciaux. Tous les autres noms peuvents être utilisés, en minuscule si
possible, et en un seul mot obligatoirement. "cond_event" est la
condition qui permet de déclarer quand l'event se réalisera grâce aux
variable "i" et "t", par exemple "t = 2" ou "i = 100" sont des
conditions valides. "code_event" est à remplacer par le code que l'on
souhaite executer.

Par exemple, nous allons écrire l'event qui fixe les conditions
initiales sur la topographie (zb) la vitesse (u.x) et sur la hauteur
d'eau (h) :*/
event init (i = 0)
{
  foreach(){
    zb[] =   (x>10)*(x-10)/(25);
    u.x[]=a*exp(-(x+12)*(x+12)) ; // pulse de vitesse pour créer la vague
    h[]=fmax(1+u.x[]-zb[],0);
    }
  boundary({zb,h,u});
}
/**
Le nom de l'event choisi, "init", est bien spécifique : il dit au
solveur que c'est le premier event "i = 0" à devoir s'éxecuter, le nom
d'event "init" permet donc de définir des events prioritaires sur les
autres. Il existe un autre nom pour les events de priorité maximal :
"defaults". La variable globale "i" définie le nombre de pas de temps
dans le solveur. La condition "i = 0" dit donc au compilateur que cet
event s'execute avant l'avancé en temps du premier pas de
temps. Notez que nous pouvons également utiliser la variable "t" pour
fixer la condition d'execution d'un event. "t" est le temps à
l'intérieur du solveur (temps virtuel). Nous nous intéressons à la
commande foreach(){} ainsi qu'à son contenu dans le paragraphe
suivant.

Exercice : Ouvrir le fichier "~/basilisk/src/test/events.c" dans un
editeur texte ou sur le wiki : [http://basilisk.fr/src/test/events.c](http://basilisk.fr/src/test/events.c).
Indiquez à quel moment les events se déclenchent. (On
pourra l'executer en rajoutant l'option "-events" lors de la
compilation pour voir le déclenchement des events, voir la partie
compilation).

## La boucle foreach(){}

Dans l'event init défini ci-dessus, nous avons utilisé la boucle
foreach(){}. Cette boucle parcourt tous les éléments du domaine. Elle
fonctionne exactement de la même manière en 1D ou en 2D, sur une
grille cartésienne, multiple ou adaptative. C'est d'ailleurs là son
intérêt majeur. A l'intérieur de cette boucle, les valeurs des champs
sont indiqués par l'ajout après le nom du champ des crochet "[]". Les
variables x et y correspondent à la valeur de l'abscisse et de
l'ordonnée du centre de l'élément. La variable Delta vaut la taille de
l'élément (rappelons qu'en cartésian simple, tous les éléments ont la
même taille).

La boucle foreach() permet également de faire des opérations sur les
éléments qui se trouvent autour de celui où la boucle se trouve. Pour
cela, il faut placer à l'intérieur des crochet la position relative de
l'élément voulu par rapport à l'élément où la boucle se trouve (que
j'appelerai "élément actuel"). Cela peut paraître compliqué à première
vue, voyons cela avec un exemple : h[-1] est la valeur que prend h[] dans
l'élément à gauche de l'élément actuel, de même h[+1] est l'élément
qui se trouve juste à droite de l'élément actuel.

Exercice : Ajoutez au code un event qui se nomme dif_fin et qui se
déclenchera à chaque pas de temps. Dans cette event, vous calculerez
la dérivée seconde spatiale du champ h[] grâce à un [schéma différence finie
d'ordre 2 centré](https://fr.wikipedia.org/wiki/Méthode_des_différences_finies)
et vous stoquerez le résultat dans le champ dh[] que nous avons
déclaré au début de ce TP.

Pour finir notre premier programme, nous ajoutons un event qui permet
de sortir le profil de hauteur d'eau ainsi que la vitesse à chaque pas
de temps, et dans un format lisible par gnuplot :
*/
event plot (t<tmax;t+=0.1) {
    printf("set title 'Vague sur une plage 1D ----- t= %.1lf '\n"
	 "p[%g:%g][-0.1:1.5]  '-' u 1:($2+$4) t'free surface (m)' w l lt 3,"
	 "'' u 1:3 t'velocity (m/s)' w l lt 4,"
	 "'' u 1:4 t'topo' w l lt 1\n",t,X0+1,X0+L0);
    foreach()
    printf (" %g %g %g %g %g \n", x, h[], u.x[], zb[], t);
    printf ("e\n"
	   "pause 0.1 \n\n");
}

/**
## La compilation

Nous pouvons compiler à la main grâce au compilateur qcc qui a les
mêmes options que le compilateur gcc (+ quelques ajouts). Par exemple,
nous avons vu précédemment que l'option -events permet de garder une
trace de l'exécution des events dans la deuxième sortie du programme
(log). Compilez votre code grâce à la commande suivante dans votre
terminal :

~~~bash
qcc CM1.c -o CM1.x -lm -O2 -Wall
~~~

Si tout se passe bien, vous devriez pouvoir ensuite exécuter votre
code à l'aide de la commande suivante (en redirigeant la sortie
principale dans le fichiers "out" ):

~~~bash
./CM1.x > out
~~~

Vous pouvez voir le résultat de votre première simulation dans gnuplot
grâce à la commande : 

~~~bash 
load 'out' 
~~~

## L'utilisation de Make

La commande "make code.tst" permet de compiler, d'executer et de
rediriger les fichiers de sortie en une seule commande. Copiez le
fichier Makefile se trouvant dans ~/basilisk/src/ dans votre
répertoire courant, puis testez make en executant la commande :

~~~bash 
make CM1.tst 
~~~

Si tout c'est bien passé, vous devriez maintenant avoir un dossier
"CM1" dans lequel vous allez trouver un fichier "out" et un fichier
"log", qui sont les deux sorties de votre programme. Vous pouvez faire
un "load 'out'" dans gnuplot pour vous assurez que le fichier "out"
est le même que précédemment. Pour se convaincre de l'utilité de la
commande make, relancez la commande :

~~~bash 
make CM1.tst 
~~~

Que vous a répondu le terminal ?

En effet, la commande make teste si le programme a changé depuis la
dernière fois qu'il a été executé. Ce qui est très utile en pratique :
on ne compte plus le nombre d'étudiants se plaignant d'un "bug" dans
leur programme alors qu'ils ne l'avaient simplement pas compilé et
qu'ils executaient une version antérieure de leur programme (qui,
elle, était buggée). Dans le doute : passez toujours par la commande
make  : si rien n'a changé, vous ne perdrez pas le temps d'une nouvelle
compilation/execution !.

## Lier les sorties avec votre programme

La commande make permet également de sortir des graphiques gnuplot grâce à la commande "make nom/plots". Ajoutez le code suivant à la fin de votre programme en commentaire et en remplaçant les "+" par des tildes :"~" :

~~~bash
+++gnuplot Animation of the free surface.
reset
set xlabel 'X'
set ylabel 'Z'
set term gif animate
set output 'movie.gif'
load './out'
+++
~~~

Maintenant, testez la commande make /plots :

~~~bash 
make CM1/plots
~~~

Vous pouvez admirer le résultat en ouvrant le fichier movie.gif dans
un explorateur internet de votre choix.  Voici le résultat chez moi :

~~~gnuplot Animation of the free surface.
reset
set xlabel 'X'
set ylabel 'Z'
set term gif animate
set output 'movie.gif'
load './out'
~~~

  Notez que vous pouvez accéder à la page source de cette page wiki à l'aide du lien "raw page source" qui se trouve en bas du menu défilant à gauche. Nous aurons l'occasion de revenir sur ce point, mais il est important de noter que dans Basilisk, le wiki et le code ne font qu'un.

[Retour Sommaire Cours](http://basilisk.fr/sandbox/geoffroy/teaching/README)
*/
