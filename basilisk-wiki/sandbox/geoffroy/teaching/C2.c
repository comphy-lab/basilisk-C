/**
# Inputs, Outputs et Maillage adaptatif 

## Introduction 

L'idée de ce TP est de vous familiariser avec les commandes d'input et
d'output de Basilisk ainsi que de vous initier au maillage adaptatif
intelligent utilisé par Basilisk. Pour cela, nous allons simuler
l'apparition d'un Tsunami au large de Tunis et nous allons étudier son
impact sur les côtes françaises. Commencez par créer un fichier
tsunami_votrenom.c dans lequel vous allez programmer la suite du TP. 

## Les Inputs

### Les include 

Nous commençons par les "include" de fichiers sources. Commencez par
inclure le fichier "spherical.h". Ce fichier spécifie que nous allons
utiliser des coordonnées sphériques. Ensuite, vous allez inclure le
fichier "saint-venant.h" qui permet de simuler des écoulements dont la
longueur caractéristique horizontale est bien plus grande que la
hauteur caractéristique. En effet, bien que la hauteur d'un Tsunami
nous paraisse impressionnante lorsque celui-ci touche nos côtes, cette
hauteur est extrémement petite comparée aux centaines de kilomêtres
d'étendue de la [mer
méditérannée](https://fr.wikipedia.org/wiki/Mer_Méditerranée).  Vous
allez également inclure le fichier "okada.h" qui permet de reproduire
l'effet d'un tremblement de terre sur la surface de la mer, ainsi que
le fichier "terrain.h". Ce dernier fichier permet d'importer dans
votre programme la topographie de la méditérannée au format Raster,
comme nous le verrons plus tard. Finissez la série d'include en
incluant le fichier "input.h". Donnez une valeur entre 2 et 12 à la
variable MAXLEVEL à l'aide la macro "#define". Nous la fixerons plus
tard précisement.


### La fonction main()

L'utilisation de coordonnées sphériques nécessite de préciser le rayon
de la terre dans la variable "Radius" : fixez cette valeur à
"6371220.". Fixez également la taille du domaine à 20 degrés 
carrés  ainsi que l'origine de la carte à l'aide de l'instruction :

~~~literatec
size(20.)   
origin (10. - L0/2., 36. - L0/2.);
~~~

Fixez G correctement. Notez qu'il s'agit là de la seule variable
dimensionnée. Nous finissons la fonction main() en fixant le nombre
d'éléments sur un axe "N" à la valeur de 2 à la puissance MAXLEVEL et
en appelant la routine run().

Exercice : Quelle est l'unité de temps dans le solveur ? Comment faire
pour utiliser comme unité de temps la minute ? Comment faire pour fixer 
comme unité de longueur le kilomêtre ?(ne pas le faire)

### Les conditions Initiales

Nous fixons les conditions initiales à l'aide d'un event que nous
créons. N'oubliez pas de lui donner un nom qui permettra de l'éxecuter
en priorité. Fixez sa condition d'execution de tel manière qu'il soit
éxecuté avant l'avancé en temps du solveur et une seule fois. Dans
cette event, vous allez utiliser la fonction terrain() afin d'importer
la topographie du bassin méditéranéen de la manière suivante : 

~~~literatec
terrain (zb, "./topo", NULL);
conserve_elevation();
~~~

Pour que cela fonctionne : téléchargez les fichiers 
["topo.kdt"](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/TopoMedi/topo.kdt),
["topo.sum"](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/TopoMedi/topo.sum) et 
["topo.pts"](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/TopoMedi/topo.pts).
Nous verrons lors
du mini-projet comment construire ces fichiers.

A l'aide d'une boucle foreach, définissez le niveau d'eau à l'altitude
zéro afin de remplir la mer méditérannée. Puis, ajoutez le bout de
code suivant qui permet de fixer la condition initiale d'un Tsunami au
large de Tunis :

~~~literatec
fault (x = 11.9, y = 37.26,
	 depth = 11e3,
	 strike = 323, dip = 40, rake = 5,
	 length = 150e3, width = 90e3,
	 U = 90,
	 iterate = adapt);
~~~

Exercice : Pourquoi avoir choisi cet endroit comme départ d'un Tsunami ?

Vous pouvez fermer l'event.

### Les conditions aux bords

Les conditions aux bords sont primordiales dans la résolution
d'équations aux dérivées partielles. Pour fixer la hauteur d'eau aux
bords à l'altitude zéro tout en la laissant libre de fluctuer autour
de cette valeur, nous allons utiliser la fonction radiation() : 

~~~literatec
u.n[left]   = - radiation(0);
u.n[right]  = + radiation(0);
u.n[bottom] = - radiation(0);
~~~

A présent, vérifiez que vous n'avez aucun bug en compilant votre
programme sans l'éxecuter en utilisant le fichier objet "kdt.o" qui se
trouve dans le dossier "~/basilisk/src/kdt". Appelez votre professeur
pour qu'il valide cette étape.

## Les Outputs

Nous allons maintenant nous intéresser aux différents fichiers de sortie. 

### Fichier "log"

Nous allons imprimer les variables relatives à la simulation dans le
deuxième fichier de sortie standard. Pour cela, créez un event qui
s'éxecutera à chaque pas de temps. A l'intérieur de cet event, nous
allons calculer différentes statistiques sur les champs de hauteur
d'eau et de vitesse à l'aide des deux commandes suivantes :

~~~literatec
stats s = statsf (h);
norm n = normf (u.x);
~~~

Ensuite, vous allez imprimez à chaque pas de temps les variables
suivantes dans le deuxième fichier de sortie à l'aide de la commande 
fprintf(stderr,...) : le temps virtuel, le pas de temps, le minimum de h, le maximum de h, le
maximum de u.x et la moyenne de u.x. Pour cela, regardez dans le
fichier "utils.h" à partir de la ligne 126 ce que contiennent les
structures "stats" et "norm", ainsi que comment elles sont calculées dans
les fonctions "statsf" et "normf" respectivement. Cette impression se
fera à l'aide de la fonction "fprintf(stderr,....)". Quand vous avez
fini, fermez l'event et compilez pour vérifier qu'il n'y a aucun bug.
Voici comment faire (copier/collé le code suivant):

~~~literatec
event logfile (i += 10)
{
  stats s = statsf (h);
  scalar v[];
  foreach()
    v[] = norm(u);
  
  stats no = statsf (v);
  if( i == 0 )
    fprintf (stderr, "#1t 2treal 3i 4h.min 5h.max 6h.sum 7u.min 8u.max 9dt\n");  
  fprintf (stderr, "%g  %d %g %g %g %g %g %g\n",
	   t, i, s.min, s.max, s.sum,  no.min, no.max, dt);
}
~~~
### Films et images.

Commençons par prendre une photo de la topographie. Pour cela, créez
un event qui s'éxecute une seule fois au début du code. Dans cette
event, déclarez une variable locale "zbmin" que vous allez initialiser
à zéro. Faites ensuite une boucle foreach() dans laquelle vous allez
tester toutes les valeur de zb[], lorsque la valeur de ce champ sera
inférieur à zbmin, fixer zbmin avec cette valeur : zbmin sera ainsi le
minimum du champ zb. Toujours dans le même event, prenez une photo de
la topographie en remplaçant "zbmax" par une "bonne" valeur dans la
commande suivante :

~~~literatec
static FILE * fzb = fopen ("topo.ppm", "w");
output_ppm (zb, fzb, min = zbmin, max = zbmax, n = 1 << MAXLEVEL, linear = true);
~~~

fermez l'event. Rajoutez l'évent suivant qui vous créera un film
montrant l'évolution de la surface libre  :

~~~literatec
event movies (t+=50) {
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  output_ppm (etam, mask = m, min = -2, max = 2, n = 512, linear = true,
	      file = "eta.mp4");
              }
~~~

Exercice : Quand cet event va-t-il s'éxecuter ?  Que représentent les couleurs rouge et bleu sur le film ? (répondre à cette question à la fin du TP)

Dans le même event, ecrire le code qui permet de faire un film de la variable "level". (cette variable est définie dans chaque cellule, comment faire pour la transformer en "scalar []" ?)

### Stations de mesure (Tide gauges)

Nous fixons plusieurs stations de mesures afin de mesurer l'évolution
de la hauteur d'eau en plusieurs points. Pour cela, nous créons la
fonction suivante :

~~~literatec
Gauge gauges[] = {
  // file   lon      lat 
  {"Nice.ga", 7.21,  43.65},
  {"Montpellier.ga", 3.93,  43.53},
  {"Perpignan.ga", 3.04,  42.70},
  {"Cannes.ga", 7.03,    43.54},
  {"Marseille.ga", 5.36,    43.30},
  {"Calvi.ga", 8.77,    42.57},
  {"Bastia.ga", 9.5,    42.70},
  {"Livorno.ga", 10.29,    43.54},
  {"Napples.ga", 14.27,    40.83},
  {"Tunis.ga", 10.34,    36.84},
  {NULL}
};
~~~

   Nous appelons cette fonction à l'aide de l'event suivant :
   
~~~literatec
event gauges1 (i++) output_gauges (gauges, {eta});
~~~

Pour finir cette partie, créez un event qui s'éxecutera au bout de 6
heures virtuelles. Cet event servira à donner un temps final au
solveur. A l'intérieur de cet event, faites simplement un printf("fin
du programme\n").

## Grille adaptative

Pour terminer, nous ajoutons la définition de MINLEVEL et ETAE 
nous ajoutons également la fonction qui va raffiner le maillage grâce
à la routine adapt_wavelet et nous appelons cette fonction à l'aide
d'un event.

~~~literatec
#define MINLEVEL 5
#define ETAE     1e-2 // error on free surface elevation (1 cm)

int adapt() {
  scalar et[];
  foreach()
    et[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({et});

  astats s = adapt_wavelet ({et}, (double[]){ETAE},
			    MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
}

event do_adapt (i++) adapt();
~~~
  

## Execution

Le programme est maintenant terminé. Avant de l'éxecuter, vous devez fixer un valeur raisonnable à MAXLEVEL (la valeur du raffinement maximal), cette valeur . Je vous conseille de commencer par une valeur petite puis d'augmenter jusqu'à avoir des temps de calcul de l'ordre de 5 min. Quelle valeur avez-vous choisi ? Cela correspond à quelle taille de cellule minimale ?

Tracez les courbes de l'évolution de la hauteur d'eau en fonction du temps pour toutes les stations de mesures. Que pouvez-vous dire pour Nice ? Pour toutes les autres stations ? Pour chaque station, donnez le temps pour lequel l'amplitude de la vague est maximale, que remarque-t-on ? Que se passe-t-il pour ce temps pour Calvi et Nice? Pourquoi ?
Regardez attentivement le film de votre simulation. Que pouvez-vous dire sur l'impact d'un tel tsunami sur la côte d'azur ? Donnez une explication qualitative de ce phénomène.

[Retour Sommaire Cours](http://basilisk.fr/sandbox/geoffroy/teaching/README)
*/
