/**
# Solveurs de Riemann et cas test

Dans ce tp, vous allez coder une bibliothèque contenant différents
solveurs de Riemann que vous testerez sur le cas test de la rupture
d'une barrage.

## Solveurs de Riemann.

Dans un programme que vous appellerez riemann_votrenom.h, codez les
solveurs de Riemann vu pendant le cours précédant la séance de TP. Ce
fichier aura la forme d'une bibliothèque de fonctions et aura pour but 
d'être appelée par d'autres programmes. Les fonctions devront avoir la même forme que
celles déjà présentes dans la bibliothèque
[riemann.h](http://basilisk.fr/src/riemann.h) de
Basilisk. Inspirez-vous en.

Appelez votre professeur lorsque vous avez fini cette étape.

## Cas test : Rupture d'un barrage 1D

### Les includes

Dans un programme nommé barrage_votrenom.c, vous allez coder la
rupture d'un barrage 1D. Commencez par faire un include de la bonne grille :
cartésienne, non-adaptative et à 1 dimension, puis du fichier
"predictor-corrector.h". 

Exercice : A quoi sert le fichier predictor-corrector.h ? Où se trouve t-il ?

### Surcharge de la fonction advance et update

Nous allons surcharger les fonctions advance et update utilisées dans
"predictor-corrector.h". Copiezet collez tout le code suivant :*/

// Nous commençons par déclarer les champs qui vont évoluer
scalar h[];
vector u[];
scalar * evolving = {h, u};

// Ainsi que quelques variables nécessaires
double G = 1.;
double dry = 1e-6;

// La fonction advance
static void advance_saint_venant (scalar * output, scalar * input, 
				  scalar * updates, double dt)
{
  scalar hi = input[0], ho = output[0], dh = updates[0];
  vector ui = vector(input[1]), uo = vector(output[1]), dhu = vector(updates[1]);

  foreach() {
    double hold = hi[];
    ho[] = hold + dt*dh[];
    if (ho[] > dry)
      foreach_dimension()
	uo.x[] = (hold*ui.x[] + dt*dhu.x[])/ho[];
    else
      foreach_dimension()
	uo.x[] = 0.;
  }
  boundary ({ho, uo});
}
/**
Ici, vous allez ajouter un include de votre bibliothèque de solveurs de Riemann.
*/
// La fonction update
double update_saint_venant (scalar * evolving, scalar * updates, double dtmax)
{
 
 
  scalar h = evolving[0];
  vector u = vector(evolving[1]);

  
  face vector Fh[];
  tensor Fq[];
  
  
  foreach_face () {
    if (h[] > dry || h[-1,0] > dry) {
      
      /** Appel du solveur de Riemann ici : */
      double fh, fu;
      rusanov (h[-1,0], h[], u.x[-1,0], u.x[], Delta, &fh, &fu, &dtmax);
      Fh.x[]   = fh;
      Fq.x.x[] = fu;
    }
    else // dry
      Fh.x[] = Fq.x.x[] = 0.;
  }
  boundary_flux ({Fh, Fq});

  scalar dh = updates[0];
  vector dhu = vector(updates[1]);

  foreach() {
    dh[] = (Fh.x[]  - Fh.x[1,0])/Delta;
    foreach_dimension()
      dhu.x[] = (Fq.x.x[]  - Fq.x.x[1,0] )/Delta;
  }
  return dtmax;
}

// Nous surchargeons les fonctions ici
event defaults (i = 0)
{
  advance = advance_saint_venant;
  update = update_saint_venant;
}
/**

### La fonction main()

Dans le fonction main(), fixez la taille du domaine à 60 m et
l'origine tel que le domaine soit centré sur X = 0. Fixez la gravité de manière à avoir le même
système d'unité que le SI. Fixez le nombre d'éléments N de façon à avoir une éxecution de l'ordre de la dizaine de secondes. Finissez la fonction main() par l'appelle de la fonction run().

### Les Conditions Initiales

Au temps t = 0, toute la partie des x négatifs a une hauteur de $h_0 = 1 m$, alors que la
partie des X positif est sèche. La vitesse est nulle partout. La
topographie est plane.

### Les Conditions aux Bords

Fixez les bonnes conditions aux bords (à l'extérieur des events).

### Solution exacte du problème

Ce problème possède une solution exacte, qui est, pour la hauteur d'eau :

si $x < -t \sqrt(G*h0)$ alors :  $h = h0$ 

si $-t \sqrt(G*h0) \leq x < 2t\sqrt(G*h0)$ alors

$$
h = \frac{x^2}{9 G t^2} - \frac{4x}{9t}\sqrt{\frac{h0}{G}} + \frac{4h0}{9}
$$

et si $2 t \sqrt(G*h0) \leq x$ alors :  $h = 0$

Déclarez deux nouveaux champs scalaires que vous nommerez "hana" et
"error". Créez un event qui s'executera à chaque pas de temps dans
lequel vous calculerez la solution analytique du problème. Vous
stoquerez cette solution dans le champ hana[]. Dans le champ error[] vous
stoquerez la valeur absolute de la différence entre la solution analytique 
et la solution numérique

                                       
### Outputs 

Nous utilisons gnuplot pour faire une animation de cette simulation.
*/
event plot (t< 10;t+=0.01) {
    printf("set title 'Rupture de barrage 1D ----- t= %.1lf '\n"
	 "p[%g:%g][-0.1:1.5]  '-' u 1:2 t'h numerical (m)' w l lt 3,"
	 "'' u 1:3 t'h analytical' w l lt 4,"
	 "'' u 1:($4*10) t'10 * error' w l lt 1\n",t,X0+1,X0+L0);
    foreach()
      printf (" %g %g %g %g %g \n", x, h[], hana[], error[], t);
    printf ("e\n\n");
}
/**
### Calcul de l'erreur

Définissez une (ou plusieurs) norme qui permet de comparer les erreurs suivant le solveur utilisé
à chaque pas de temps. Calculez-là et imprimez-là dans le deuxième
fichier de sortie. N'oubliez pas d'imprimer également le temps.

## Résultats

Faites autant de run de votre programme que vous avez codé de solveurs de Riemann en en utilisant un différent à chaque fois.
Sur un même graphique, analysez les erreurs pour les différents solveurs
de Riemann. Lequel est le meilleur pour ce cas ?

## Remerciements et crédits
Ce TP est très largement inspiré de [celui de Pierre Yves Lagrée](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/code_C_saintvenant.pdf) sur le même sujet, avec son aimable autorisation.


## [Retour Sommaire Cours](http://basilisk.fr/sandbox/geoffroy/teaching/README)
*/

