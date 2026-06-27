/**
# Écoulement granulaire sec sur une pente à \(35^\circ\)

Ce programme simule en une dimension un tas de billes de verre initialement
retenu par une porte, puis libéré sur une pente inclinée suivie d'une partie
horizontale. La géométrie et les dimensions initiales correspondent au cas
sec à \(35^\circ\) de Poulain *et al.* (2023), pour des billes de verre de
diamètre \(d=4\) mm.

Le solveur est celui des équations de Saint-Venant de Basilisk, auquel on
ajoute une loi de frottement basal de Pouliquen--Forterre. Les paramètres de
cette loi sont ceux utilisés dans l'article pour les billes de 4 mm.

> **Important.** Malgré le nom historique du fichier, ce programme n'est pas
> le code HySEA. HySEA est un code bidimensionnel à deux couches
> eau--granulaire. Ici, on résout seulement le cas sec, en 1D, avec la
> formulation cartésienne Saint-Venant de Basilisk. Il constitue donc un
> modèle réduit adapté à la comparaison avec la partie granulaire sèche de
> l'article.

## Équations résolues

Les inconnues sont l'épaisseur verticale \(h(x,t)\) et la vitesse moyenne
horizontale \(u(x,t)\). Le système de Saint-Venant avec topographie
\(z_b(x)\) s'écrit

$$
\frac{\partial h}{\partial t}
+ \frac{\partial (hu)}{\partial x} = 0,
$$

$$
\frac{\partial (hu)}{\partial t}
+ \frac{\partial}{\partial x}
\left(hu^2 + \frac{G h^2}{2}\right)
= -Gh\frac{\partial z_b}{\partial x}
  -Gh\,\mu(h,Fr)\,\operatorname{sgn}(u).
$$

Basilisk traite le premier terme source, lié à la topographie, avec le
solveur `saint-venant.h`. Le frottement est appliqué après chaque pas de
temps par un découplage d'opérateurs, sous la forme

$$
\frac{\partial u}{\partial t}
= -G\,\mu(h,Fr)\,\operatorname{sgn}(u).
$$

Le nombre de Froude granulaire utilisé dans la loi de frottement est

$$
Fr = \frac{|u|}{\sqrt{G h\cos\theta_{\mathrm{loc}}}},
\qquad
\cos\theta_{\mathrm{loc}}
= \frac{1}{\sqrt{1+(\partial_x z_b)^2}}.
$$

## Loi de frottement de Pouliquen--Forterre

On définit \(\mu_i=\tan\delta_i\), \(L=1.3d\), puis

$$
\mu_{\mathrm{stop}}(h)
= \mu_1 + \frac{\mu_2-\mu_1}{1+h/L},
\qquad
\mu_{\mathrm{start}}(h)
= \mu_3 + \frac{\mu_2-\mu_1}{1+h/L}.
$$

La loi comporte trois régimes :

$$
\mu(h,Fr)=
\begin{cases}
\displaystyle
\mu_1 + \frac{\mu_2-\mu_1}
{1+\beta h/(LFr)}, & Fr\geq\beta,\\[0.8em]
\displaystyle
\mu_{\mathrm{start}}
+ \left(\mu_{\mathrm{stop}}-\mu_{\mathrm{start}}\right)
\left(\frac{Fr}{\beta}\right)^\xi,
& 0<Fr<\beta,\\[0.8em]
\displaystyle
\min\left(\mu_{\mathrm{start}},
\left|\partial_x(z_b+h)\right|\right), & Fr=0.
\end{cases}
$$

Très près de \(Fr=0\), une petite régularisation numérique est ajoutée afin
que le raccord entre l'état statique et le régime intermédiaire soit continu.
Le paramètre `epsilon_pf` ne fait donc **pas** partie de la loi expérimentale
originale : il sert uniquement à contrôler la robustesse numérique.

## Adimensionnement

Les variables de Basilisk sont adimensionnées selon

$$
x=L_{\mathrm{ref}}x^*,\qquad
h=L_{\mathrm{ref}}h^*,\qquad
t=T_{\mathrm{ref}}t^*,\qquad
u=\frac{L_{\mathrm{ref}}}{T_{\mathrm{ref}}}u^*,
$$

avec

$$
T_{\mathrm{ref}}=\sqrt{\frac{L_{\mathrm{ref}}}{g}}.
$$

Le code prend \(L_{\mathrm{ref}}=1\) m. Ainsi, les longueurs numériques
restent égales aux valeurs exprimées en mètres, mais la gravité
adimensionnée vaut exactement \(G=gT_{\mathrm{ref}}^2/L_{\mathrm{ref}}=1\).
Les fichiers de sortie reconvertissent explicitement les données en unités
physiques : mètres et secondes.

## Géométrie et condition initiale

Le fond est horizontal en aval, puis devient une pente de \(35^\circ\) en
remontant vers la porte. La cassure pente--horizontal est remplacée par un
raccord quadratique de largeur `transition_width_phys`, continu en altitude
et en pente. Cette régularisation évite de créer une singularité numérique au
changement brutal de topographie ; une stratégie analogue est utilisée dans
Poulain *et al.* (2023).

Le tas est placé entre \(x=x_{\rm gate}\) et \(x=x_{\rm gate}+L_i\), au
repos initialement. Sa surface libre est horizontale.

## Code
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "grid/cartesian1D.h"
#include "saint-venant.h"

/**
## Paramètres physiques, numériques et géométriques

Toutes les quantités portant le suffixe `_phys` sont exprimées en unités SI.
Les variables sans suffixe sont les variables adimensionnées effectivement
manipulées par Basilisk.
*/

/* Échelles de référence et accélération de la pesanteur. */
double Lref = 1.;                 /* m */
double Tref;                      /* s */
const double gravity_phys = 9.81; /* m s^-2 */

/* Géométrie expérimentale : cas sec à 35 degrés, billes de 4 mm. */
const double theta_phys = 35.*pi/180.;
const double xgate_phys = 0.;
const double zrelease_phys = 0.15;       /* hauteur du fond à la porte [m] */
const double pile_length_phys = 0.131;   /* Li [m] */
const double pile_height_phys = 0.103;   /* Hi [m] */
const double transition_width_phys = 0.10; /* largeur du raccord [m] */

/* Paramètres de la loi de Pouliquen--Forterre. */
const double delta1 = 22.1*pi/180.;
const double delta2 = 31.8*pi/180.;
const double delta3 = 23.3*pi/180.;
const double beta_pf = 0.136;
const double xi_pf = 1e-3;
const double grain_diameter_phys = 4e-3; /* d [m] */

/*
   Pour tester la calibration empirique de Poulain et al. dans le modèle
   cartésien, il suffit de remplacer les trois angles ci-dessus par :

   delta1 = 27.1*pi/180.;
   delta2 = 36.8*pi/180.;
   delta3 = 28.3*pi/180.;
*/

/* Paramètres de sortie, tous exprimés initialement en unités physiques. */
const double tmax_phys = 3.;
const double output_dt_phys = 0.02;
const double h_front_threshold_phys = 1e-4;
const char * animation_gif = "animate_dry35_muv_hysea.gif";
const char * animation_mp4 = "animate_dry35_muv_hysea.mp4";

/* Variables adimensionnées construites dans main(). */
double theta, slope_tan;
double xgate, zrelease, pile_length, pile_height;
double xbreak, transition_width, transition_start, transition_end;
double Lgrain, mu1, mu2, mu3;
double epsilon_pf = 1e-2;
double h_front_threshold;
double tmax;

FILE * animation_fp = NULL;

/**
## Topographie lissée

Si \(x_b\) est la position théorique de la cassure et \(\ell_t\) la largeur
du raccord, la topographie est

$$
z_b(x)=
\begin{cases}
0, & x\leq x_b-\ell_t/2,\\
\tan\theta\,\ell_t\,s^2/2,
& x_b-\ell_t/2 < x < x_b+\ell_t/2,\\
\tan\theta\,(x-x_b), & x\geq x_b+\ell_t/2,
\end{cases}
\qquad
s=\frac{x-(x_b-\ell_t/2)}{\ell_t}.
$$

La pente est donc nulle à l'entrée du raccord et égale à
\(\tan\theta\) à sa sortie.
*/

static double smooth_bottom (double xp)
{
  if (xp <= transition_start)
    return 0.;

  if (xp >= transition_end)
    return slope_tan*(xp - xbreak);

  double s = (xp - transition_start)/transition_width;
  return slope_tan*transition_width*sq(s)/2.;
}

static double smooth_bottom_slope (double xp)
{
  if (xp <= transition_start)
    return 0.;

  if (xp >= transition_end)
    return slope_tan;

  double s = (xp - transition_start)/transition_width;
  return slope_tan*s;
}

/**
## Loi de Pouliquen--Forterre régularisée

`smoothstep()` vaut zéro en \(s=0\), un en \(s=1\), et possède une dérivée
nulle aux deux extrémités. Il permet de raccorder en douceur l'état statique
à la branche intermédiaire dans l'intervalle
\(0<Fr<\epsilon_{\rm PF}\).
*/

static double smoothstep (double s)
{
  s = max(0., min(1., s));
  return sq(s)*(3. - 2.*s);
}

static double mu_pouliquen_forterre (double Fr, double hh,
                                     double dzb_dx, double dh_dx)
{
  double h_sur_L = hh/(Lgrain + 1e-30);

  double mu_stop_h =
    mu1 + (mu2 - mu1)/(1. + h_sur_L);

  double mu_start_h =
    mu3 + (mu2 - mu1)/(1. + h_sur_L);

  /* Critère statique fondé sur la pente de la surface z_b + h. */
  double surface_slope = fabs(dzb_dx + dh_dx);
  double mu_static = min(mu_start_h, surface_slope);

  if (Fr <= 0.)
    return mu_static;

  /* Régime intermédiaire : 0 < Fr < beta. */
  double mu_intermediate =
    mu_start_h + (mu_stop_h - mu_start_h)*pow(Fr/beta_pf, xi_pf);

  /* Raccord numérique : seulement actif pour Fr < epsilon_pf. */
  if (Fr < epsilon_pf)
    return mu_static + smoothstep(Fr/epsilon_pf)*
      (mu_intermediate - mu_static);

  if (Fr < beta_pf)
    return mu_intermediate;

  /* Régime dynamique : Fr >= beta. */
  return mu1 + (mu2 - mu1)/(1. + beta_pf*hh/(Lgrain*Fr));
}

/**
## Écriture d'un profil

Chaque ligne du fichier envoyé sur la sortie d'erreur contient :

```
x [m]    h [m]    t [s]    z_b [m]    u [m/s]
```

La sortie d'erreur est volontairement réservée aux profils afin de pouvoir
séparer simplement les deux jeux de données lors de l'exécution du code :
la sortie standard contient uniquement l'évolution du front.
*/

static void write_profile (double time_phys)
{
  foreach()
    fprintf (stderr, "%.12g %.12g %.12g %.12g %.12g\n",
             x*Lref,
             h[]*Lref,
             time_phys,
             zb[]*Lref,
             u.x[]*Lref/Tref);
}

/**
## Programme principal

Le premier argument de la ligne de commande fixe la résolution \(N\).
Le second, facultatif, fixe `epsilon_pf`.
*/

int main (int argc, char * argv[])
{
  /* Passage des données physiques aux variables adimensionnées. */
  Tref = sqrt(Lref/gravity_phys);
  G = gravity_phys*sq(Tref)/Lref; /* donc G = 1 avec le choix de Tref */

  theta = theta_phys;
  slope_tan = tan(theta);

  xgate = xgate_phys/Lref;
  zrelease = zrelease_phys/Lref;
  pile_length = pile_length_phys/Lref;
  pile_height = pile_height_phys/Lref;

  transition_width = transition_width_phys/Lref;
  xbreak = xgate - zrelease/slope_tan;
  transition_start = xbreak - transition_width/2.;
  transition_end = xbreak + transition_width/2.;

  Lgrain = 1.3*grain_diameter_phys/Lref;
  mu1 = tan(delta1);
  mu2 = tan(delta2);
  mu3 = tan(delta3);

  h_front_threshold = h_front_threshold_phys/Lref;
  tmax = tmax_phys/Tref;

  /* Domaine physique [-0.65, 0.20] m. */
  X0 = -0.65/Lref;
  L0 = 0.85/Lref;
  N = argc > 1 ? atoi(argv[1]) : 1024;

  if (argc > 2)
    epsilon_pf = atof(argv[2]);

  if (N <= 0) {
    fprintf (stderr, "Erreur : N doit etre strictement positif.\n");
    return 1;
  }

  if (epsilon_pf <= 0. || epsilon_pf >= beta_pf) {
    fprintf (stderr,
            "Erreur : epsilon_pf doit verifier 0 < epsilon_pf < beta_pf = %g.\n",
            beta_pf);
    return 1;
  }

  /* Pas de temps maximal adimensionné : 10^-4 Tref. */
  DT = 1e-4;

  run();
}

/**
## Condition initiale

Le tas est au repos derrière la porte. Sa surface libre initiale est située à
l'altitude \(z_{\rm release}+H_i\). Comme le modèle est cartésien,
`h[]` représente une hauteur verticale et non l'épaisseur normale à la pente.
*/

event init (i = 0)
{
  foreach() {
    zb[] = smooth_bottom(x);

    double in_pile = (x >= xgate && x <= xgate + pile_length);
    double surface0 = zrelease + pile_height;

    h[] = in_pile ? max(surface0 - zb[], 0.) : 0.;
    u.x[] = 0.;
  }

  boundary ({zb, h, u});
}

/**
## Frottement basal

À chaque itération, on calcule le nombre de Froude local, puis le coefficient
\(\mu(h,Fr)\). La vitesse est ensuite diminuée de
\(G\mu\,\mathrm{d}t\), sans jamais changer son signe. Cette dernière
précaution évite une oscillation numérique autour de \(u=0\).
*/

event friction (i++)
{
  foreach() {
    if (h[] <= dry) {
      u.x[] = 0.;
      continue;
    }

    /* La pente du fond est analytique ; celle de h est discrétisée. */
    double dzb_dx = smooth_bottom_slope(x);
    double dh_dx = (h[1] - h[-1])/(2.*Delta);

    double cos_local = 1./sqrt(1. + sq(dzb_dx));
    double speed = fabs(u.x[]);
    double Fr = speed/(sqrt(G*h[]*cos_local) + 1e-30);

    double mu_eff =
      mu_pouliquen_forterre(Fr, h[], dzb_dx, dh_dx);

    double new_speed = max(speed - G*mu_eff*dt, 0.);
    u.x[] = speed > 1e-30 ? new_speed*u.x[]/speed : 0.;
  }

  boundary ({u});
}

/**
## Profils aux temps de comparaison

Les profils sont enregistrés à \(t=0\), \(0.2\), \(0.4\) et \(3\) s,
comme dans les instantanés de comparaison avec l'expérience. Chaque événement
est explicitement programmé en temps adimensionné afin que les temps inscrits
dans le fichier restent les temps physiques demandés.
*/

event profile_t0 (t = 0.)
{
  write_profile(t*Tref);
}

event profile_t02 (t = 0.2/Tref)
{
  write_profile(t*Tref);
}

event profile_t04 (t = 0.4/Tref)
{
  write_profile(t*Tref);
}

event profile_t30 (t = 3./Tref)
{
  write_profile(t*Tref);
}

/**
## Position du front

L'écoulement se propage vers les \(x\) décroissants. Le front est donc le
plus petit \(x\) tel que \(h>h_{\rm seuil}\). La sortie standard contient
simplement :

```
t [s]    run-out depuis la porte [m]
```
*/

event front (t = 0.; t <= tmax; t += output_dt_phys/Tref)
{
  double xf = xgate;

  foreach (reduction(min:xf))
    if (h[] > h_front_threshold)
      xf = min(xf, x);

  fprintf (stdout, "%.12g %.12g\n", t*Tref, (xgate - xf)*Lref);
}

/**
## Animation Gnuplot

Une image est ajoutée toutes les `output_dt_phys` secondes. Le GIF est ensuite
converti automatiquement en MP4 avec `ffmpeg`. Les messages éventuels de
Gnuplot et de FFmpeg sont redirigés vers des fichiers journaux, afin de ne pas
polluer les fichiers de données écrits sur `stdout` et `stderr`.
*/

event animatedplot (t = 0.; t <= tmax; t += output_dt_phys/Tref)
{
  if (!animation_fp)
    animation_fp = popen("gnuplot > gnuplot_animation.log 2>&1", "w");

  if (!animation_fp)
    return;

  if (t == 0.)
    fprintf (animation_fp,
      "set term gif animate delay 5 enhanced font 'Liberation Sans,14'\n"
      "set output '%s'\n"
      "set grid\n",
      animation_gif);

  fprintf (animation_fp,
    "set title 'Saint-Venant sec + Pouliquen-Forterre - theta = 35 deg - t = %.2f s' font 'Liberation Sans,18'\n"
    "set xlabel 'x (m)' font 'Liberation Sans,16'\n"
    "set ylabel 'z (m)' font 'Liberation Sans,16'\n"
    "plot [%g:%g][%g:%g] "
    "'-' u 1:2:3 t 'sol' w filledcurves lc rgb '#E2D1C4', "
    "'-' u 1:2:3 t 'tas granulaire' w filledcurves lc rgb '#6B3E26', "
    "'-' u 1:2 t 'topographie' w l lc rgb '#000000' lw 1.2\n",
    t*Tref,
    X0*Lref, (X0 + L0)*Lref,
    -0.05, zrelease_phys + pile_height_phys + 0.10);

  foreach()
    fprintf (animation_fp, "%g %g %g\n",
             x*Lref, -0.05, zb[]*Lref);
  fprintf (animation_fp, "e\n");

  foreach() {
    if (h[] > h_front_threshold)
      fprintf (animation_fp, "%g %g %g\n",
               x*Lref, zb[]*Lref, (zb[] + h[])*Lref);
    else
      fprintf (animation_fp, "\n");
  }
  fprintf (animation_fp, "e\n");

  foreach()
    fprintf (animation_fp, "%g %g\n", x*Lref, zb[]*Lref);
  fprintf (animation_fp, "e\n\n");

  fflush(animation_fp);
}

static void finalize_animation (void)
{
  if (animation_fp) {
    fprintf (animation_fp, "unset output\n");
    pclose(animation_fp);
    animation_fp = NULL;
  }

  char command[512];
  snprintf(command, sizeof(command),
           "ffmpeg -hide_banner -loglevel error -y -i '%s' "
           "-movflags +faststart -pix_fmt yuv420p '%s' "
           "2> ffmpeg_animation.log",
           animation_gif, animation_mp4);

  system(command);
}

/**
## Fin du calcul

La fermeture explicite du tube Gnuplot est indispensable : elle écrit la
dernière image du GIF avant la conversion en MP4.
*/

event stop (t = tmax)
{
  finalize_animation();
  return 1;
}

/**
## Compilation et exécution

~~~bash
qcc -O2 -Wall dry_landslide_muv_hysea_documented.c \\
    -o dry_landslide_muv_hysea_documented -lm

mkdir -p ../../Data/1D/article_poulain

./dry_landslide_muv_hysea_documented 2048 1e-2 \\
  > ../../Data/1D/article_poulain/front_dry35_muv_pf.dat \\
  2> ../../Data/1D/article_poulain/profiles_dry35_muv_pf.dat
~~~

Le premier argument est le nombre de mailles et le deuxième est
`epsilon_pf`. Pour tester l'influence de la régularisation, on peut relancer
avec, par exemple, `5e-3`, `1e-3` et `1e-4`.

Les fichiers suivants sont créés dans le répertoire d'exécution :

* `animate_dry35_muv_hysea.gif` : animation GIF ;
* `animate_dry35_muv_hysea.mp4` : même animation au format MP4 ;
* `gnuplot_animation.log` et `ffmpeg_animation.log` : journaux éventuels.

## Références

* P. Poulain *et al.*, *Performance and limits of a shallow-water model for
  landslide-generated tsunamis: from laboratory experiments to simulations
  of flank collapses at Montagne Pelée (Martinique)*, Geophysical Journal
  International, 233, 796--825, 2023.
* O. Pouliquen & Y. Forterre, *Friction law for dense granular flows:
  application to the motion of a mass down a rough inclined plane*, Journal
  of Fluid Mechanics, 453, 133--151, 2002.
*/
