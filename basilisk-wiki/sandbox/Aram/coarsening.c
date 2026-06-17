/** 
Simulation du coarsening des roll waves granulaires 
   avec les équations :

   $$ \frac{\partial h}{\partial t} + \frac{\partial (h \bar{u})}{\partial x} = 0 $$

   $$ \frac{\partial (h \bar{u})}{\partial t}
     + \frac{\partial (h \bar{u}^2)}{\partial x}
     = h g \sin(\zeta)
       - \mu h g \cos(\zeta)
       - \frac{\partial}{\partial x}\left(\frac{1}{2} h^2 g \cos(\zeta)\right)
       + \frac{\partial}{\partial x}\left(\nu h^{3/2} \frac{\partial \bar{u}}{\partial x}\right) $$

   Avec un bruit initial.
*/

#include "grid/cartesian1D.h"
#include "saint-venant.h"

#define G_GRAV 9.81              // Accélération de la pesanteur (m/s^2)
#define ZETA (35.1 * M_PI / 180.0) // Angle d'inclinaison du canal en radians (35.1°)
#define H0 0.0042                // Épaisseur initiale uniforme de la couche (4.2 mm)

// Paramètres de la friction de Pouliquen & Forterre
#define MU1 tan(32.9 * M_PI / 180.0) // Friction minimale \mu_1 = tan(32.9°)
#define MU2 tan(42.0 * M_PI / 180.0) // Friction maximale \mu_2 = tan(42.0°)
#define BETA 0.65                // Constante empirique \beta
#define L_CHAR 0.001             // Échelle de longueur caractéristique \mathcal{L} (1 mm)

// Paramètre visqueux granulaire
#define NU_VISC 2.4e-3           // Coefficient \nu (m^3/2 / s)

int main() {
  L0 = 50.0;
  X0 = 0.0;
  N = 2048;
  
  periodic(right);
  G = G_GRAV * cos(ZETA);
  
  run();
}

event init (t = 0) {
  double u0 = (BETA * pow(H0, 1.5) * sqrt(G) / L_CHAR) * ((tan(ZETA) - MU1) / (MU2 - tan(ZETA)));
  
  foreach() {
    double bruit = 0.02 * ((double)rand() / RAND_MAX - 0.5);
    h[] = H0 * (1.0 + bruit);
    u.x[] = (H0 * u0) / h[]; 
  }
  boundary({h, u.x});
}

event source (i++) {
  face vector v_flux[];
  
  foreach_face() {
    double h_f = 0.5 * (h[] + h[-1]);
    double du_dx = (u.x[] - u.x[-1]) / Delta;
    v_flux.x[] = NU_VISC * pow(h_f, 1.5) * du_dx;
  }
  
  foreach() {
    if (h[] > 1e-5) {
      double u_mag = fabs(u.x[]) + 1e-10;
      double mu = MU1 + (MU2 - MU1) / (1.0 + BETA * pow(h[], 1.5) * sqrt(G) / (L_CHAR * u_mag));
      double source_grav_fric = G * (tan(ZETA) - mu);
      double source_visqueuse = (v_flux.x[1] - v_flux.x[]) / (h[] * Delta);
      
      u.x[] += dt * (source_grav_fric + source_visqueuse);
    } else {
      u.x[] = 0.;
    }
  }
  boundary({u.x});
}

event outputs (t += 5; t <= 5000) {
  
  // 1. Écriture du profil complet (hauteur et vitesse)
  char name[80];
  sprintf(name, "profil_t_%g.dat", t);
  FILE * fp = fopen(name, "w");
  foreach() {
    fprintf(fp, "%g %g %g\n", x, h[] * 1000.0, u.x[]);
  }
  fclose(fp);

  // 2. Comptage des maxima locaux
  int nb_vagues = 0;
  boundary({h}); 
  
  foreach(reduction(+:nb_vagues)) {
    if (h[] > h[-1] && h[] > h[1] && h[] > H0 * 1.02) {
      nb_vagues++;
    }
  }

  // 3. Sauvegarde du nombre de vagues
  FILE * fp_waves = fopen("evolution_vagues.dat", (t == 20) ? "w" : "a");
  if (fp_waves) {
    fprintf(fp_waves, "%g %d\n", t, nb_vagues);
    fclose(fp_waves);
  }

  printf("Fichier %s sauvegardé | t = %g s | Vagues détectées = %d\n", name, t, nb_vagues);
}

// Événement Gnuplot 1 : profil h(x) à chaque pas de sauvegarde
event plot_profil (t += 20; t <= 5000) {
  char name_png[80], name_dat[80];
  sprintf(name_png, "profil_rollwaves_t%g.png", t);
  sprintf(name_dat, "profil_t_%g.dat", t);

  FILE * gp = popen("gnuplot", "w");
  if (gp) {
    fprintf(gp,
      "set terminal pngcairo size 800,400 enhanced font 'Verdana,10'\n"
      "set output '%s'\n"
      "set title 'Profil des Ondes de Roulis Granulaires - t = %g s'\n"
      "set xlabel 'Position x (m)'\n"
      "set ylabel 'Epaisseur h (mm)'\n"
      "set xrange [0:50]\n"
      "set yrange [2:12]\n"
      "set grid\n"
      "plot '%s' using 1:2 with lines lw 2 lc rgb 'blue' title 'h(x) a t = %g s'\n",
      name_png, t, name_dat, t);
    pclose(gp);
  }
}

// Événement Gnuplot 2 : courbe de coarsening N(t), tracée en fin de simulation
event plot_coarsening (t = end) {
  FILE * gp = popen("gnuplot", "w");
  if (gp) {
    fprintf(gp,
      "set terminal pngcairo size 700,400 enhanced font 'Verdana,10'\n"
      "set output 'courbe_coarsening.png'\n"
      "set title 'Evolution du nombre de vagues au cours du temps (Coarsening)'\n"
      "set xlabel 'Temps t (s)'\n"
      "set ylabel 'Nombre de vagues N'\n"
      "set xrange [0:5000]\n"
      "set yrange [0:100]\n"
      "set grid\n"
      "plot 'evolution_vagues.dat' using 1:2 with linespoints lw 2 lc rgb 'red'"
      " title 'N(t)'\n");
    pclose(gp);
  }
}

// Événement Gnuplot 3 : animation GIF
event plot_animation (t = end) {
  FILE * gp = popen("gnuplot", "w");
  if (gp) {
    fprintf(gp,
      "set terminal gif animate delay 20 loop 0 size 800,400 enhanced font 'Verdana,10'\n"
      "set output 'rollwaves_evolution.gif'\n"
      "set xlabel 'Position x (m)'\n"
      "set ylabel 'Epaisseur h (mm)'\n"
      "set xrange [0:50]\n"
      "set yrange [2:12]\n"
      "set grid\n");

    for (double t_anim = 5; t_anim <= 5000; t_anim += 5) {
      fprintf(gp,
        "set title 'Evolution des Ondes de Roulis Granulaires - t = %g s'\n"
        "plot 'profil_t_%g.dat' using 1:2 with lines lw 2 lc rgb 'blue'"
        " title 'h(x) a t = %g s'\n",
        t_anim, t_anim, t_anim);
    }
    pclose(gp);
  }
}


/**
~~~gnuplot Coarsening

set terminal pngcairo size 700,400 enhanced font 'Verdana,10'
set output 'courbe_coarsening.png'

set title "Évolution du nombre de vagues au cours du temps (Coarsening)"
set xlabel "Temps t (s)"
set ylabel "Nombre de vagues N"

set xrange [0:2500]
set yrange [0:100]
set grid

plot "evolution_vagues.dat" using 1:2 with linespoints lw 2 lc rgb "red" title "N(t)"

~~~
*/




/**
~~~gnuplot Animation Coarsening
# Configuration du terminal pour générer un GIF animé
# delay 40 = 0.4 seconde entre chaque image | loop 0 = répétition infinie
set terminal gif animate delay 30 loop 0 size 800,400 enhanced font 'Verdana,10'
set output 'rollwaves_evolution.gif'

set xlabel "Position x (m)"
set ylabel "Épaisseur h (mm)"
set xrange [0:50]
set yrange [2:12]
set grid

# Boucle pour lire les fichiers de t = 100s à t = 1000s
do for [t=20:2500:20] {
    
    set title sprintf("Évolution des Ondes de Roulis Granulaires - t = %d s", t)
    fichier_data = sprintf("profil_t_%d.dat", t)
    
    # Trace la courbe pour l'instant t actuel
    plot fichier_data using 1:2 with lines lw 2 lc rgb "blue" title sprintf("h(x) à t = %d s", t)
}


~~~
*/

/**
## Bibliography
  
 * D. Razis, A. N. Edwards, J. M. N. T. Gray, and Ko van der Weele
 "Arrested coarsening of granular roll waves"
*/