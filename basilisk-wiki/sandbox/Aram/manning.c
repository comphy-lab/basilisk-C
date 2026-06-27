/** 
## Simulation des trains d'ondes (roll-waves) en eau claire

Pour décrire les instabilités d'un écoulement à surface libre sur plan incliné on utilise les équations de Saint-Venant:
   $$
   \partial_t h + \partial_x (hu) = 0
   $$
   $$
   \partial_t (hu) + \partial_x \left( hu^2 + \frac{1}{2} g h^2 \right) = -gh      \, \partial_x z - \frac{g n^2}{h^{1/3}} u|u|
   $$  
   avec pour le frottement:
$$
S_0 = g n^2 \frac{Fr^2}{h_0^{1/3}}
$$
   et un coefficient de frottement de Manning pour modéliser le frottement turbulent:
$$
n = 0.02
$$
Pour un nombre de Froude supérieur à 2, on montre que l'écoulement est instable et qu'un train d'ondes apparait, appelé "roll-waves" en anglais.


## Code
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"

#ifndef FR
#  define FR 2.2
#endif

#define EPS    0.005          // amplitude perturbation ε = 0.5%
#define K_MANN (5.*pi)        // nombre d'onde k = 5π
#define Q0     0.001          // débit q (m²/s)
#define N_MANN 0.02           // coefficient de Manning n (m^{-1/3}.s)

double h0, u0, S0;

FILE *fp_profils = NULL;
FILE *fp_diag    = NULL;

int main()
{
  X0 = 0.;
  L0 = 2.;
  G  = 9.81;
  N  = 2500;    // Δx = 8e-4 m (identique à Delestre)

  h0 = pow(Q0 / (sqrt(G) * FR), 2./3.);
  u0 = Q0 / h0;
  S0 = G * N_MANN * N_MANN * FR * FR / pow(h0, 1./3.);

  DT = 0.5 * (L0/N) / (u0 + sqrt(G * h0));
  
  char fname_profils[64], fname_diag[64];
  snprintf(fname_profils, sizeof(fname_profils), "out_mann_Fr%.4g.dat", (double)FR);
  snprintf(fname_diag,    sizeof(fname_diag),    "diag_mann%.4g.dat",   (double)FR);

  fp_profils = fopen(fname_profils, "w");
  fp_diag    = fopen(fname_diag,    "w");

  if (!fp_profils || !fp_diag) {
    fprintf(stderr, "Erreur : impossible d'ouvrir les fichiers de sortie\n");
    return 1;
  }

  fprintf(fp_profils, "# x  h  u  zb  t\n");
  fprintf(fp_diag,    "# t  log(amplitude)\n");

  fprintf(stderr, "=== PARAMETRES MANNING ===\n");
  fprintf(stderr, "Fr    = %.4f\n", (double)FR);
  fprintf(stderr, "h0    = %.8f m\n", h0);
  fprintf(stderr, "u0    = %.6f m/s\n", u0);
  fprintf(stderr, "S0    = %.6f\n", S0);
  fprintf(stderr, "n     = %.4f\n", N_MANN);
  fprintf(stderr, "DT    = %.6f s\n", DT);
  fprintf(stderr, "Check Fr = u0/sqrt(g*h0) = %.6f\n", u0/sqrt(G*h0));
  fprintf(stderr, "Critere Craya : Fr %s 1.5 => %s\n",
          (double)FR > 1.5 ? ">" : ((double)FR < 1.5 ? "<" : "="),
          (double)FR > 1.5 ? "INSTABLE" : ((double)FR < 1.5 ? "STABLE" : "CRITIQUE"));
  fprintf(stderr, "Sortie profils : %s\n", fname_profils);
  fprintf(stderr, "Sortie diag    : %s\n", fname_diag);
  fprintf(stderr, "==================\n\n");

  periodic(right);
  run();

  fclose(fp_profils);
  fclose(fp_diag);

  fprintf(stderr, "Fichiers ecrits : %s  %s\n", fname_profils, fname_diag);
  return 0;
}

event init (i = 0)
{
  foreach() {
    zb[]  = 0.;
    h[]   = h0 * (1. + EPS * sin(K_MANN * x));
    u.x[] = Q0 / h[];   // débit constant : u = q/h
  }
}

event acceleration (i++) {
  foreach()
    u.x[] += dt * G * S0;
}

event friction (i++) {
  foreach() {
    if (h[] > dry) {
      double Cf = G * N_MANN * N_MANN * fabs(u.x[]) / pow(h[], 4./3.);
      foreach_dimension()
        u.x[] /= (1. + dt * Cf);
    }
  }
}

event output (t = 0; t <= 50; t += 10) {
  fprintf(fp_profils, "# t = %g\n", t);
  foreach()
    fprintf(fp_profils, "%.12g %.12g %.12g %.12g %.12g\n",
            x, h[], u.x[], zb[], t);
  fprintf(fp_profils, "\n\n");
  fflush(fp_profils);
}

event diag (t = 0; t <= 50; t += 0.5) {
  double amp = 0.;
  foreach()
    if (fabs(h[] - h0) > amp) amp = fabs(h[] - h0);
  if (amp > 1e-15) {
    fprintf(fp_diag, "%.4f %.8f\n", t, log(amp));
    fflush(fp_diag);
  }
}


/**
## Results

Développement d'un train d'onde, l'échelle de la hauteur est amplifiée par rapport à la direction de l'écoulement 

~~~pythonplot Fr = 2.2
import numpy as np
import matplotlib.pyplot as plt
import os


fname = f"out_mann_Fr2.2.dat"

T_LIST = [0, 10, 20, 30, 40, 50]


if not os.path.isfile(fname):
    print(f"  Fichier manquant : {fname}")
else:
    blocs = []
    courant = []
    
    with open(fname) as f:
        for ligne in f:
            s = ligne.strip()
            if s.startswith('#') or s == '':
                if courant:
                    blocs.append(np.array(courant))
                    courant = []
            else:
                vals = list(map(float, s.split()))
                if len(vals) >= 2:
                    courant.append(vals)
    if courant:
        blocs.append(np.array(courant))
    if not blocs:
        print("  Aucun bloc de données trouvé dans le fichier.")
    else:
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.set_title(f"Roll-waves frottement Manning — Fr=2.2 — profils h(x)", fontsize=13)

        # On trace au maximum les 6 premiers blocs temporels extraits
        for idx, data in enumerate(blocs[:6]):
            x = data[:, 0]
            h = data[:, 1]
            t = T_LIST[idx] if idx < len(T_LIST) else idx * 10
            ax.plot(x, h, lw=1.8, label=f't={t}s')

        ax.set_xlabel('x (m)', fontsize=12)
        ax.set_ylabel('h (m)', fontsize=12)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        
        plt.savefig("profils_superposes_Fr2.2.png")
       
~~~    
*/

/**
## Bibliography
  
 * [J. Needham and J. H. Merkin](http://www.jstor.org/stable/2397934)
 On roll waves down an open inclined channel
 Proc. R. Soc. Lond. A 394, 259-278 (1984)
 
 * [Olivier Delestre](https://theses.hal.science/tel-00531377v3)
 "Simulation du ruissellement d’eau de pluie sur des surfaces agricoles"
 
 * [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 Sorbonne-Université
 
*/