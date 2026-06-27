/** 
## Simulation des trains d'ondes (roll-waves) en eau claire
Pour décrire les instabilités d'un écoulement à surface libre sur plan incliné on utilise les équations de Saint-Venant:
   $$
   \partial_t h + \partial_x (hu) = 0
   $$

   $$
   \partial_t (hu)
   + \partial_x \left( hu^2 + \frac{g h^2}{2} \right)
   = -gh\,\partial_x z - \frac{f}{8}u|u|
   $$
   
avec pour le frottement:
$$
S_0 = \frac{f}{8}\,\frac{u^2}{g h} = \frac{f}{8} Fr^2
$$
   et un coefficient de frottement de Darcy--Weisbach pour modéliser le frottement turbulent:
$$
f = 0.048
$$
Pour un nombre de Froude supérieur à 1.5, on montre que l'écoulement est instable et qu'un train d'ondes apparait, appelé "roll-waves" en anglais.



## Code
*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"

#ifndef FR
#  define FR 2.5    // on peut changer Fr pour voir, le cas d'un écoulement instable (Fr<2) et aussi le cas critique (Fr=2)
#endif
#define EPS   0.005
#define K_DW  (5.*pi)
#define Q0    0.001
#define F_DW  0.048
#define S0    ((F_DW/8.)*sq(FR))

double h0, u0;

FILE *fp_profils = NULL;
FILE *fp_diag    = NULL;

int main()
{
  X0 = 0.;
  L0 = 2.;
  G  = 9.81;
  N  = 2500;    // Δx = 8e-4 m

  h0 = pow(Q0 / (sqrt(G) * FR), 2./3.);
  u0 = Q0 / h0;
  DT = 0.5 * (L0/N) / (u0 + sqrt(G * h0));

  char fname_profils[64], fname_diag[64];
  snprintf(fname_profils, sizeof(fname_profils), "out_dw_Fr%.4g.dat", (double)FR);
  snprintf(fname_diag,    sizeof(fname_diag),    "diag_Fr%.4g.dat",   (double)FR);

  fp_profils = fopen(fname_profils, "w");
  fp_diag    = fopen(fname_diag,    "w");

  if (!fp_profils || !fp_diag) {
    fprintf(stderr, "Erreur : impossible d'ouvrir les fichiers de sortie\n");
    return 1;
  }

  fprintf(fp_profils, "# x  h  u  zb  t\n");
  fprintf(fp_diag,    "# t  log(amplitude)\n");

  fprintf(stderr, "=== PARAMETRES ===\n");
  fprintf(stderr, "Fr = %.4f\n",    (double)FR);
  fprintf(stderr, "h0 = %.8f m\n",  h0);
  fprintf(stderr, "u0 = %.6f m/s\n", u0);
  fprintf(stderr, "S0 = %.6f\n",    (double)S0);
  fprintf(stderr, "DT = %.6f s\n",  DT);
  fprintf(stderr, "Dx = %.6f m\n",  L0/N);
  fprintf(stderr, "Check Fr = %.6f\n", u0/sqrt(G*h0));
  fprintf(stderr, "Critere : Fr %s 2 => %s\n",
          (double)FR > 2. ? ">" : "<",
          (double)FR > 2. ? "INSTABLE" : "STABLE");
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
    h[]   = h0 * (1. + EPS * sin(K_DW * x));
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
      double Cf = (F_DW / 8.) * fabs(u.x[]) / h[];
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

Pour Fr = 2.5, accord entre méthode numérique et solution analytique de Dressler.

~~~pythonplot Dressler (analytique) et résultat numérique
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.signal import find_peaks
from scipy.integrate import simpson

# ── 1. Paramètres Physiques ──────────────────────────────────────────────────
G = 9.81; f = 0.048; Fr = 2.5; Q0 = 0.001; L0 = 2.0
T_PLOT = 50.

m  = (f/8.) * Fr**2
H0 = (Q0 / (np.sqrt(G) * Fr))**(2./3.)
K  = np.sqrt(G) * H0**1.5

c  = (1. + np.sqrt(8.*m/f)) * np.sqrt(G * H0)

# ── 2. Lecture des données Basilisk ──────────────────────────────────────────
def read_profiles(fname):
    blocks = {}
    t_cur, xs_b, hs_b = None, [], []
    try:
        with open(fname) as fh:
            for line in fh:
                s = line.strip()
                if s.startswith("# t ="):
                    if t_cur is not None and xs_b:
                        blocks[t_cur] = (np.array(xs_b), np.array(hs_b))
                    t_cur = float(s.split("=")[1])
                    xs_b, hs_b = [], []
                elif s and not s.startswith("#"):
                    p = s.split()
                    if len(p) >= 2:
                        try: xs_b.append(float(p[0])); hs_b.append(float(p[1]))
                        except ValueError: pass
        if t_cur is not None and xs_b:
            blocks[t_cur] = (np.array(xs_b), np.array(hs_b))
    except FileNotFoundError: return None
    return blocks

blocks = read_profiles("out_dw_Fr2.5.dat")
if not blocks:
    raise FileNotFoundError("Aucun fichier de données Basilisk trouvé.")

times  = sorted(blocks.keys())
t_sel  = min(times, key=lambda t: abs(t - T_PLOT))
xs, hs = blocks[t_sel]
idx = np.argsort(xs); xs, hs = xs[idx], hs[idx]

# ── 3. Mesure de λ (Longueur d'onde) ─────────────────────────────────────────
peaks, _ = find_peaks(hs, prominence=0.1*(hs.max() - hs.min()))
N_WAVES  = len(peaks)
lam_meas = L0 / N_WAVES if N_WAVES > 0 else 0.4

# ── 4. Calcul des racines Dressler (Éq. 4.17) ────────────────────────────────
b   = H0 - (c**2 * f) / (8.*G*m)
ccc = (f * H0**2) / (8.*m)
D   = b**2 - 4.*ccc

if D < 0:
    raise ValueError(f"D < 0 ({D:.4g}). c={c} est trop élevé pour ces paramètres.")

sqD = np.sqrt(D)
HA  = max((-b + sqD)/2., (-b - sqD)/2.)
HB  = min((-b + sqD)/2., (-b - sqD)/2.)
Ac  = (HA**2 + HA*H0 + H0**2) / (m * (HA - HB))
Bc  = (HB**2 + HB*H0 + H0**2) / (m * (HA - HB))

def zeta(H):
    H = np.asarray(H, dtype=float)
    return ((H - H0)/m + Ac*np.log(np.abs((H - HA)/(H0 - HA))) 
            - Bc*np.log(np.abs((H - HB)/(H0 - HB))))

def Ha_of_Hb(Hb):
    return 0.5 * (np.sqrt(Hb**2 + 8.*K**2 / (G*Hb)) - Hb)

def lam_func(Hb):
    return zeta(Hb) - zeta(Ha_of_Hb(Hb))

# ── 5. Résolution de Hb ──────────────────────────────────────────────────────
Hb_arr = np.linspace(H0*1.001, H0*15, 5000)
diff   = np.array([lam_func(h) for h in Hb_arr]) - lam_meas
sc     = np.where(np.diff(np.sign(diff)))[0]

if len(sc) == 0:
    raise ValueError("Impossible de trouver Hb pour ce lambda. c est probablement trop élevé.")

Hb = brentq(lambda h: lam_func(h) - lam_meas, Hb_arr[sc[0]], Hb_arr[sc[0]+1])
Ha = Ha_of_Hb(Hb)

# ── 6. Construction du profil et Recalage Vertical ───────────────────────────
N_pts = 4000; eps = H0*1e-8
H1 = np.linspace(Ha, H0 - eps, N_pts)
H2 = np.linspace(H0 + eps, Hb, N_pts)
z_ref = float(zeta(Ha))
z_period = np.concatenate([zeta(H1), zeta(H2)]) - z_ref
H_period = np.concatenate([H1, H2])

# Calcul de la hauteur moyenne de Dressler pour corriger l'offset vertical
h_dressler_avg = simpson(y=H_period, x=z_period) / lam_meas
h_basilisk_avg = np.mean(hs)
v_offset = h_basilisk_avg - h_dressler_avg

print(f"H_avg Basilisk: {h_basilisk_avg:.6f}")
print(f"H_avg Dressler: {h_dressler_avg:.6f}")
print(f"Correction verticale appliquée: {v_offset:.6f}")

# Calage de phase
x_pic0 = xs[peaks[0]] if len(peaks) > 0 else 0.
x0 = (x_pic0 - lam_meas) % lam_meas

# ── 7. Tracé ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 5))
ax.plot(xs, hs, color="#e6b450", lw=2.5, label=f"Basilisk (Simulation)")

for k in range(-1, N_WAVES + 2):
    x_k = z_period + x0 + k*lam_meas
    mask = (x_k >= -0.1) & (x_k <= L0 + 0.1)
    # On applique v_offset pour que les masses correspondent
    ax.plot(x_k[mask], H_period[mask] + v_offset, color="#39d353", lw=2, ls="--",
            label="Dressler" if k==0 else "")
    
    xc = lam_meas + x0 + k*lam_meas
    if 0 <= xc <= L0:
        ax.vlines(xc, Ha + v_offset, Hb + v_offset, colors="#39d353", lw=2, ls="--")

ax.set_xlim(0, L0)
ax.set_ylim(min(hs)*0.95, max(hs)*1.05)
ax.set_xlabel("x [m]")
ax.set_ylabel("h [m]")
ax.set_title(f"Comparaison Fr={Fr}")
ax.legend()
ax.grid(alpha=0.2)
plt.tight_layout()
plt.savefig("dressler.png")

~~~    
*/


/**

Développement d'un train d'onde, l'échelle de la hauteur est amplifiée par rapport à la direction de l'écoulement.

~~~pythonplot Fr = 2.5
import numpy as np
import matplotlib.pyplot as plt
import os


fname = f"out_dw_Fr2.5.dat"

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
        ax.set_title(f"Roll-waves Darcy-Weisbach — Fr=2.5 — profils h(x)", fontsize=13)

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
        
        plt.savefig("profils_superposes_Fr2.5.png")
       
~~~    
*/

/**
## Bibliography
  
 * [J. Needham and J. H. Merkin](http://www.jstor.org/stable/2397934)
 "On roll waves down an open inclined channel"
 Proc. R. Soc. Lond. A 394, 259-278 (1984)
 
 * [Olivier Delestre](https://theses.hal.science/tel-00531377v3)
 "Simulation du ruissellement d’eau de pluie sur des surfaces agricoles"
 
 * [Lagrée P-Y](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
"Equations de Saint Venant et application, Ecoulements en milieux naturels" Cours MSF12, M1 Sorbonne-Université

*/