/**
#   python Shallow Water with hydrostatic reconstruction for M1 TP 
#   same as https://basilisk.fr/sandbox/M1EMN/PYTHON/svdbWB.c 
#   note the nice table of graphs at the end
import numpy as np
import matplotlib.pyplot as plt

# paramètres
dt = 0.01
tmax = 32
dx = 0.025
nx = 2000
Cf = 0.0
h0 = 2
alpha = .1
xbosse = 0.0
xtas = -2

times_to_save = np.arange(0, 64, 1)  # 
saved_solutions = []
saved_times = []

# tableaux
x = np.zeros(nx+2)
h = np.zeros(nx+2)
u = np.zeros(nx+2)
Q = np.zeros(nx+2)
Z = np.zeros(nx+2)
dZ = np.zeros(nx+2)

fp = np.zeros(nx+2)
fd = np.zeros(nx+2)

hg = np.zeros(nx+2)
hd = np.zeros(nx+2)
hig = np.zeros(nx+2)
hid = np.zeros(nx+2)
ug = np.zeros(nx+2)
ud = np.zeros(nx+2)

hn = np.zeros(nx+2)
un = np.zeros(nx+2)

# reconstruction hydrostatique
def reconsetat(hl, hr, dz):
    hil = max(0.0, hl - max(0.0, dz))
    hir = max(0.0, hr - max(0.0, -dz))
    return hil, hir


# flux HLL pour h
def FHLL1(ug, ud, hg, hd):
    c1 = min(ug - np.sqrt(hg), ud - np.sqrt(hd))
    c2 = max(ug + np.sqrt(hg), ud + np.sqrt(hd))

    if c1 >= 0:
        f1 = hg * ug
    elif c1 < 0 and 0 < c2:
        f1 = (c2*hg*ug - c1*hd*ud)/(c2-c1) + c1*c2*(hd-hg)/(c2-c1)
    else:
        f1 = hd * ud
    return f1

# flux HLL pour hu
def FHLL2(ug, ud, hg, hd):
    c1 = min(ug - np.sqrt(hg), ud - np.sqrt(hd))
    c2 = max(ug + np.sqrt(hg), ud + np.sqrt(hd))

    if c1 >= 0:
        f2 = hg*ug*ug + hg*hg/2
    elif c1 < 0 and 0 < c2:
        f2 = (c2*(hg*ug*ug+hg*hg/2) - c1*(hd*ud*ud+hd*hd/2))/(c2-c1) \
             + c1*c2*(hd*ud-hg*ug)/(c2-c1)
    else:
        f2 = hd*ud*ud + hd*hd/2
    return f2


# initialisation
for i in range(nx+2):
    x[i] = (i - nx/2) * dx
    if x[i] > xbosse:
        Z[i] = alpha * (x[i] - xbosse)
    if x[i] < xtas:
        h[i] = h0

for i in range(nx+1):
    dZ[i] = Z[i+1] - Z[i]

t = 0
it = 0

while t <= tmax:

    t += dt
    it += 1

    hd[:] = h
    hg[:] = h
    ud[:] = u
    ug[:] = u

    # reconstruction
    for i in range(1, nx+1):
        hid[i-1], hig[i] = reconsetat(hd[i-1], hg[i], dZ[i-1])

    # flux de h et q=hu
    cfl2 = 0
    for i in range(1, nx+1):
        fp[i] = FHLL1(ud[i-1], ug[i], hid[i-1], hig[i])
        fd[i] = FHLL2(ud[i-1], ug[i], hid[i-1], hig[i])
 
    # mise à jour de h et q=hu
    for i in range(1, nx):

        hn[i] = h[i] - dt*(fp[i+1]-fp[i])/dx

        if h[i] > 0:
            q = h[i]*u[i] - dt*(fd[i+1]-fd[i])/dx \
                - dt*(hig[i]**2/2 - hg[i]**2/2
                      + hd[i]**2/2 - hid[i]**2/2)/dx

            if q != 0 and Cf > 0 :
                q = q *h[i]**2/(h[i]**2 + dt*Cf)

            un[i] = q/hn[i]
        else:
            un[i] = 0

    # conditions aux limites, cellules fantomes 
    hn[0] = hn[1]
    un[0] = -un[1]

    hn[nx+1] = hn[nx]
    un[nx+1] = un[nx]

    # mise à jour des champs
    for i in range(nx+2):
        h[i] = hn[i]
        u[i] = un[i]
        Q[i] = hn[i]*un[i]

    # sauvegarde périodique
    if it % 5 == 0:
        np.savetxt(
            "solpy.OUT",
            np.column_stack((x[:nx], h[:nx], Q[:nx], Z[:nx]))
        )
    # sauvegarde pour tracé futur
    if len(saved_times) < len(times_to_save):
        if t >= times_to_save[len(saved_times)]:
            saved_solutions.append(h.copy())
            saved_times.append(t)
            print( it, " t =", t)
            
# tracé final sur une matrice de figures, merci le chat qui p
fig, axes = plt.subplots(6, 5, figsize=(10, 10))

for k, ax in enumerate(axes.flat):

    if k < len(saved_solutions):

        h_plot = saved_solutions[k]

        ax.plot(x[:nx], h_plot[:nx] + Z[:nx])
        ax.plot(x[:nx], Z[:nx])
        ax.set_ylim(0, 4) 
        ax.set_title(f"t={k}",loc='left',fontsize=6)
        ax.tick_params(axis='both', labelsize=6)
        ax.set_xlabel("x",fontsize = 5)
        ax.set_ylabel("h+Z, h",fontsize = 5)

plt.subplots_adjust(
    left=0.06,
    right=0.98,
    bottom=0.06,
    top=0.95,
    wspace=0.25,
    hspace=0.35
)
plt.savefig("resultats_dambreak.pdf")
plt.show()

*/