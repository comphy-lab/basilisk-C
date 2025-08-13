"""
# Compute the theoretical solutions of 3D oscillating drop/bubble

This script computes various theoretical solutions for a 3D
oscillating drop/bubble problem, including the simplest [Lamb solution](#lamb1932) 
based on normal-mode analysis, and the [Prosperetti solution](#prosperetti1980)
based on solving the initial-value problem.

"""

import numpy as np
import matplotlib.pyplot as plt
from mpmath import invertlaplace, sqrt, besseli, besselk
from math import ceil

""" 
## Lamb solution
"""

class SolLamb():
    def __init__(self, **kwargs):
        self.r0 = kwargs.get("r0", 1.)
        self.a0 = kwargs.get("a0", 0.025)
        self.rho_l = kwargs.get("rho_l", 10.)
        self.rho_g = kwargs.get("rho_g", 1.)
        self.mu_l = kwargs.get("mu_l", 5.e-2)
        self.mu_g = kwargs.get("mu_g", 5.e-4)
        self.sigma = kwargs.get("sigma", 10.)
        self.drop = kwargs.get("drop", True)
        self.nu_l = self.mu_l / self.rho_l
        self.no = kwargs.get("no", 2)
        self.da0 = kwargs.get("da0", 0.)

        self.omega2 = 0.
        self.omega = 0.
        self.b0 = 0.

        if self.drop:
            self.rho_in, self.rho_out = self.rho_l, self.rho_g
            self.mu_in, self.mu_out = self.mu_l, self.mu_g
        else:
            self.rho_in, self.rho_out = self.rho_g, self.rho_l
            self.mu_in, self.mu_out = self.mu_g, self.mu_l

        self.nu_in = self.mu_in / self.rho_in
        self.nu_out = self.mu_out / self.rho_out

        self._set_coef(self.no)
        self.La = 2. * self.rho_l * self.sigma * self.r0 / self.mu_l**2
        print("Laplace Number:{:f}".format(self.La))

    def _set_coef(self, no):
        # p. 475 in Lamb's book, Eq. (10)
        # P. 641 in Lamb's book, Eqs. (12) (13)
        self.no = no
        if self.drop:
            self.b0 = (no - 1) * (2 * no + 1) * self.mu_l \
                      / self.rho_l / self.r0**2
        else:
            self.b0 = (no + 2) * (2 * no + 1) * self.mu_l \
                      / self.rho_l / self.r0**2

        self.omega2 = no * (no - 1) * (no + 1) * (no + 2) * self.sigma \
                      / ((no + 1) * self.rho_in + no * self.rho_out) / self.r0**3
        self.omega = np.sqrt(self.omega2)
        self.beta0 = no * (no + 2) / (2 * no + 1)

    def get_solution(self, no, te, ts=0., nt=10):
        self._set_coef(no)
        if isinstance(te, (np.ndarray, list, tuple)):
            t_s = np.array(te)
        else:
            t_s = np.linspace(ts, te, nt)
        self.t = t_s
        self.a_t = np.exp(-self.b0 * t_s) * np.cos(self.omega * t_s) * self.a0
        return self.t, self.a_t


""" 
## Prosperetti solution with free-surface assumption

"""

class SolProsperetti(SolLamb):
    def __init__(self, **kwargs):
        super(SolProsperetti, self).__init__(**kwargs)

    def _set_coef(self, no):
        # Prosperetti (1980) Eqs.(3) - (8)
        self.no = no
        if self.drop:
            self.omega2 = no * (no - 1) * (no + 2) * self.sigma \
                          / (self.rho_l * self.r0**3)
            self.b0 = (no - 1) * (2 * no + 1) * self.mu_l \
                      / self.rho_l / self.r0**2
            self.beta0 = (no - 1) * (no + 1) / (2 * no + 1)
        else:
            self.omega2 = (no - 1) * (no + 1) * (no + 2) * self.sigma \
                          / (self.rho_l * self.r0**3)
            self.b0 = (no + 2) * (2 * no + 1) * self.mu_l \
                      / self.rho_l / self.r0**2
            self.beta0 = no * (no + 2) / (2 * no + 1)

        self.omega = np.sqrt(self.omega2)

    def _get_f_ap(self, a=1.):
        # Prosperetti (1980) Eq.(13).
        def f_ap(p):
            p /= a
            q = sqrt(p / self.nu_l) * self.r0
            if self.drop:
                gamma_f = q * besseli(self.no + 0.5, q) / besseli(self.no + 1.5, q)
                qp = 1. / (1 - 0.5 * gamma_f)
            else:
                gamma_f = q * besselk(self.no + 0.5, q) / besselk(self.no - 0.5, q)
                qp = -1. / (1 + 0.5 * gamma_f)

            cn = p * self.da0 - self.omega2 * self.a0
            cd = p**2 + 2. * self.b0 * p + self.omega2 \
                 + 2. * self.beta0 * self.b0 * qp * p
            return (self.a0 + cn / cd) / p / a
        return f_ap

    def get_solution(self, no, te, nt=10, scale=False):
        self._set_coef(no)

        if isinstance(te, (np.ndarray, list, tuple)):
            t_s = np.array(te)
        else:
            t_s = np.linspace(0., te, nt)

        # scaling the time does not change anything.
        if scale:
            ca = (2. * np.pi) / self.omega
        else:
            ca = 1.

        f_p = self._get_f_ap(ca)
        a_t = []
        for _t in t_s:
            _t /= ca
            if _t < 1.e-6:
                a_t.append(self.a0)
            else:
                a = invertlaplace(f_p, _t, method="cohen")
                a_t.append(a)
        self.t = t_s
        self.a_t = np.array(a_t)
        return self.t, self.a_t


""" 
## Prosperetti solution with the effects of both fluids.
"""

class SolProsperettiTP(SolProsperetti):
    def _set_coef(self, no):
        
        # Prosperetti (1980) Eqs.(10) - (12), Eq.(14).
        self.no = no
        self.omega2 = no * (no - 1) * (no + 1) * (no + 2) * self.sigma \
                      / ((no + 1) * self.rho_in + no * self.rho_out) / self.r0**3
        self.omega = np.sqrt(self.omega2)

    def _get_f_ap(self, a=1.):
        # Prosperetti (1980) Eq.(14).
        def f_ap(p):
            _no = self.no
            p /= a

            q_in = sqrt(p / self.nu_in) * self.r0
            q_out = sqrt(p / self.nu_out) * self.r0
            gamma_n = _no * self.rho_out + (_no + 1) * self.rho_in

            gamma_in = q_in * besseli(_no + 0.5, q_in) \
                       / besseli(_no + 1.5, q_in)
            gamma_out = q_out * besselk(_no + 0.5, q_out) \
                        / besselk(_no - 0.5, q_out)

            dmu = self.mu_out - self.mu_in

            c_in = (2. * _no + 1) * self.mu_in * gamma_in \
                   + 2. * _no * (_no + 2) * dmu
            c_out = (2. * _no + 1) * self.mu_out * gamma_out \
                    - 2. * (_no - 1)* (_no + 1) * dmu
            dp = c_in * c_out / \
                 (self.mu_out * gamma_out + self.mu_in * gamma_in + 2. * dmu)

            cn = p * self.da0 - self.omega2 * self.a0
            cd = p**2 + self.omega2 + dp * p / gamma_n
            return (self.a0 + cn / cd) / p / a
        return f_ap

"""
## Example
"""

if __name__ == "__main__":
    fontsize_legend = 15
    fontsize_label = 15
    fontsize_title = 15

    save_data = not True
    save_figs = not True

    # parameters from my simulation
    # $$ \textrm{Oh} = \frac{\mu_l}{\rho_l \sigma D} $$
    Oh = 0.01
    # amplitude
    eps = 0.01
    para = {"r0": 1., "a0": eps, "rho_l": 1., "rho_g": 1.,
            "mu_l": Oh * 2**0.5, "mu_g": Oh * 2**0.5, "sigma": 1.}

    # droplet case
    lam_sol = SolLamb(**para)
    pro_sol = SolProsperetti(**para)
    pro_tp_sol = SolProsperettiTP(**para)

    # bubble case
    lam_bub_sol = SolLamb(**para, drop=False)
    pro_bub_sol = SolProsperetti(**para, drop=False)
    pro_bub_tp_sol = SolProsperettiTP(**para, drop=False)

    per = 2. * np.pi / lam_sol.omega
    per_bub = 2. * np.pi / lam_bub_sol.omega

    # number of periods
    n_per = 5
    # sample points within one period
    nt = 32 * n_per

    # sample points
    t_sample = t_s = np.linspace(0., n_per * per, nt)

    t, a = lam_sol.get_solution(2, t_sample)
    t_bub, a_bub = lam_bub_sol.get_solution(2, t_sample)

    t_p, a_p = pro_sol.get_solution(2, t_sample)
    t_p_bub, a_p_bub = pro_bub_sol.get_solution(2, t_sample)

    t_tp, a_tp = pro_tp_sol.get_solution(2, t_sample)
    t_tp_bub, a_tp_bub = pro_bub_tp_sol.get_solution(2, t_sample)

    La = ceil(lam_sol.La)

    # save the figures
    for i in range(2):
        fig, ax = plt.subplots(1, 1)
        if i == 0:
            ax.plot(t / per, a / eps, "r-")
            ax.plot(t_p / per, a_p / eps, "b--")
            ax.plot(t_tp / per, a_tp / eps, "g-.")
            case_name = "droplet"
        else:
            ax.plot(t_bub / per_bub, a_bub / eps, "r-")
            ax.plot(t_p_bub / per_bub, a_p_bub / eps, "b--")
            ax.plot(t_tp_bub / per_bub, a_tp_bub / eps, "g-.")
            case_name = "bubble"

        ax.set_xlim((0, n_per))
        ax.set_ylim((-1.2, 1.2))
        ax.set_ylabel(r"$a(t) / a_0$", fontsize=fontsize_label)
        ax.set_xlabel(r"$t / T_0$", fontsize=fontsize_label)
        ax.tick_params(labelsize=fontsize_label)
        ax.set_title("{:s}, La = {:d}".format(case_name, La))
        ax.legend(["Lamb", "Prosperetti-FS", "Prosperetti"],
                   fontsize=fontsize_legend, loc="upper right")
        ax.grid(which="both")
        fig.set_tight_layout(True)

        if save_figs:
            path = "{:s}_ana_{:d}.pdf".format(case_name, La)
            fig.savefig(path)

    # output the solutions
    if save_data:
        # Three solutions: Lamb, Prosperetti free-surface, Prosperetti solutions
        path = "droplet_ana_{:d}.dat".format(La)
        with open(path, mode="w") as f:
            for ind, _t in enumerate(t):
                f.write("{:.8e} {:.8e} {:.8e} {:.8e}\n".
                        format(_t, float(a[ind]), float(a_p[ind]), float(a_tp[ind])))

        path = "bubble_ana_{:d}.dat".format(La)
        with open(path, mode="w") as f:
            for ind, _t in enumerate(t_bub):
                f.write("{:.8e} {:.8e} {:.8e} {:.8e}\n".
                        format(_t, float(a_bub[ind]), float(a_p_bub[ind]), float(a_tp[ind])))

    plt.show()

"""
## References

~~~bib
@Book{lamb1932,
  title={Hydrodynamics},
  author={Lamb, Horace},
  year={1932},
  publisher={Dover}
}

@article{prosperetti1980,
  title={Free oscillations of drops and bubbles: the initial-value problem},
  author={Andrea Prosperetti},
  journal={J. Fluid Mech.},
  volume={100},
  pages={333-347},
  year={1980},
  publisher={Cambridge University Press},
  doi={https://doi.org/10.1017/S0022112080001188}
}
~~~
"""
