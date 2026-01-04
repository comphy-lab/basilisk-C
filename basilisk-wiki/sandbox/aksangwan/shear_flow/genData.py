# ------------------------------------------------------------------------------
# Script for Generating Initial Conditions for Two-Phase Shear Flow Simulation
# ------------------------------------------------------------------------------
# This Python script is used to generate the initial velocity field and interface
# shape for the two-phase shear flow simulation in the Basilisk C framework.
# The script calculates the initial conditions based on wave parameters such as 
# amplitude, wavenumber, wave speed, and other simulation parameters.
#
# Key Features:
# - Generates velocity data for both the x and y components (u.x and u.y).
# - Computes interface shape using linear stability analysis with a given wave 
#   number (k), amplitude (a), and other physical parameters.
# - Saves the data to be used by the main simulation code (`single.c`).
# - Can output the velocity and interface data to a `.dat` file for use in 
#   the simulation.
#
# Dependencies:
# - `numpy` for numerical computations.
# - `scipy` for interpolation and mathematical functions.
# - `matplotlib` for plotting the initial interface (optional).
#
# Usage:
# - To run the script, use:
#     python3 genData.py -a  -cr <real_wave_speed> -ci <imaginary_wave_speed> 
#                         -k  -U0  -g  -d <density_ratio> 
#                         -ha  -eta
#
# Outputs:
# - `veldata.dat`: Contains the velocity field in two columns (x and y components).
# - `eta.dat`: Contains the initial interface height profile.
#
# Example command:
#     python3 genData.py -a 0.02 -cr 61.25600429 -ci 22.55092783 -k 0.2920977808 
#                         -U0 199.5675324 -g 980.0 -d 0.5 -ha 2.54 -eta
#
# Author: A. K. Sangwan
# Date: 2025
# ------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import sympy as sm
import scipy as sc
from typing import Tuple
import cmath as cm
from argparse import ArgumentParser


def arguments():
    parser = ArgumentParser()
    parser.add_argument('-a', type=float)
    parser.add_argument('-cr', type=float)
    parser.add_argument('-ci', type=float)
    parser.add_argument('-k', type=float)
    parser.add_argument('-H', type=float)
    parser.add_argument('-ha', type=float)
    parser.add_argument('-U0', type=float)
    parser.add_argument('-g', type=float)
    parser.add_argument('-d', type=float)
    parser.add_argument('-vel', action='store_true')
    parser.add_argument('-eta', action='store_true')
    return parser.parse_args()

def vel_field(y, params, d=0):
    U0 = params['U0']
    Us = params['Us']
    ha = params['ha']
    hw = params['hw']
    if not isinstance(y, np.ndarray):
        y = np.array([y])
    x = sm.symbols('x')
    Ua_expr = sm.diff(U0 - (U0 - Us) * sm.exp(-x / ha), (x, d))
    Uw_expr = sm.diff(Us * sm.exp(x / hw), (x, d))
    Ua = sm.lambdify(x, Ua_expr, 'numpy')
    Uw = sm.lambdify(x, Uw_expr, 'numpy')
    return np.where(y >= 0, Ua(y), Uw(y))

def velocity(x, c: complex, k: float, params: dict) -> Tuple[np.ndarray, np.ndarray]:
    U0 = params['U0']
    Us = params['Us']
    ha = params['ha']
    hw = params['hw']

    Aa = (U0-Us)/(U0-c)
    a = params['a']
    Aw = Us/c
    ka = k*ha
    kw = k*hw
    aa = ka - np.sqrt(1+ka**2)
    aw = kw - np.sqrt(1+kw**2)
    ba = ka + np.sqrt(1+ka**2)
    bw = kw + np.sqrt(1+kw**2)
    ua = (a*(c-Us)/sc.special.hyp2f1(aa, ba, 2*ka+1, Aa)
          * (-k*np.exp(-k*x)*sc.special.hyp2f1(aa, ba, 2*ka+1, Aa*np.exp(-x/ha))
             + (Aa*np.exp(-x/ha - k*x)/(ha*(2*ka+1))
                * sc.special.hyp2f1(aa+1, ba+1, 2*ka+2, Aa*np.exp(-x/ha)))))
    uw = (a*(c-Us)/sc.special.hyp2f1(aw, bw, 2*kw+1, Aw)
          * (k*np.exp(k*x)*sc.special.hyp2f1(aw, bw, 2*kw+1, Aw*np.exp(x/hw))
             - (Aw*np.exp(x/hw + k*x)/(hw*(2*kw+1))
                * sc.special.hyp2f1(aw+1, bw+1, 2*kw+2, Aw*np.exp(x/hw)))))
    va = (-1j*a*k*(c-Us)*np.exp(-k*x)/sc.special.hyp2f1(aa, ba, 2*ka+1, Aa)
          * sc.special.hyp2f1(aa, ba, 2*ka+1, Aa*np.exp(-x/ha)))
    vw = (-1j*a*k*(c-Us)*np.exp(k*x)/sc.special.hyp2f1(aw, bw, 2*kw+1, Aw)
          * sc.special.hyp2f1(aw, bw, 2*kw+1, Aw*np.exp(x/hw)))
    return ua, uw, va, vw

#   np.where(y >= 0, ua, uw), np.where(y >= 0, va, vw)

def main():
    arg = arguments()
    # c = np.array([float(i) for i in arg.c.strip('][').split(',')])
    c = arg.cr + 1j*arg.ci #complex(arg.cr, arg.ci)
    k = arg.k
    h = arg.H
    A = arg.a
    vo = arg.U0
    g = arg.g
    b = arg.ha
    params = {
        'U0': vo,
        'Us': 0.0,
        'ha': b,
        'hw': 0.1,
        'delta': arg.d,
        'a': A,
        'g': arg.g,
        'T': 72.0
    }
    if arg.vel:
        data = np.loadtxt("data.csv", delimiter=",", dtype=float)
        velocities = np.zeros((data.shape[0], 2))
        # data[data[:,2] >= 0.5,2] = 1.0
        # data[data[:,2] < 0.5,2] = 0.0
        air_u, water_u, air_w, water_w = velocity( data[:, 1], c, k, params)
        velocities[:, 0] = ((water_u*data[:, 2] 
                            + air_u*(1-data[:, 2]))*np.exp(1j*k*data[:, 0]) 
                            + vel_field(data[:, 1], params)).real
        velocities[:, 1] = ((water_w*data[:, 2] 
                            + air_w*(1-data[:, 2]))
                            *np.exp(1j*k*data[:, 0])).real
        np.savetxt("veldata.dat", velocities, delimiter=",")
    
    if arg.eta:
        data = np.loadtxt("etadata.csv", delimiter=",", dtype=float)
        data = np.unique(data)
        eta = np.zeros((data.shape[0], 2))
        eta[:, 0] = data
        eta[:, 1] = (arg.a*np.exp(1j*k*data)).real
        # eta = np.round(eta, decimals=8)
        np.savetxt("eta.dat", eta, delimiter="\t")


if __name__ == "__main__":
    main()
