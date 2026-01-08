#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 20:28:46 2026
        (adapted from a previous program in Matlab)

Program designed to illustrate the chapter 5 of the course
  "Introduction to hydrodynamical instabilities"
  M2 DET, Toulouse

@author: David Fabre 

This program is protected by Creative Commons Licence CC-BY-NC
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
plt.close('all')

# Définition of the dynamical function for Lorenz system
def lorenz(XYZ, *, r=28, P=10,  beta=8/3.):
    X,Y,Z = XYZ
    X_dot = P * (Y - X)
    Y_dot = r * X - Y - X * Z
    Z_dot = X * Y - beta * Z
    return np.array([X_dot, Y_dot, Z_dot])

# Function for selection of parameters by user wirth default value
def myinput(question,default):
    answer = input(f"{question} [default : {default}] : " )
    if not answer:
        return default 
    else:
        return answer

# Selection or parameters by user     
print("\n Integration of the Lorenz system (demonstration program \n")
print(" Select the kind of computation you want to do : ")
print("       0  : Basic plot of the Lorenz attractor")
print("       1  : animation showing the trajectory and reconstruction of convection flow")
print("       2  : animation showing two trajectories initially close")
kind = int(myinput(" >>> Your choice ? ",0));

# Physical parameters
r = float(myinput("\n Enter the value of the control parameter r :",28))
P = 10;beta = 2.667;
print("    ( NB other parameters are set to P = "+str(P)+" ; beta  = "+str(beta)+" )")

# Conditions initiales (on peut les modifier)
print(" Select Initial condition :") 
Xinit = float(myinput("           X_init = ",0.5*np.sqrt(r)))
Yinit = float(myinput("           Y_init = ",-0.5*np.sqrt(r)))
Zinit = float(myinput("           Z_init = ",r/10))
init = Xinit,Yinit,Zinit

# Shema Runge-Kutta 4 used for time integration
def rk4_step(f, state, dt, **kwargs):
    # f : fonction dérivée (lorenz)
    k1 = f(state, **kwargs)
    k2 = f(state + 0.5*dt*k1, **kwargs)
    k3 = f(state + 0.5*dt*k2, **kwargs)
    k4 = f(state + dt*k3, **kwargs)
    return state + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)

# Paramètres de simulation
dt = 0.01           # pas de temps
num_steps = 10000   # nombre de pas

# Tableau pour stocker les valeurs
data = np.empty((num_steps + 2, 3))
data = np.zeros((num_steps, 3))
data[0] = init

# Intégration par RK4
for i in range(num_steps-1):
    data[i + 1] = rk4_step(lorenz, data[i], dt,r=r,P=P,beta=beta)
   
    
if (kind==0):
    # -------------------------
    # Version du programme qui trace juste l'attracteur
    # -------------------------
    
    # Création du graphique 3D
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    # Tracé de la trajectoire
    ax.plot(data[:,0],data[:,1],data[:,2], lw=0.5)  # xyzs.T transpose (3, N)
    ax.set_xlabel("X Axis")
    ax.set_ylabel("Y Axis")
    ax.set_zlabel("Z Axis")
    ax.set_title("Attracteur de Lorenz pour r = "+str(r))
    plt.show()


if(kind==2):
    # -------------------------
    # Version du programme qui trace deux trajectoires    
    # -------------------------
    
    print("\n For plot of two trajectories : distance between initial conditions") 
    epsilon = float(myinput("           epsilon = ",1e-3))
    
    # Integrate with RK4, second trajectory
    data2 = np.zeros((num_steps, 3))
    data2[0] = data[0] + [epsilon,epsilon,epsilon]
    for i in range(num_steps - 1):
        data2[i+1] = rk4_step(lorenz, data2[i], dt,r=r,P=P,beta=beta)
    
    # Prepare figure
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    line, = ax.plot([], [], [], lw=0.8,color='red')
    line2, = ax.plot([], [], [], lw=0.8,color='blue')
    ax.set_xlim((-20, 20))
    ax.set_ylim((-30, 30))
    ax.set_zlim((0, 50))  
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    
    # Animation function
    def update_2traj(frame):
        line.set_data(data[:frame, 0], data[:frame, 1])
        line.set_3d_properties(data[:frame, 2])
        line2.set_data(data2[:frame, 0], data2[:frame, 1])
        line2.set_3d_properties(data2[:frame, 2])
        
        return line,line2
    
    # Create animation
    ani = FuncAnimation(fig, update_2traj, frames=num_steps, interval=5, blit=True)
    plt.show()


if(kind==1):
    # -------------------------
    # Version du programme pour visualiser le champ de temp/vitesse
    # -------------------------
    
    # préparation maillage 2D pour contour + quiver
    n = 20  # résolution du quiver (peut être plus petit pour performance)
    x_line = np.linspace(0, 2, n)
    y_line = np.linspace(0, 1, n)
    x_grid, y_grid = np.meshgrid(x_line, y_line)
    
    def PlotConvection(XX,x_grid,y_grid):
        # Fonction qui sert a representer le champ de temperature et de vitesse 
        X,Y,Z = XX  
        ampX = .02;ampY = 0.02;ampZ = .005; pi = np.pi;     
        UU = -X*ampX*np.cos(y_grid*pi)*np.sin(x_grid*pi)
        VV = X*ampX*np.sin(y_grid*pi)*np.cos(x_grid*pi)
        TT = ((1-y_grid)*np.ones(x_grid.shape)) \
                    - ampZ*Z*(np.sin(2*pi*y_grid)*np.ones(x_grid.shape)) \
                    + ampY*Y*(np.sin(y_grid*pi)*np.cos(x_grid*pi))
        return TT,UU,VV
    
    # création figure + axes
    fig = plt.figure(figsize=(12,6))
    ax3d = fig.add_subplot(1,2,1, projection='3d')
    ax2d = fig.add_subplot(1,2,2)
    
    # axe 3D
    line3d, = ax3d.plot([], [], [], color='blue', lw=0.8)
    ax3d.set_xlim(-20, 20); ax3d.set_ylim(-30, 30); ax3d.set_zlim(0,50)
    ax3d.set_title("Attracteur de Lorenz 3D")
    ax3d.set_xlabel("X"); ax3d.set_ylabel("Y"); ax3d.set_zlabel("Z")
    
    # axe 2D
    TT,UU,VV= PlotConvection(init,x_grid,y_grid)
    Shape = np.ones_like(x_grid)
    contour_plot = ax2d.contourf(x_grid, y_grid, TT, levels=20, cmap='viridis')
    quiver_plot = ax2d.quiver(x_grid, y_grid, UU,VV,color="white")
    ax2d.set_title("Iso-contours T + quiver [U,V]")
    ax2d.set_xlabel("X"); ax2d.set_ylabel("Y")
    ax2d.axis('equal')
    cbar = fig.colorbar(contour_plot, ax=ax2d, shrink=0.8)
    
    # fonction d’update
    def update(frame):
        # mise à jour 3D
        line3d.set_data(data[:frame, 0], data[:frame, 1])
        line3d.set_3d_properties(data[:frame, 2])
        # Construction des champs à plotter
        TT,UU,VV= PlotConvection(data[frame, :],x_grid,y_grid)
        # suppression de l'ancien contourf
        if hasattr(ax2d, 'contour_plot'):
            for coll in ax2d.contour_plot.collections:
                coll.remove()
        # tracer le nouveau contourf
        ax2d.contour_plot = ax2d.contourf(x_grid, y_grid, TT, levels=20, cmap='viridis')
        # suppression de l'ancien quiver
        if hasattr(ax2d, 'quiver_plot'):
            ax2d.quiver_plot.remove()    
        # tracer le nouveau quiver
        ax2d.quiver_plot = ax2d.quiver(x_grid, y_grid, UU, VV, color='white', scale=np.sqrt(r)/2)
        return line3d, ax2d.contour_plot, ax2d.quiver_plot
    
    # animation
    ani = FuncAnimation(
        fig,
        update,
        frames=range(0, num_steps, 10),  # mise à jour tous les 10 pas
        interval=40,
        blit=False
    )
    
    plt.tight_layout()
    plt.show()



