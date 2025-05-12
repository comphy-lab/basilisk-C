import numpy as np

sigma = np.sqrt(1e6)
mu_bubble = 1/1000; #mu'
mu_fluid = 1 #mu
mu_ratio = mu_fluid/mu_bubble
mu_ratio = mu_fluid/mu_bubble
Pe = 0.08269396*2e6
epsilon = 1e-2

B = -(3*sigma**6 - 3*sigma + 2*sigma**6*mu_ratio+3*sigma*mu_ratio)/ \
    (4*sigma**6-9*sigma**5+10*sigma**3-9*sigma+4*sigma**6*mu_ratio \
     -6*sigma**5*mu_ratio+6*sigma*mu_ratio+4*mu_ratio)
D = (2*B*sigma**3 + 4*B - 6*B*sigma-3*sigma)/ \
    (4*sigma**6-10*sigma**3+6*sigma)
C = 1/2-2*B/3/sigma-5*D*sigma**2/3
A = -B-C-D


#def d_psi_d_r(r, theta):
#    return np.sin(theta)**2*(-A/r**2+B+2*C*r+4*D*r**3)

#def d_psi_d_theta(r, theta):
#    return (A/r + B*r + C*r**2+D*r**4)*2*np.sin(theta)*np.cos(theta)

def u_r(r,theta):
    return -(A*r**(-3) + B*r**(-1) + C + D*r**2)*2*np.cos(theta)


def u_theta(r,theta):
    return np.sin(theta)*(-A*r**(-3)+B*r**(-1)+2*C+4*D*r**2)

from scipy.integrate import solve_bvp

r = np.linspace(1, 15, 1000000)

def fun(x, c):
    return np.vstack((c[1], Pe/2*u_r(x,0)*c[1] - 2/x*c[1]))

def bc(ya, yb):

    return np.array([ya[0]-1, yb[0]])


y = np.ones((2, r.size))
y[0,0] = 1
y[1,0] = -20

res = solve_bvp(fun, bc, r, y, max_nodes=1e8, verbose = 2)

print(res.y[1,0]*2)

import tqdm

def odefunc(c, theta, pbar, state):
    
    
    # state is a list containing last updated time t:
    # state = [last_t, dt]
    # I used a list because its values can be carried between function
    # calls throughout the ODE integration
    last_t, dt = state
    
    # let's subdivide t_span into 1000 parts
    # call update(n) here where n = (t - last_t) / dt
    n = int((theta - last_t)/dt)
    
    # we need this to take into account that n is a rounded number.
    state[0] = last_t + dt * n
    
    
    dcdtheta = np.zeros(X.shape)
    
    #BOUNDARY CONDITIONS
    dcdtheta[0] = 0 # constant at boundary condition mantains a dirichlet condition where the value equal to the initial value
    c[-1] = c[-2] #make outer boundary symmetric
    dcdtheta[-1] = 0
    
    
    #descretisation with linear upwind differencing for inner nodes, use first order upwind differencing for node 1, -2
    backward = (2/Pe*X[1]*(c[2]-2*c[1]+c[0])/h**2 + (4/Pe-X[1]*u_r(X[1],theta))*(c[1]-c[0])/h)/u_theta(X[1],theta)
    forward  = (2/Pe*X[1]*(c[2]-2*c[1]+c[0])/h**2 + (4/Pe-X[1]*u_r(X[1],theta))*(c[2]-c[1])/h)/u_theta(X[1],theta)
    dcdtheta[1] = ((4/Pe-X[1]*u_r(X[1],theta)) > 0)*forward + ((4/Pe-X[1]*u_r(X[1],theta)) < 0)*backward
    
    backward = (2/Pe*X[-2]*(c[-1]-2*c[-2]+c[-3])/h**2 + (4/Pe-X[-2]*u_r(X[-2],theta))*(c[-2]-c[-3])/h)/u_theta(X[-2],theta)
    forward  = (2/Pe*X[-2]*(c[-1]-2*c[-2]+c[-3])/h**2 + (4/Pe-X[-2]*u_r(X[-2],theta))*(c[-1]-c[-2])/h)/u_theta(X[-2],theta)
    dcdtheta[-2] = ((4/Pe-X[-2]*u_r(X[-2],theta)) > 0)*forward + ((4/Pe-X[-2]*u_r(X[-2],theta)) < 0)*backward
    
    backward = (2/Pe*X[2:-2]*(c[3:-1]-2*c[2:-2]+c[1:-3])/h**2 + (4/Pe-X[2:-2]*u_r(X[2:-2],theta))*(3*c[2:-2]-4*c[1:-3]+c[:-4])/(2*h))/u_theta(X[2:-2],theta)
    forward  = (2/Pe*X[2:-2]*(c[3:-1]-2*c[2:-2]+c[1:-3])/h**2 + (4/Pe-X[2:-2]*u_r(X[2:-2],theta))*(-c[4:] + 4*c[3:-1]-3*c[2:-2])/(2*h))/u_theta(X[2:-2],theta)
    dcdtheta[2:-2] = ((4/Pe-X[2:-2]*u_r(X[2:-2],theta)) > 0)*forward + ((4/Pe-X[2:-2]*u_r(X[2:-2],theta)) < 0)*backward

    pbar.set_postfix({'theta': theta, 'max': np.max(dcdtheta), 'min': np.min(dcdtheta)})
    pbar.update(n)
    #print(theta, np.max(dcdtheta), np.min(dcdtheta))
    
    return dcdtheta

from scipy.integrate import odeint


N = 20000  # number of points to discretize
start = 1
end = 8
X = np.linspace(start, end, N) # position along the rod
h = (end - start) / (N - 1)
init = res.sol(X)[0]

tspan = np.linspace(epsilon, np.pi-epsilon, 1001)



with tqdm.tqdm(total=1000, unit="‰",  position=0, leave=True) as pbar:
    sol = odeint(odefunc, 
                 init, 
                 tspan,
                 args = (pbar, [epsilon, (np.pi-2*epsilon)/1000])
                )

np.savetxt('Local_Sherwood.csv', sol, delimiter=",")