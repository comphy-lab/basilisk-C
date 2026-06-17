r"""
An analytical solution can be found

    $\overline{u}(\sigma,t) = B_n(\sigma,t) \sum_{j=0}^{\infty} \frac{A_n^{j+1}(\sigma,t)}{j!}
                                                I_2(j+1, \Theta_1, \Theta_2)$

with $B_n=\frac{\omega}{k_n^{2}\sigma L} e^{- \sigma H_0 k}$, $A_n=\sigma
ak_n$, $I_2$ is computed using the reduction relation (integration by part):

    I_2(m)= \frac{1}{m} \cos^{m-1} \theta \sin \theta + \frac{m-1}{m}I_2(m-2)
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.integrate as integrate
g=9.81
'''
The functions to integrate numerically
'''
def f_eta(x, k, a):
    return a*np.sin(k*x)

def f_u(x, z, k, a):
    om = np.sqrt(g*k)
    return a*om*np.exp(k*z)*np.sin(k*x)

def f_ul(x, d, k, a):    
    return f_u(x, f_eta(x, k, a)-d, k, a)
'''
The function needed for the analytical solution
'''
def first_term_reduction(m, theta):
    return np.cos(theta)**(m-1)*np.sin(theta)
'''
The reduction of cos^n
'''
def reduction_cosn(n, theta1, theta2):
    result = 0.
    if n//2==0:
        m=0
    else:
        m=1
    if n==1:
        m=1
    while m<=n:
        if m==0:
            result += theta2-theta1
        elif m==1:
            result += np.sin(theta2) - np.sin(theta1)
        else:
            # if m!=n:
            #     result = (m+1)/(m+2) * result # here m is already the next iteration ! m+2-1=m+1
            result +=( 
                    1/m*( first_term_reduction(m, theta2)
                        - first_term_reduction(m, theta1) )
                    + (m-1)/m * result 
                    )
            
        m = m+2
    return result

def f_Bn(sigma, t, omega, kn, H0, L):
    if sigma <=0 or sigma > 1:
        raise Exception(f'sigma is defined on ]0:1], it is sigma={sigma}')
    return omega/(kn**2*sigma*L)*np.exp(-sigma*H0*kn)

def f_An(sigma, a, kn):
    if sigma <=0 or sigma > 1:
        raise Exception(f'sigma is defined on ]0:1], it is sigma={sigma}')
    return sigma*a*kn

def analytical_u_profile_mode_m(n, sigma, t, omega, m, H0, L, a):
    '''
    Integrate on the segment [0,L]. Sum of first 'n' terms, for mode 'm'
    '''
    km = 2*m*np.pi/L
    Bn = f_Bn(sigma, t, omega, km, H0, L)
    An = f_An(sigma, a, km)
    result = 0
    tht1 = -omega*t
    tht2 = km*L - omega*t
    for j in range(n):
        fact = np.prod(np.arange(1, j+1)) # +1 is to take into account j
        result += An**(j+1)/fact * reduction_cosn(j+1, tht1, tht2)
    result = Bn*result
    return result

"""
A function to plot a verification test that the serie converges to the 
numerically integrated function
"""

def verif_first_terms_recursion():
    tht1 = 0
    tht2 = 2*np.pi
    

    true_terms = np.zeros(4)

    # m = 0
    true_terms[0] = tht2-tht1
    # m = 1
    true_terms[1] = np.sin(tht2)-np.sin(tht1)
    # m = 2
    true_terms[2] = ( (np.cos(tht2)*np.sin(tht2)
                        - np.cos(tht1)*np.sin(tht1))
                    + 1/2*(tht2-tht1) )
    # m = 3
    true_terms[3] = ( 
                      1/3*(np.cos(tht2)**2*np.sin(tht2)
                           - np.cos(tht1)**2*np.sin(tht1))
                    + 2/3*( true_terms[2] ) 
                     )


    terms = np.zeros(4)
    for k in range(len(terms)):
        terms[k] = reduction_cosn(k, tht1, tht2)
    
    for k in range(4): # len(terms)
        print(f'k={k} true = {true_terms[k]}, term = {terms[k]}')
    

def verif_analytical_solution(nmax):
    dzl = 0.05
    H0 = 100        # m, depth of water
    a = 1           # m, amplitude of the wave
    L = 200.        # m, wavelength
    g = 9.81        # m.s-2
    k = 2*np.pi/L   # m-1, wavenumber
    
    omega = np.sqrt(g*k) # linear dispersion relation

    # 1) numerical integration
    # Integration over 1 wavelength for a few layers
    # (see "Layer definitions")
    myX = np.arange(-L/2,L/2,L/1000)
    total_H = f_eta(myX, k, a) + H0
    H_frac = np.arange(0,1,0.05)
    z_layers = np.zeros((len(myX),len(H_frac)))
    # -> building an array with layer height(x)
    for l in range(len(H_frac)):
        z_layers[:,l] = total_H*H_frac[l]
    a_profile = np.zeros(len(H_frac))
    # -> we integrate
    borne1 = 0.
    borne2 = L
    for l in range(len(H_frac)):
        result, err = integrate.quad(lambda x: f_ul(x,d=(f_eta(x,k,a)+H0)*H_frac[l], k=k, a=a),
                                     borne1, borne2)
        result = 1/L * result
        if err>1e-6:
            print(result, err)
            raise Exception('Error in integration of function')
        a_profile[l] = result

    # 2) Analytical solution
    zl = np.arange(dzl,1.,dzl)
    a2_profile = np.zeros(zl.shape)
    cmap = mpl.colormaps['plasma']
    colors = cmap(np.linspace(0, 1, nmax))
    fig, ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=100)
    for k in range(nmax):
        profile = np.zeros(zl.shape)
        for kz in range(len(zl)):
            a2_profile[kz] = analytical_u_profile_mode_m(n=k, 
                                                    sigma=zl[kz], 
                                                    t=0., 
                                                    omega=omega, 
                                                    m=1, 
                                                    H0=H0, 
                                                    L=L, 
                                                    a=a)
        ax.plot(a2_profile, zl, label=f'n={k}', c=colors[k])
    ax.plot(a_profile, H_frac, label='integ num')
    ax.set_xlabel(r'$\overline{u}$')
    ax.set_ylim([1,0])
    ax.set_ylabel('layer index')
    ax.legend()
    
