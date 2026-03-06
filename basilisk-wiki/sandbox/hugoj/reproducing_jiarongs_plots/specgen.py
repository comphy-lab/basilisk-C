import numpy as np
from scipy.interpolate import griddata


# Function to convert polar to cartesian and interpolate
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

""" Separate functions that use different omni-directional spectra """
def spectrum_PM (P, kp, kmod):
    F_kmod = P*kmod**(-2.5)*np.exp(-1.25*(kp/kmod)**2)
    return F_kmod

def spectrum_JONSWAP (alpha, kp, kmod):
    F_kmod = alpha*kmod**(-3)*np.exp(-1.25*(kp/kmod)**2)
    return F_kmod

def spectrum_Gaussian (G, span, kp, kmod):
    F_kmod = (G/span)*np.exp(-0.5*((kmod-kp)**2/span**2))
    return F_kmod
 
""" Generate linealy spaced kx, ky according to specified # of modes, and interpolated F(kx,ky) """
def spectrum_gen_linear(shape, N_mode=32, N_power=5, L=200):
    
    """ The function to generate a kx-ky spectrum based on uni-directional spectrum and a 
        cos^N (theta) directional spreading.
        Arguments: 
            shape: the spectrum shape (a function)
            N_mode: # of modes (Right now it's N_mode for kx and N_mode+1 for ky; 
                    has to much what's hard-coded in the spectrum.h headerfile.
            L: physical domain size
            """
    
    """ # of modes for the uni-directional spectrum 
        (doesn't matter as much because of interpolation anyway) """
    N_kmod = 64; N_theta = 64 # Uniform grid in kmod and ktheta, can be finer than N_mode 
    thetam = 0 # midline direction
    kmod = np.linspace(2*np.pi/L,1.41*100*2*np.pi/L,N_kmod) # Change 
    theta = np.linspace(-0.5*np.pi, 0.5*np.pi, N_theta) + thetam # Centered around thetam
    kmod_tile, theta_tile = np.meshgrid(kmod,theta)

    """ Pick the spectrum shape """
    F_kmod = shape (kmod) # includes the spectral shape, peak, and energy level of choice
    D_theta = np.abs(np.cos(theta-thetam)**N_power) 
    dtheta = theta[1]-theta[0]
    D_theta = D_theta/np.trapz(D_theta, theta)  # Normalize so the sum equals one
    F_kmod_tile, D_theta_tile = np.meshgrid(F_kmod,D_theta) 
    F_kmodtheta_tile = F_kmod_tile*D_theta_tile/kmod_tile # Notice!! Normalize by k
    
    """ Uniform grid in kx,ky """
    kx = np.arange(1,N_mode+1)*2*np.pi/L # based on the grid, interval can't go smaller then pi/L
    ky = np.arange(-N_mode/2,N_mode/2+1)*2*np.pi/L
    kx_tile, ky_tile = np.meshgrid(kx,ky)
    kxp_tile, kyp_tile = pol2cart(kmod_tile, theta_tile)
    
    """ Project from uniform k to uniform kx,ky """
    F_kxky_tile = griddata((kxp_tile.ravel(), kyp_tile.ravel()), F_kmodtheta_tile.ravel(), 
                           (kx_tile, ky_tile), method='linear', fill_value=0) 
    
    return kmod, F_kmod, kx, ky, F_kxky_tile

""" Generate log spaced kx, ky according to specified # of modes, and interpolated F(kx,ky).
    This has not been implemented yet and ideally kx and ky should also be logarithmically spaced.
    But then we need to make modification to the c code as well (dkx and dky not uniform anymore).
    Also right now this is hard-coded PM spectrum """

# def spectrum_gen_log(kp=2*np.pi/100, P=0.45, L=500, N_power=5):
#     thetam = 0
#     N_mode = 32; # Uniform grid in kx and ky : how many we want to put in the initialization
#     N_kmod = 64; N_theta = 64 # Uniform grid in kmod and ktheta
#     kmod = np.logspace(np.log(2*np.pi/L),np.log(200*2*np.pi/L),N_kmod)
#     theta = np.linspace(-0.5*np.pi, 0.5*np.pi, N_theta)
#     kmod_tile, theta_tile = np.meshgrid(kmod,theta)
    
#     """ JONSWAP without shape modification gamma """
#     F_kmod = P*kmod**(-3)*np.exp(-2.5*(kp/kmod)**2) + P*400*kmod**(-4)*np.exp(-2.5*(kp*10/kmod)**2) 
#     D_theta = np.cos(theta-thetam)**N_power
#     dtheta = theta[1]-theta[0]
#     """ Is this necessary? Yes! """
#     D_theta = D_theta/np.trapz(D_theta, theta)  # Normalize so the sum equals one
#     F_kmod_tile, D_theta_tile = np.meshgrid(F_kmod,D_theta) 
#     F_kmodtheta_tile = F_kmod_tile*D_theta_tile/kmod_tile # Notice!! Normalize by k
#     """ Uniform grid in kx,ky """
#     kx = np.arange(1,N_mode+1)*2*np.pi/L # based on the grid, interval can't go smaller then pi/L
#     ky = np.arange(-N_mode/2,N_mode/2+1)*2*np.pi/L
#     kx_tile, ky_tile = np.meshgrid(kx,ky)
#     kxp_tile, kyp_tile = pol2cart(kmod_tile, theta_tile)
    
#     """ Project from uniform k to uniform kx,ky """
#     F_kxky_tile = griddata((kxp_tile.ravel(), kyp_tile.ravel()), F_kmodtheta_tile.ravel(), 
#                            (kx_tile, ky_tile), method='linear', fill_value=0) 
#     return kmod, F_kmod, kx, ky, F_kxky_tile



''' Add modes together to generate initial eta field. Random phase.
    For visualizing purpose but essentially the same idea as the initialization in spectrum.h '''

def eta_random(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile):

    ''' Function to generate a random field given kx-ky spectrum and the x-y-array. Only
        works with uniformly spaced kx-ky so far.
    '''
    
    np.random.seed(0) 
    # It is only for focusing case but we don't set tb, xb (time and location of breaking)
#     tb = 40; xb = 0; yb = 0
#     phase_tile = -kx_tile*xb-ky_tile*yb+np.random.random_sample(kx_tile.shape)*2*np.pi 
    phase_tile = np.random.random_sample(kx_tile.shape)*2*np.pi 
    eta_tile = np.zeros(x_tile.shape)
    
    kmod_cart_tile, theta_cart_tile = cart2pol(kx_tile,ky_tile)
    # frequency based on kx or kmod (might need to change for capillary waves)
    omega_tile = (9.8*kmod_cart_tile)**0.5
    dkx = kx_tile[0,1]-kx_tile[0,0]; dky = ky_tile[1,0]-ky_tile[0,0]
    N_grid = x_tile.shape[0]

    # Simple summation. To-do: parallelize this
    for i1 in range(0, N_grid):
        for i2 in range(0, N_grid):
            ampl = (2*F_kxky_tile*dkx*dky)**0.5
            a = (kx_tile*x_tile[i1,i2]+ky_tile*y_tile[i1,i2])-omega_tile*t+phase_tile
            mode = ampl*(np.cos(a)) # uniform spacing in kx and ky
            eta_tile[i1,i2] = np.sum(mode)    
    return eta_tile, phase_tile









