/**
# Test Quantities for Batchelor Vortex

This test computes and saves various quantities related to the Batchelor vortex
correction, including base flow and first-order perturbations.

~~~pythonplot correction
import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('results_quantities_Uc0.0.txt', delimiter='\t', usecols=[0,1,2])
data2 = np.loadtxt('results_quantities_Uc0.1.txt', delimiter='\t', usecols=[0,1,2])
data3 = np.loadtxt('results_quantities_Uc1.0.txt', delimiter='\t', usecols=[0,1,2])

fig, ax1 = plt.subplots(ncols=1, sharex=False, figsize=(3.5,3.5))

ax1.plot(data1[:,0], data1[:,1]/data1[:,0], label='$U_c=0.0$', ls='-');
ax1.plot(data2[:,0], data2[:,1]/data2[:,0], label='$U_c=0.1$', ls='--');
ax1.plot(data3[:,0], data3[:,1]/data3[:,0], label='$U_c=1.0$', ls='-.'); 
ax1.legend()
ax1.set_xlabel('$r$')
ax1.set_ylabel('$\\hat{\\psi}^{(1)}/r$')
ax1.set_xlim([0,4])
ax1.set_ylim([0,1.4])
plt.tight_layout()
plt.savefig('test_correction.svg')
~~~

~~~pythonplot correction
import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('results_quantities_Uc0.0.txt', delimiter='\t', usecols=[0,4,5,6])
data2 = np.loadtxt('results_quantities_Uc0.1.txt', delimiter='\t', usecols=[0,4,5,6])
data3 = np.loadtxt('results_quantities_Uc1.0.txt', delimiter='\t', usecols=[0,4,5,6])

fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharex=False, figsize=(3.5*3,3.5))

labels = [r'$u_r^{(0)}$', r'$u_\theta^{(0)}$', r'$u_\xi^{(0)}$']

for i, ax in enumerate((ax1, ax2, ax3)):
  ax.plot(data1[:,0], data1[:,i+1], label='$U_c=0.0$', ls='-');
  ax.plot(data2[:,0], data2[:,i+1], label='$U_c=0.1$', ls='--');
  ax.plot(data3[:,0], data3[:,i+1], label='$U_c=1.0$', ls='-.'); 
  ax.legend()
  ax.set_xlabel('$r$')
  ax.set_ylabel(labels[i])
  ax.set_xlim([0,4])
plt.tight_layout()
plt.savefig('test_correction2.svg')
~~~

~~~pythonplot correction
import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('results_quantities_Uc0.0.txt', delimiter='\t', usecols=[0,7,8,9])
data2 = np.loadtxt('results_quantities_Uc0.1.txt', delimiter='\t', usecols=[0,7,8,9])
data3 = np.loadtxt('results_quantities_Uc1.0.txt', delimiter='\t', usecols=[0,7,8,9])

fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharex=False, figsize=(3.5*3,3.5))

labels = [r'$\omega_r^{(0)}$', r'$\omega_\theta^{(0)}$', r'$\omega_\xi^{(0)}$']

for i, ax in enumerate((ax1, ax2, ax3)):
  ax.plot(data1[:,0], data1[:,i+1], label='$U_c=0.0$', ls='-');
  ax.plot(data2[:,0], data2[:,i+1], label='$U_c=0.1$', ls='--');
  ax.plot(data3[:,0], data3[:,i+1], label='$U_c=1.0$', ls='-.'); 
  ax.legend()
  ax.set_xlabel('$r$')
  ax.set_ylabel(labels[i])
  ax.set_xlim([0,4])
plt.tight_layout()
plt.savefig('test_correction3.svg')
~~~

~~~pythonplot correction
import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('results_quantities_Uc0.0.txt', delimiter='\t', usecols=[0,10,11,12])
data2 = np.loadtxt('results_quantities_Uc0.1.txt', delimiter='\t', usecols=[0,10,11,12])
data3 = np.loadtxt('results_quantities_Uc1.0.txt', delimiter='\t', usecols=[0,10,11,12])

fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharex=False, figsize=(3.5*3,3.5))

labels = [r'$u_r^{(1)}$', r'$u_\theta^{(1)}$', r'$u_\xi^{(1)}$']

for i, ax in enumerate((ax1, ax2, ax3)):
  ax.plot(data1[:,0], data1[:,i+1], label='$U_c=0.0$', ls='-');
  ax.plot(data2[:,0], data2[:,i+1], label='$U_c=0.1$', ls='--');
  ax.plot(data3[:,0], data3[:,i+1], label='$U_c=1.0$', ls='-.'); 
  ax.legend()
  ax.set_xlabel('$r$')
  ax.set_ylabel(labels[i])
  ax.set_xlim([0,4])
plt.tight_layout()
plt.savefig('test_correction4.svg')
~~~

~~~pythonplot correction
import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('results_quantities_Uc0.0.txt', delimiter='\t', usecols=[0,13,14,15])
data2 = np.loadtxt('results_quantities_Uc0.1.txt', delimiter='\t', usecols=[0,13,14,15])
data3 = np.loadtxt('results_quantities_Uc1.0.txt', delimiter='\t', usecols=[0,13,14,15])

fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharex=False, figsize=(3.5*3,3.5))

labels = [r'$\omega_r^{(1)}$', r'$\omega_\theta^{(1)}$', r'$\omega_\xi^{(1)}$']

for i, ax in enumerate((ax1, ax2, ax3)):
  ax.plot(data1[:,0], data1[:,i+1], label='$U_c=0.0$', ls='-');
  ax.plot(data2[:,0], data2[:,i+1], label='$U_c=0.1$', ls='--');
  ax.plot(data3[:,0], data3[:,i+1], label='$U_c=1.0$', ls='-.'); 
  ax.legend()
  ax.set_xlabel('$r$')
  ax.set_ylabel(labels[i])
  ax.set_xlim([0,4])
plt.tight_layout()
plt.savefig('test_correction5.svg')
~~~

*/

#include "grid/octree.h"
#include "view.h"
#include "acastillo/input_fields/filaments.h"

void compute_and_save(double U_c, const char* filename) {
  double rho;
  double a = 1.0;
  double kappa = 0.05; 

  if (pid() == 0){
    FILE *file = fopen(filename, "w");
    
    for (rho = 0.0; rho <= 4.0; rho += 0.05) {
      integration_results results_A;
      // Compute A, A', A''
      compute_A_with_derivatives(rho, a, U_c, &results_A);

      double A = results_A.A;
      double dA = results_A.A_p;
      double d2A = results_A.A_pp;
      
      // Variables for quantities
      double u0_r, u0_th, u0_xi;
      double w0_r, w0_th, w0_xi;
      double u1_r, u1_th, u1_xi;
      double w1_r, w1_th, w1_xi;


      double S0 = 0.0;
      double denom = 0.0;
      double dS0_dr = 0.0;
      if (rho < 1e-9) {
        // Handle singularity at rho = 0
        // At rho=0, u0_th -> 0, w0_th -> 0
        u0_r = 0.0; u0_th = 0.0; u0_xi = U_c;
        w0_r = 0.0; w0_th = 0.0; w0_xi = 2.0/(a*a);
        
        // Perturbations might be singular or zero depending on A behavior near 0
        // A ~ rho^2 near 0? If so A/rho^2 is finite.
        // For now set to 0 to avoid NaN if not handled by limits
        u1_r = 0.0; u1_th = 0.0; u1_xi = 0.0;
        w1_r = 0.0; w1_th = 0.0; w1_xi = 0.0;
      } else {
        double g_rho = exp(-sq(rho / a));
        double rho2 = sq(rho);
        double a2 = sq(a);

        // --- Zeroth-Order Base Flow (u0, w0) ---
        u0_r = 0.0;
        u0_th = (1.0 / rho) * (1.0 - g_rho);
        u0_xi = U_c * g_rho;

        w0_r = 0.0;
        w0_th = U_c * (2.0 * rho / a2) * g_rho;
        w0_xi = (2.0 / a2) * g_rho;

        // --- Intermediate terms needed for first-order terms ---
        // --- Swirl number at leading order
        S0 = rho * w0_th / u0_th; 
        denom = 1.0 - g_rho;
        dS0_dr = (3.0 * rho2 - 2.0 * sq(rho2) / a2) * g_rho / denom
                      - (2.0 * sq(rho2) / a2) * g_rho / (denom * denom);
        dS0_dr *= 2.0 * U_c / a2;

        u1_r  = -A/rho2;
        u1_th = u0_th - dA/rho;
        u1_xi = u0_xi + S0*A/rho2;

        w1_r  = -S0*A/cube(rho);
        w1_th = w0_th + S0*A/cube(rho) - S0*dA/rho2 - dS0_dr*A/rho2;
        w1_xi = u0_th/rho + w0_xi - d2A/rho - dA/rho2 + A/cube(rho);

        u1_r  *= rho * kappa;
        u1_th *= rho * kappa; 
        u1_xi *= rho * kappa;

        w1_r  *= rho * kappa;
        w1_th *= rho * kappa;
        w1_xi *= rho * kappa;
      }

      fprintf(file, "%.4f\t%.4e\t%.4e\t%.4e\t", rho,A,dA,d2A);
      fprintf(file, "%.4e\t%.4e\t%.4e\t", u0_r, u0_th, u0_xi);
      fprintf(file, "%.4e\t%.4e\t%.4e\t", w0_r, w0_th, w0_xi);
      fprintf(file, "%.4e\t%.4e\t%.4e\t", u1_r, u1_th, u1_xi);
      fprintf(file, "%.4e\t%.4e\t%.4e\t", w1_r, w1_th, w1_xi);
      fprintf(file, "%.4e\t%.4e\t%.4e\t", S0, denom, dS0_dr);
      fprintf(file, "\n");
    }
    fclose(file);
  }
}

int main() {
  compute_and_save(0.0, "results_quantities_Uc0.0.txt");
  compute_and_save(0.1, "results_quantities_Uc0.1.txt");
  compute_and_save(1.0, "results_quantities_Uc1.0.txt");
  return 0;
}
