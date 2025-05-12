/** Riemann solver for St. Venant.  Uses wave propagation algorithm.

Compute waves, speeds and fluctuations at each interface for the 
   scalar advection problem. 
*/


double grav = 1;

void rp1_swe(int meqn, int mwaves, double *ql, double *qr, double *waves, 
                double *speeds, double *amdq, double *apdq)
{
    if (meqn != 2) 
    {
        printf("Error : meqn != 2\n");
        exit(0);
    }
    double  hl = ql[0],  hr = qr[0];
    double hul = ql[1], hur = qr[1];

    double delta[2];  
    for(int m = 0; m < meqn; m++)
    {
        delta[m] = qr[m] - ql[m];
    }

    if (hl <= 0 || hr <= 0) 
    {
        printf("hr, hl <= 0;   %g %g\n", hl, hr);
        exit(0);
    }

    double ubar, cbar, divsqrt;
    divsqrt = sqrt(hl) + sqrt(hr);

    ubar = (hul/sqrt(hl) + hur/sqrt(hr))/divsqrt;
    cbar = sqrt(0.5*grav*(hr + hl));

    double a1, a2;
    a1 = 0.5*(-delta[1] + (ubar + cbar) * delta[0])/cbar;
    a2 = 0.5*( delta[1] - (ubar - cbar) * delta[0])/cbar;

    int mw = 0;
    waves[0 + mw*meqn] = a1;
    waves[1 + mw*meqn] = a1*(ubar - cbar);
    speeds[mw] = ubar - cbar;

    mw = 1;
    waves[0 + mw*meqn] = a2;
    waves[1 + mw*meqn] = a2*(ubar + cbar);
    speeds[mw] = ubar + cbar;

    for(int m = 0; m < meqn; m++)
    {
        amdq[m] = 0;
        apdq[m] = 0;
    }

    for(int mw = 0; mw < mwaves; mw++)
        for(int m = 0; m < meqn; m++)
            if (speeds[mw] < 0) 
                amdq[m] += speeds[mw]*waves[m+mw*meqn];
            else
                apdq[m] += speeds[mw]*waves[m+mw*meqn];               
}
