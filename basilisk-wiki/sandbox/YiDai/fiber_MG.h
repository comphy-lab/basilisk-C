/**
## Implementation of grass fiber
We follow [Huang 2007](https://doi.org/10.1016/j.jcp.2007.07.002) to imeplement flexible fibers. 
*/

#define ALPHA (-2E3)
#define BETA (-10.)
#define GAMMA (1E-2)
// #define GAMMA 0.
bool BC_CLAMP;
bool GRAVITY;
bool INTERA;
coord GCAST = {0., -10};
double DensityD = 1.; // Rho_fiber / Rho_air

typedef struct FNode
{
    coord post0;
    coord post1;
    coord post2;                    // the position of t0, t1, t2 are record as diag position
    coord Xstari;
    double Tt0;            
    double Tt1;                     // Tension force at current and previous time step
    coord Ui;                       // the velocity of the nodes
    coord Uib;                      // the interpolation velocity from fluid to fiber -- later on will be delete
    coord accU;                     // accumulative velocity difference, later on will be delete
    coord LagForce;                 // Lag force
    coord BendForce;                // Bending force
    #if _MPI
        int pid;
    #endif
    Cache stencil;                  // stencil for locating the neighbor points
} FNode;

typedef struct lagfiber
{
    int NN;                         // the number of sections
    int nlp;                        // the number of lag points
    double Lr;                      // length of the fiber
    FNode *nodes;                   // the array of lag nodes
    double D_s;                     // Delta_s = the spacing of arclength (assumption)
    bool isactive;
} lagfiber;

#ifndef NFBS
    #define NFBS 1
#endif

typedef struct lagrass
{
    lagfiber fb[NFBS];
    int nfbs;
} lagrass;
lagrass fbs;
#define FB(i) (fbs.fb[i])


/**
define the stencil size
 */
#if dimension < 3
#define STENCIL_SIZE 25
#else
#define STENCIL_SIZE 125
#endif

#define ACROSS_PERIODIC(a, b) (fabs(a - b) > L0 / 2.)
#define PERIODIC_1DIST(a, b) (fabs(a - L0 - b) > L0 / 2. ? a + L0 - b : a - L0 - b)
#define GENERAL_1DIST(a, b) (ACROSS_PERIODIC(a, b) ? PERIODIC_1DIST(a, b) : a - b)

// some MACRO for calculation
#if dimension < 3
#define cnorm(a) (sqrt(sq(a.x) + sq(a.y)))
#define cdot(a, b) (a.x * b.x + a.y * b.y)
#define Vdiff(a, b) ((coord){(a).x - (b).x, (a).y - (b).y})
#define star_Vdiff(a, b) ((coord){(a).x * 2. - (b).x, (a).y * 2. - (b).y})
#define sdot(a) (a.x * a.x + a.y * a.y)
#else
#define cnorm(a) (sqrt(sq(a.x) + sq(a.y) + sq(a.z)))
#define cdot(a, b) (a.x * b.x + a.y * b.y + a.z * b.z)
#define Vdiff(a, b) ((coord){(a).x - (b).x, (a).y - (b).y, (a).z - (b).z})
#define sdot(a) (a.x * a.x + a.y * a.y + a.z * a.z)
#endif

void initialize_fiber(lagfiber * fiber)
{
    // filling the details of fiber here
    double Lr = 1.;                              // length of the fiber
    int NN = 99;              // section on the fiber (N)
    int nlp = 100;             // nodes on the fiber (N+1)
    double D_s = Lr / NN;                 // the spacing of the node
    fiber->NN = NN;
    fiber->nlp = nlp;
    fiber->Lr = Lr;
    fiber->D_s = D_s;
    fiber->nodes = malloc(nlp * sizeof(FNode));

    double k = 0.1 * pi;
    for (int i = 0; i < nlp; i++)
    {
        fiber->nodes[i].post0.x = X0 + L0 / 2. - sin(k) + i * D_s * sin(k);
        fiber->nodes[i].post0.y = Y0 + 4. - cos(k) + i * D_s * cos(k);
        fiber->nodes[i].post1.x = X0 + L0 / 2. - sin(k) + i * D_s * sin(k);
        fiber->nodes[i].post1.y = Y0 + 4. - cos(k) + i * D_s * cos(k);
        fiber->nodes[i].Tt1 = 1;
        fiber->nodes[i].Xstari = star_Vdiff(fbs.fb[0].nodes[i].post1, fbs.fb[0].nodes[i].post0);
        foreach_dimension()
        {
            fbs.fb[0].nodes[i].post2.x = 2 * fbs.fb[0].nodes[i].post1.x - fbs.fb[0].nodes[i].post0.x;
            fiber->nodes[i].Ui.x = 0.;
            fiber->nodes[i].accU.x = 0.;
            fiber->nodes[i].LagForce.x = 0;
            fiber->nodes[i].BendForce.x = 0;
        }
        fiber->nodes[i].stencil.n = STENCIL_SIZE;
        fiber->nodes[i].stencil.nm = STENCIL_SIZE;
        fiber->nodes[i].stencil.p = malloc(STENCIL_SIZE * sizeof(Index));
    }
    BC_CLAMP = false;
    GRAVITY = true;
    INTERA = false;
}

event init(i = 0)
{
    initialize_fiber(&FB(0));
}

void free_fiber(lagfiber *fiber)
{
    for (int i = 0; i < fiber->nlp; i++)
        free(fiber->nodes[i].stencil.p);
    free(fiber->nodes);
}

void free_fibers(lagrass *gfbs)
{
    for (int i = 0; i < gfbs->nfbs; i++)
        if (fbs.fb[i].isactive)
                free_fiber(&(gfbs->fb[i]));
}

/**
 The bending force at X*
*/
void bendingforce(lagfiber *fiber)
{
    int nlp = fiber->nlp;
    double D_s = fiber->D_s;
    coord Xstari[nlp];
    for (int k = 0; k < nlp; k++)
        Xstari[k] = fiber->nodes[k].Xstari;

    // check Fb(X*) --- BC at free end (i = 0): eq24; at fixed end (i = N): approximate as 0;
    // the exception inner BC i = 1 and i = nlp - 2
    coord bc_clp = {0., 1.};
    double Lgamma = -GAMMA * pow(D_s, -4);
    coord DDXstari[nlp];
    for (int i = 1; i < nlp - 1; i++)
        foreach_dimension()
            DDXstari[i].x = Xstari[i + 1].x - 2. * Xstari[i].x + Xstari[i - 1].x;

    foreach_dimension()
    {
        fiber->nodes[0].BendForce.x = Lgamma * (DDXstari[2].x - DDXstari[1].x);
        fiber->nodes[1].BendForce.x = Lgamma * (DDXstari[2].x - 2. * DDXstari[1].x);
        if (BC_CLAMP == true)
        {
            fiber->nodes[nlp - 2].BendForce.x = Lgamma * (2. * bc_clp.x * D_s - 2. * (Xstari[nlp - 1].x - Xstari[nlp - 2].x) -
                                                          2. * DDXstari[nlp - 2].x + DDXstari[nlp - 3].x);
        }
        else
        {
            fiber->nodes[nlp - 2].BendForce.x = Lgamma * (-2. * DDXstari[nlp - 2].x + DDXstari[nlp - 3].x);
        }
    }
    for (int i = 2; i < nlp - 2; i++)
    {
        foreach_dimension()
        {
            fiber->nodes[i].BendForce.x = Lgamma * (DDXstari[i + 1].x - 2. * DDXstari[i].x + DDXstari[i - 1].x);
        }
    }
}

/**
we solve the poisson equation in the form of dT^{2}_{x+1/2}/ds^{2} = f(s)
use multi-grid method
 */

void relaxT_GS(double *X, double *f, coord *diffXstari, int l_NN)
{
    double Lr = fbs.fb[0].Lr;
    int l_nlp = l_NN + 1;
    double l_D_s = Lr / l_NN;
    for (int k = 1; k < l_nlp - 2; k++)
    {
        X[k] = 0.5 / sdot(diffXstari[k + 0]) *
                (X[k + 1] * cdot(diffXstari[k + 1], diffXstari[k + 0]) +
                X[k - 1] * cdot(diffXstari[k + 0], diffXstari[k - 1]) -
                pow(l_D_s, 4) * f[k]);
    }
    X[0] = 1. / 3. / sdot(diffXstari[0]) *
            (X[1] * cdot(diffXstari[1], diffXstari[0]) - pow(l_D_s, 4) * f[0]);
    X[l_nlp - 2] = 1. / sdot(diffXstari[l_nlp - 2]) *
                    (X[l_nlp - 3] * cdot(diffXstari[l_nlp - 2], diffXstari[l_nlp - 3]) -
                    pow(l_D_s, 4) * f[l_nlp - 2]);
}

double resid(double *X, double *f, coord *diffXstari, double *residual, int l_NN)
{
    double resb = -1.;
    double Lr = fbs.fb[0].Lr;
    int l_nlp = l_NN + 1;
    double l_D_s = Lr / l_NN;
    // calculate res =  dT^{2}_{x+1/2}/ds^{2} - f(s)
    residual[0] = 3. * X[0] * sdot(diffXstari[0]) -
               (X[1] * cdot(diffXstari[1], diffXstari[0]) -
                pow(l_D_s, 4) * f[0]);
    residual[l_nlp - 2] = X[l_nlp - 2] * sdot(diffXstari[l_nlp - 2]) -
               (X[l_nlp - 3] * cdot(diffXstari[l_nlp - 2], diffXstari[l_nlp - 3]) -
                pow(l_D_s, 4) * f[l_nlp - 2]);
    for (int k = 1; k < l_nlp - 2; k++)
    {
        residual[k] = 2. * X[k] * sdot(diffXstari[k + 0]) -
                   (X[k + 1] * cdot(diffXstari[k + 1], diffXstari[k + 0]) +
                    X[k - 1] * cdot(diffXstari[k + 0], diffXstari[k - 1]) -
                    pow(l_D_s, 4) * f[k]);
    }
    for (int ki = 0; ki < l_nlp - 1; ki++){
        if (fabs(residual[ki]) > resb)
            resb = fabs(residual[ki]);
    }
    return resb;
}

void V_cycle(double *X, double *f, coord *diffXstari, int l_NN, int num_iter){
    int iter = 0;
    // int num_iter = 10;
    do
    { // GS
        relaxT_GS(X, f, diffXstari, l_NN);
        iter++;
    } while (iter < num_iter);
    double residual[l_NN];
    resid(X, f, diffXstari, residual, l_NN);
    // restriction
    int NNc = l_NN / 2 + 1;
    double fC[NNc];
    coord DXC[NNc];
    double XC[NNc];
    for (int ri = 0; ri < NNc; ri++)
    {
        fC[ri] = residual[ri * 2];
        XC[ri] = 0.;
        foreach_dimension() DXC[ri].x = diffXstari[ri * 2].x;
    }
    int NN = fbs.fb[0].NN;
    if (NNc < NN * 3 / 4){
        do
        { // GS
            relaxT_GS(XC, fC, DXC, NNc);
            iter++;
        } while (iter < num_iter);
    }
    else{
        V_cycle(XC, fC, DXC, NNc, num_iter);
    }
    double e[l_NN];
    for (int ie = 0; ie < l_NN; ie += 2){
        int ind = ie / 2;
        e[ie] = XC[ind];
        X[ie] = X[ie] - e[ie];
    }
    for (int io = 1; io < l_NN; io += 2)
    {
        int ind = io / 2;
        e[io] = (XC[ind] + XC[ind + 1]) / 2.;
        X[io] = X[io] - e[io];
    }
}

event Tensionpoisson(i++)
{
    dt = dtnext(DT);
    int nlp = fbs.fb[0].nlp;
    double D_s = fbs.fb[0].D_s;
    coord diffXstari[nlp-1];
    for (int k = 0; k < nlp-1; k++)
    {
        diffXstari[k] = Vdiff(fbs.fb[0].nodes[k+1].Xstari, fbs.fb[0].nodes[k].Xstari);
    }
    // calculate the LagForce and bending force
    
    // interactionforce(&FB(0));
    bendingforce(&FB(0));
    
    // we provide RHS1,2 with i = [0, N-1]; 
    double RHS1[nlp - 1];
    double RHS2[nlp - 1];
    double RHS3[nlp - 1];
    double f[nlp - 1];
    for (int i = 0; i < nlp - 1; i++)
    {
        double dotXi1_p1 = sdot(Vdiff(fbs.fb[0].nodes[i + 1].post1, fbs.fb[0].nodes[i].post1));
        double dotXi1_p0 = sdot(Vdiff(fbs.fb[0].nodes[i + 1].post0, fbs.fb[0].nodes[i].post0));
        RHS1[i] = 0.5 / sq(dt) * (1 - 2. * dotXi1_p1 / sq(D_s) + dotXi1_p0 / sq(D_s));
        RHS2[i] = -1. / sq(D_s) * sdot(Vdiff(fbs.fb[0].nodes[i + 1].Ui, fbs.fb[0].nodes[i].Ui));
    }
    // RHS3 with i = [1, N-2]
    for (int j = 0; j < nlp - 2; j++)
    {
        coord diffBendi1 = Vdiff(fbs.fb[0].nodes[j + 1].BendForce, fbs.fb[0].nodes[j].BendForce);
        coord diffLagi1  = Vdiff(fbs.fb[0].nodes[j + 1].LagForce, fbs.fb[0].nodes[j].LagForce);
        RHS3[j] = -1. / sq(D_s) * cdot(diffXstari[j], Vdiff(diffBendi1, diffLagi1));
        f[j] = RHS1[j] + RHS2[j] + RHS3[j];
    }
    coord diffBendLagN_1 = Vdiff(fbs.fb[0].nodes[nlp - 2].LagForce, fbs.fb[0].nodes[nlp - 2].BendForce);
    coord Gra_vector = {GRAVITY * GCAST.x * DensityD, GRAVITY * GCAST.y * DensityD};
    f[nlp - 2] = RHS1[nlp - 2] + RHS2[nlp - 2] + cdot(diffXstari[nlp - 2], Vdiff(diffBendLagN_1, Gra_vector)) / sq(D_s);

    double Tt1[nlp - 1];
    for (int jj = 0; jj < nlp - 1; jj++)
    {
        Tt1[jj] = fbs.fb[0].nodes[jj].Tt1;
    }
    int NN = fbs.fb[0].NN;
    double residual[NN];
    int Mnum_cyc = 300;
    int num_iter = 5;
    if (t == 0)
        num_iter = 10;
    double tolerance = 1e-12;
    double resb, resa;
    resb = resa = resid(Tt1, f, diffXstari, residual, NN);
    for (int ic = 0; 
         ic < Mnum_cyc && (ic < 1 || resa > tolerance);
         ic++)
    {
        V_cycle(Tt1, f, diffXstari, NN, num_iter);
        resa = resid(Tt1, f, diffXstari, residual, NN);
        if (resa > tolerance)
        {
            if (resb / resa < 2 && num_iter < 300)
                num_iter++;
            else if (resb / resa > 10 && num_iter > 2)
                num_iter--;
        }
        resb = resa;
        // if (!((int)(t / 0.01) % 2))
        //     printf("%g %d %g\n", t, ic, resa);
    }
    for (int kk = 0; kk < nlp - 1; kk++)
    {
        fbs.fb[0].nodes[kk].Tt1 = Tt1[kk];
    }
}

event advectfiber(i++, last)
{
    int nlp = fbs.fb[0].nlp;
    double D_s = fbs.fb[0].D_s;
    coord Xstari[nlp];
    for (int k = 0; k < nlp; k++)
        Xstari[k] = fbs.fb[0].nodes[k].Xstari;

    double Tt1[nlp - 1];
    for (int jj = 0; jj < nlp - 1; jj++)
    {
        Tt1[jj] = fbs.fb[0].nodes[jj].Tt1;
    }

    coord res_X[nlp];
    for (int i = 0; i < nlp - 1; i++)
    {
        foreach_dimension()
        {
            res_X[i].x = fbs.fb[0].nodes[i].BendForce.x; // because of the second derivative, the bending force switch sign
            res_X[i].x += GRAVITY * GCAST.x;
            res_X[i].x += -fbs.fb[0].nodes[i].LagForce.x;
        }
    }

    foreach_dimension()
        fbs.fb[0].nodes[0].post2.x = 1. / (1. / sq(dt) + 2 * Tt1[0] / sq(D_s)) *
                       (2 * Tt1[0] / sq(D_s) * fbs.fb[0].nodes[1].post2.x +
                        Xstari[0].x / sq(dt) + res_X[0].x);
    for (int j = 1; j < nlp - 1; j++)
    {
        foreach_dimension()
        {
            fbs.fb[0].nodes[j].post2.x = 1. / (sq(D_s / dt) + Tt1[j] + Tt1[j - 1]) *
                                         (Tt1[j] * fbs.fb[0].nodes[j + 1].post2.x +
                                          Tt1[j - 1] * fbs.fb[0].nodes[j - 1].post2.x +
                                          sq(D_s / dt) * Xstari[j].x + sq(D_s) * res_X[j].x);
        }
    }
}

// update the fiber
event update_fiber(i++, last){
    // set the field for next time step
    int nlp = fbs.fb[0].nlp;
    for (int jj = 0; jj < nlp; jj++)
    {
        foreach_dimension()
        {
            fbs.fb[0].nodes[jj].post0.x = fbs.fb[0].nodes[jj].post1.x;
            fbs.fb[0].nodes[jj].post1.x = fbs.fb[0].nodes[jj].post2.x;
            fbs.fb[0].nodes[jj].post2.x = 2 * fbs.fb[0].nodes[jj].post1.x - fbs.fb[0].nodes[jj].post0.x;
        }
    }
    for (int k = 0; k < nlp; k++)
    {
        fbs.fb[0].nodes[k].Xstari = star_Vdiff(fbs.fb[0].nodes[k].post1, fbs.fb[0].nodes[k].post0);
    }
}

event logdata(i++){
    char loc_y[60];
    snprintf(loc_y, 60, "loc_y");
    static FILE *fp_loc_y = fopen(loc_y, "w");
    fprintf(fp_loc_y, "%g\t%g\n", t, fbs.fb[0].nodes[0].post1.x);
    fflush(fp_loc_y);
}

event logerror(i++)
{
    double D_s = fbs.fb[0].D_s;
    int nlp = fbs.fb[0].nlp;
    double er = -1.;
    for (int ki = 0; ki < nlp-1; ki++)
    {
        double sig = fabs(sq((fbs.fb[0].nodes[ki + 1].post1.x - fbs.fb[0].nodes[ki].post1.x) / D_s) + sq((fbs.fb[0].nodes[ki + 1].post1.y - fbs.fb[0].nodes[ki].post1.y) / D_s) - 1.);
        if (sig > er)
            er = sig;
    }
    char sig_name[60];
    snprintf(sig_name, 60, "sig_er");
    static FILE *fp_sig_name = fopen(sig_name, "w");
    fprintf(fp_sig_name, "%g\t%g\n", t, er);
    fflush(fp_sig_name);
}

/**
## Reference
Huang, W.X., Shin, S.J. and Sung, H.J., 2007. Simulation of flexible filaments in a uniform flow by the immersed boundary method. Journal of computational physics, 226(2), pp.2206-2228.
*/