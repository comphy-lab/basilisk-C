# Need to include MYLIBS dir to python path 
from py_functionsBA.AP import *

Stud = 'RS'
RESTORE     = 0

mu_r        = 1.
rho_r       = 1./0.9
PHI         = 0.1
Bo          = 1
Ga          = 1
N           = 1
Dim         = 3
ndc         = 10

def run_run():
    name        = f'/Dim_{Dim}/N_{N}/PHI_{PHI}/rho_r_{rho_r:.2f}/Bo_{Bo}/mu_r_{mu_r}'
    name        = f'TEST/'
    PS = AP('Ga',[100],name,name_of_C_file=Stud)
    # Main parameters 
    PS.mu_r         = mu_r
    PS.rho_r        = rho_r
    PS.Bo           = Bo
    PS.PHI          = PHI
    PS.N            = N
    PS.Ga           = Ga
    # set time scale and HowManyTime Scale
    PS.TimeScale    = 'Tg'
    PS.HMT          = 500
    PS.TMAX         = 100
    PS.dimension    = Dim
    # PS.SAVESTEP     = PS.TMAX/3
    PS.MOVIES       = 1
    PS.ndcmin       = ndc
    PS.nProc        = 4
    PS.Nbsup        = 100
    PS.RESTORE      = RESTORE
    PS.CORREC       = 1
    PS.DEBUG        = 0
    PS.COAL         = 0
    PS.run()

run_run()
