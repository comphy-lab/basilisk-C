import numpy as np
import sys
import math 
import os
import filecmp
import shutil
class AP():
    def __init__(self,para_name: str,parameters_list: list,namestud: str,name_of_C_file='RS'):
        self.para_name = para_name
        self.parameters_list = parameters_list
        self.resultsdir = "../results/"
        self.namestud=namestud
        self.namedir 	= self.resultsdir+namestud.strip()+'/'
        self.dumpfile = 'last-dump'
        self.name_of_C_file = name_of_C_file
        self.MPI = True
        self.nProc = 9
        ## Defaults adimensional parameters 
        # physical parameters
        self.Bo  = 0.5              # Bond number
        self.Ga = 55                # Galileo number
        self.rho_r = 1.1668599      # ratio of rho
        self.mu_r = 0.4237805       # ratio between mu
        self.PHI = 0.15             # concentration phi
        self.N = 10                 # sqrt of the number of bubbles 
        self.Nbsup = 2              # number of sup bubble that need to be in case there is breakage
        # numerical parameters 
        self.ndcmin         = 40
        self.LDn            = self.ndcmin        # number of layers per diameter of drops for continuity waves computation
        self.TMAX           = 2000         # end time 
        self.dtprint        = 0.1          # writing time step
        self.FIRSTSTEP      = self.dtprint * 1
        self.SAVESTEP       = self.TMAX/3
        self.MOVIES         = 1
        # geometrical parameters 
        self.dimension = 2                # set the dimension of the simulation 
        self.nopara    = ''
        
        ## constant dimensional parameters 
        self.g = 1.                 # gravity 
        self.r1 = 1.
        self.D = 1.                 # diameter of bubbles
        self.rho_d = 1.             # density of the dispersed phase 
        
        ## depending parameters 
        # Physiical parameters
        if self.dimension == 2:                                                   # if it is in 2d
            self.Nb = self.N**2                              # number of cells 
            self.Ls = self.N*self.D*math.sqrt(math.pi/(4.*self.PHI))        # size of the domain
            self.Vol = math.pi * self.D**2/4.
        if self.dimension == 3:                                                   # if it is in 3d
            self.Nb = self.N**3                         # number of cells 
            self.Ls = self.N*self.D*math.pow(math.pi/(6.*self.PHI),1./3.)   # size of the domain
            self.Vol = math.pi * (self.D/2.)**3 * 4./3.
        self.mu_d = math.sqrt(abs(1-self.rho_r)*self.rho_r*self.g*self.D**3)*self.rho_d/(self.Ga*self.mu_r);    
        self.mu_f = self.mu_d*self.mu_r
        self.rho_f = self.rho_d*self.rho_r
        self.sig = abs(1-self.rho_r)*self.rho_d*self.g*self.D**2/self.Bo
        # Geometrical parameters
        self.DP = self.Ls /self.N           # size of a unit cell
        self.RND = 0.49*(self.DP-self.D)    # maximum displacement of the buuble from the center of the unit cell
        # numerical parameters
        self.LEVEL =int( math.ceil(np.log(self.ndcmin*self.Ls/self.D)/np.log(2.))  )             # number of cells per diameters of a bubble
        self.ndc = 2**self.LEVEL/self.Ls*self.D
        self.ESP2 = 1e-6
        self.Dx = self.Ls / 2**self.LEVEL
        # sedimentational velocity
        self.Ug = 2./9. * self.D**2/4. / self.mu_f *(self.rho_d - self.rho_f) *self.g
        # others non dimensional numbers 
        self.Oh = self.mu_d / math.sqrt(self.rho_d*self.D/2.*self.sig)
        self.St = self.rho_d * self.D * self.Ug/self.mu_f/18.
        self.Re = self.rho_f * self.D * self.Ug/self.mu_f
        # active artificial coalescence 
        self.COAL = 0
        self.RAND = 0
        

        
        # restore
        self.CORREC = 1 
        self.DEBUG = 0
        self.RESTORE = 0 
        self.RESTORE_DIR  = self.name_of_C_file
        
        # period of the smallest capilarity wave
        self.h = self.Ls / 2**self.LEVEL
        self.rho_m = (self.rho_d + self.rho_f)/2.
        self.Tsig = math.sqrt(self.rho_m*self.D**3/(math.pi * self.sig))
        self.TimeScale =0
        self.HMT =1
        
        # Acceleration scale 
        self.Asig = self.sig / (self.rho_m * pow(self.D,self.dimension))
        self.Fsig = self.sig * self.D
        
        # Viscus and inertial time 
        self.Drho = abs(self.rho_d-self.rho_f)
        
        self.Ug = np.sqrt(self.g * self.D * self.Drho /self.rho_f)
        self.Tg = self.D / self.Ug 
        
        
        self.Uv = self.Drho * self.g /self.mu_f * self.D**2         
        self.Tv = self.D / self.Uv
        
        self.Used = self.D**2 / self.mu_f * abs(self.rho_d - self.rho_f) *self.g
        self.Tsed = self.D / self.Used
        
        self.F = self.Drho * self.g * self.Vol
        self.Mhg = self.D * self.F 
        # first update of the parameters 
        self.parameters = dict()
        self.update_parameters_dictionnary()
        
    def run(self):
        # initialize dump dir
        os.system('rm -rf dump')
        if(self.RESTORE): 
            os.system('cp -r '+self.RESTORE_DIR+' dump')
            try:
                self.parameters['LOI'] = int(open("dump/number_of_vof.txt","r").read())
            except: 
                print("COULD NOT RESTORE !")

        os.system('sed -i "s/#define dimension.*/#define dimension '+str(self.dimension)+'/" '+str(self.name_of_C_file)+'.c')
        os.system('make clean')
        os.system('rm -rf .qcc*') 
        os.system('mkdir -p '+self.resultsdir)
        os.system('mkdir -p '+self.namedir) 
        print('\n')
        print('Study on ',self.para_name,self.parameters_list)
        print('in ',self.namedir)
        print('\n')
        for P in self.parameters_list:
            os.system('rm -rf '+self.name_of_C_file)         
            self.set_parameters_dictionnary(P)
            self.print_parameters()
            if self.nProc > 1:
                print('RUNINNG parallel')
                cmd = 'CC="mpicc -D_MPI='+str(int(self.nProc)) +'"'+ ' make '+self.name_of_C_file+'.tst'
                print(cmd)
                os.system(cmd)
            else:    
                print('RUNINNG simple')
                os.system('make '+self.name_of_C_file+'.tst')
            self.print_parameters_after_sim()
            self.cp_dir_and_files(P)
    
    def set_parameters_dictionnary(self,parameter_value):
        self.update_parameters_dictionnary()
        self.parameters[self.para_name] = parameter_value
        self.calcul_parameters()
        self.update_parameters_dictionnary()
        self.parameters[self.para_name] = parameter_value
        
        if self.TimeScale:
            self.TMAX = self.parameters[self.TimeScale] * self.HMT
            self.SAVESTEP = self.TMAX/3.
        self.update_parameters_dictionnary()
        self.parameters[self.para_name] = parameter_value
    
    def calcul_parameters(self):
        if self.parameters['dimension'] == 2:                                                 
            self.Nb = self.parameters['N']**2
            Vd = math.pi*self.parameters['D']**2./4. 
            self.Ls = self.parameters['N']*math.sqrt(Vd/self.parameters['PHI'])
            self.Vol = math.pi * (self.parameters['D']/2.)**2
        if self.parameters['dimension'] == 3:                                                  
            self.Nb = self.parameters['N']**3               
            Vd = 4./3.*math.pi*(self.parameters['D']/2.)**3. 
            self.Ls = self.parameters['N']*(Vd/self.parameters['PHI'])**(1./3.)
            self.Vol = math.pi * (self.parameters['D']/2.)**3 * 4./3.
        self.mu_d = math.sqrt(abs(1-self.parameters['rho_r'])*self.parameters['rho_r']*self.parameters['g']*self.parameters['D']**3)*self.parameters['rho_d']/(self.parameters['Ga']*self.parameters['mu_r']);    
        self.mu_f = self.mu_d*self.parameters['mu_r']
        self.rho_f = self.parameters['rho_d']*self.parameters['rho_r']
        self.rho_d  = self.parameters['rho_d']
        self.D = self.parameters['D']
        self.g = self.parameters['g']
        self.sig = abs(1-self.parameters['rho_r'])*self.parameters['rho_d']*self.parameters['g']*self.parameters['D']**2/self.parameters['Bo']
        # parameters['G']eometrical parameters
        self.DP = self.Ls /self.parameters['N']           # size of a unit cell
        self.RND = 0.4*(self.DP-self.parameters['D'])    # maximum displacement of the buuble from the center of the unit cell
        # numerical parameters
        self.LEVEL =int( math.ceil(np.log(self.parameters['ndcmin']*self.Ls/self.parameters['D'])/np.log(2.)))             # number of cells per diameters of a bubble
        self.ndc = 2**self.LEVEL/self.Ls*self.parameters['D']
        # Viscous time and inertial time 

        # others non dimensional numbers 
        self.Oh = self.mu_d / math.sqrt(self.parameters['rho_d']*self.parameters['D']/2.*self.sig)
        self.St = self.parameters['rho_d'] * self.parameters['D'] * self.parameters['Ug']/self.mu_f/18.
        self.Re = self.rho_f * self.parameters['D'] * self.parameters['Ug']/self.mu_f
        # period of the smallest capilarity wave
        self.h = self.Ls / 2**self.LEVEL
        self.rho_m = (self.parameters['rho_d'] + self.rho_f)/2.
        self.Tsig = math.sqrt(self.rho_m * self.parameters['D']**3/(math.pi * self.sig))
        # Set the end time of simulation 
        self.Asig = self.sig / (self.rho_m * pow(self.parameters['D'] , self.dimension))
        self.Fsig = self.sig * self.parameters['D']
        # Compute the inertial time scale 
        self.Used = self.D**2 / self.mu_f *abs(self.rho_d - self.rho_f) *self.g
        self.Tsed = self.D / self.Used
        
        self.Drho = abs(self.rho_d-self.rho_f)
        self.Ug = np.sqrt(self.g * self.D * self.Drho /self.rho_f)

        self.Tg = self.D / self.Ug 
        self.F = self.Drho * self.g * self.Vol
        self.Mhg = self.D * self.F 
        
        self.Uv = self.Drho * self.g /self.mu_f * self.D**2         
        self.Tv = self.D / self.Uv
        self.LDn    = self.parameters['ndcmin']        # number of layers per diameter of drops
        self.Dx = self.Ls / 2**self.LEVEL
        
    def update_parameters_dictionnary(self):
        self.parameters["CORREC"]           = self.CORREC
        self.parameters["Tsig"]             = self.Tsig
        self.parameters["Asig"]             = self.Asig
        self.parameters["Fsig"]             = self.Fsig
        self.parameters["TimeScale"]        = self.TimeScale
        self.parameters["DEBUG"]            = self.DEBUG
        self.parameters["RAND"]             = self.RAND
        self.parameters["dumpfile"]         = self.dumpfile
        self.parameters["RESTORE"]          = self.RESTORE
        self.parameters["RESTORE_DIR"]      = self.RESTORE_DIR
        self.parameters["Uv"]               = self.Uv
        self.parameters["Tv"]               = self.Tv
        self.parameters["F"]                = self.F
        self.parameters["Mhg"]               = self.Mhg
        self.parameters["HMT"]              = self.HMT
        self.parameters["Tg"]               = self.Tg
        self.parameters["Ug"]               = self.Ug
        self.parameters["COAL"]             = self.COAL
        self.parameters["St"]               = self.St
        self.parameters["Re"]               = self.Re
        self.parameters["Oh"]               = self.Oh
        self.parameters["Nbsup"]            = self.Nbsup
        self.parameters["para_name"]        = self.para_name
        self.parameters["namestud"]         = self.namestud
        self.parameters["name_of_C_file"]   = self.name_of_C_file
        self.parameters["MPI"]              = self.MPI
        self.parameters["Ug"]               = self.Ug
        self.parameters["Drho"]             = self.Drho
        self.parameters["nProc"]            = self.nProc
        self.parameters["Bo"]               = self.Bo              # Bond number
        self.parameters["Ga"]               = self.Ga                # Galileo number
        self.parameters["rho_r"]            = self.rho_r      # ratio of rho
        self.parameters["mu_r"]             = self.mu_r       # ratio between mu
        self.parameters["PHI"]              = self.PHI             # concentration phi
        self.parameters["Nb"]               = self.Nb               # number of bubbles (Must be a perfect square or it will be rounded) 
        self.parameters["ndc"]              = self.ndc               # number of cells per diameters of a bubble
        self.parameters["dtprint"]          = self.dtprint          # writing time step
        self.parameters["dimension"]        = self.dimension      
        self.parameters["g"]                = self.g
        self.parameters["D"]                = self.D
        self.parameters["rho_d"]            = self.rho_d
        self.parameters["LDn"]              = self.LDn
        self.parameters["N"]                = self.N
        self.parameters["Ls"]               = self.Ls
        self.parameters["mu_d"]             = self.mu_d
        self.parameters["mu_f"]             = self.mu_f
        self.parameters["LEVEL"]            = self.LEVEL
        self.parameters["rho_f"]            = self.rho_f
        self.parameters["sig"]              = self.sig
        self.parameters["DP"]               = self.DP
        self.parameters["r1"]               = self.r1
        self.parameters["RND"]              = self.RND
        self.parameters["TMAX"]             = self.TMAX
        self.parameters["FIRSTSTEP"]        = self.FIRSTSTEP
        self.parameters["nopara"]           = self.nopara
        self.parameters["SAVESTEP"]         = self.SAVESTEP
        self.parameters["MOVIES"]           = self.MOVIES
        self.parameters["ndcmin"]           = self.ndcmin
        self.parameters["Tsed"]             = self.Tsed
        self.parameters["Used"]             = self.Used
        self.parameters["Dx"]               = self.Dx
        
    def print_parameters(self):
        with open('parameters.h','w') as the_file:
            for key,value in self.parameters.items():
                if(type(value) != str and type(value) != bool  ):
                    the_file.write('#define '+key+' '+str(round(value,7))+'\n')
                elif(type(value) == str):
                    the_file.write('const char * '+key+'= "'+str(value)+'";\n')
        with open('parameters.csv','w') as the_file:
            for key,value in self.parameters.items():
                if(type(value) != str and type(value) != bool  ):
                    the_file.write(f'{key},')
            the_file.write('\n')
            for key,value in self.parameters.items():
                if(type(value) != str and type(value) != bool  ):
                    the_file.write(str(round(value,7))+',')
            the_file.write('\n')
                    
        with open('parameters.py','w') as the_file:
            the_file.write('parameters = {\n')
            for key,value in self.parameters.items():
                if(type(value) != str):
                    the_file.write('"'+key+'":'+str(value)+',\n')
                else:    
                    the_file.write('"'+key+'":"'+str(value)+'",\n')
            the_file.write('}')
                    
    def print_parameters_after_sim(self):
        with open(self.name_of_C_file+'/parameters.h','w') as the_file:
            for key,value in self.parameters.items():
                if(type(value) != str and type(value) != bool  ):
                    the_file.write('#define '+key+' '+str(round(value,5))+'\n')
        with open(self.name_of_C_file+'/parameters.py','w') as the_file:
            the_file.write('parameters = {\n')
            for key,value in self.parameters.items():
                if(type(value) != str):
                    the_file.write('"'+key+'":'+str(value)+',\n')
                else:    
                    the_file.write('"'+key+'":"'+str(value)+'",\n')
            the_file.write('}')
            
        with open(self.name_of_C_file+'/parameters.csv','w') as the_file:
            for key,value in self.parameters.items():
                if(type(value) != str and type(value) != bool  ):
                    the_file.write(f'{key},')
            the_file.write('\n')
            for key,value in self.parameters.items():
                if(type(value) != str and type(value) != bool  ):
                    the_file.write(str(round(value,7))+',')
            the_file.write('\n')
        
    def cp_dir_and_files(self,P):
        namedir_simulation = self.namedir+self.para_name+'_'+str(P)
        try:
            shutil.rmtree(namedir_simulation)
            os.mkdir(namedir_simulation)
        except:
            os.mkdir(namedir_simulation)
        os.system('cp -r --backup=simple '+self.name_of_C_file+'/*'+' '+namedir_simulation)
        print(os.listdir(''+self.name_of_C_file+'/.'))
        
