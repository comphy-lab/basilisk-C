#!/usr/bin/env python


# System imports
from   __future__ import print_function
#from   numpy      import *
from   optparse   import *
import sys

#
# Parse the command-line arguments
parser = OptionParser()
parser.add_option("-t", "--testharness", action="store_true",
                  dest="testharness", default=False,
                  help="test local build modules; prevent loading system-installed modules")
parser.add_option("-v", "--verbosity", type="int", dest="verbosity", default=2,
                  help="set the verbosity level [default 2]")
parser.add_option("-n", "--numelem", type="int", dest="numelem", default=20,
                  help="set the number of elements [default 100]")
parser.add_option("--sigma", type="float", dest="sigma", default=0.0,
                  help="set alpha")
parser.add_option("--scale", type="float", dest="scale", default=1.0,
                  help="set scale")
parser.add_option("--iters", type="int", dest="maxNewtonIters", default=100,
                  help="set maxNewtonIters")
parser.add_option("-N", type="int", dest="N", default=64,   # TODO: what is
                  help="set N")                              # correct default?
options,args = parser.parse_args()


#
# Under normal usage, simply use 'from PyTrilinos import Epetra'.  For testing,
# we want to be able to control whether we import from the build directory or
# from a system-installed version of PyTrilinos.
from PyTrilinos import Teuchos   
from PyTrilinos import Epetra    
from PyTrilinos import EpetraExt 
from PyTrilinos import NOX       
from PyTrilinos import LOCA      

import pickle
import numpy as np
import pypg as bas


# ######################################################################
# class for the Jacobian operator
class ChanJac(Epetra.Operator):
    tranpose_ = False
    # these are used for shiftedmatrix
    alpha_ = 1.0
    beta_ = 0.0
    x_jac = []
    epetra_map = []

    def __init__(self, comm, N = 1):
        Epetra.Operator.__init__(self)
        self.__comm  = comm
        self.__label = "Chan_Linear_Operator"
        self.epetra_map = Epetra.Map(N,0,comm)
#        self.fileout = "sol.dat"

    def Comm(self):
        return self.__comm

    def OperatorDomainMap(self):
        return self.epetra_map

    def OperatorRangeMap(self):
        return self.epetra_map

    def SetUseTranspose(self,trans):
        try:
            self.transpose_ = bool(trans)
        except:
            return -1
        return 0

    def UseTranspose(self):
        return self.transpose_

    def Label(self):
        return self.__label

    def Apply(self, x, y):
        try:
            # print('--> alpha, beta = ', self.alpha_,', ',self.beta_)
            (nvec,nvar) = x.shape
            for nv in range(0,nvec):
              bas.pystep_lt(self.x_jac[:],y[nv,:],x[nv,:])
#            y.Multiply(2.0, x, self.x_jac, 0.0)
            y.Update(self.alpha_, y, self.beta_, x,0.0)
            return 0
        except Exception as e:
            print("A python exception was raised in ChanJac>.Apply:")
            print(e)
            return -1

    def ApplyInverse(self, x, y):
        print("NOX Error")
        return 1

    def NormInf(self):
        print("NOX Error")
        return 1

    def HasNormInf(self):
        return False

    def setShMat(self, alpha, beta):
        self.alpha_ = alpha
        self.beta_  = beta

# ######################################################################

class ChanProblemInterface(LOCA.Epetra.Interface.TimeDependent,
                           NOX.Epetra.Interface.Jacobian,
                           LOCA.Epetra.Interface.Required):
  history = []
  Jaco_ = []
  def __init__(self, comm, N, Jacob_ = None):
    # Initialize base class first
    print('--> ChanProblemInterface.__init__, N = ',N)
    LOCA.Epetra.Interface.TimeDependent.__init__(self)
    NOX.Epetra.Interface.Jacobian.__init__(self)

    self.__comm = comm
    self.n = N
    self.sigma_ = 0.0

    vec = np.zeros(N)
    for i in range(0,self.n):
      vec[i] = 0.
#      vec[i] = i*(self.n-1-i)/((self.n-1.)*(self.n-1.)) + 0.001

    map_ = Epetra.Map(N,0,comm)
    self.initGuess = Epetra.Vector(map_,vec)

    self.Jaco_ = Jacob_

  def printSolution(self, x, conParam):
      #self.history.append([conParam,x[1]])
      self.history.append(x[:])
      print('x = ', x, type(x), conParam, self.sigma_)
      with open("sol.dat",'a') as fout:
        np.savetxt(fout,x[:].reshape(1, x.shape[0]))

      with open("cont_diag.dat",'a') as fout:
        np.savetxt(fout,[conParam])

  def computeF(self, x, F, fillFlag):
      # we ignore the fillFlag argument and always perform all computation
      try:
#          n = self.n
#          xx = np.zeros(n)
#          ff = np.zeros(n)
#          xx[:] = x
          bas.pyadjust_contpar(self.sigma_)
          bas.pystep(x[:],F[:])
#          F[:] = ff

          # for i in range(0,n):
          #     F[i] = x[i]**2 - self.sigma_ * float(i)
          return True
      except Exception as e:
          print("Proc", self.__myPID,
                "ChanProblemInterface.computeF() has thrown an exception:")
          print(str(type(e))[18:-2] + ":", e)
          return False

  def InitialGuess(self):
      return self.initGuess

  def computeJacobian(self, x, Jac):
      try:
          self.Jaco_.x_jac = x.copy()
          self.Jaco_.setShMat(1.0,0.0);
          return True
      except Exception as e:
          print("Proc", self.__myPID, "ChanProblemInterface.computeJacobian() "
                "has thrown an exception:")
          print(str(type(e))[18:-2] + ":", e)
          return False

  def setParams(self, p):
    self.sigma_ = p.getValue("sigma")

  def setParameters(self, p):
    self.sigma_ = p.getValue("sigma")

  def computeShiftedMatrix(self, alpha, beta, x, A):
      pass

  def applyShiftedMatrix(self, alpha, beta, input, result):
      print('--> applyShiftedMatrix')
      self.Jaco_.Apply(input, result)
      result.update(alpha,result,beta,input,0.0);

  def setJacobianOperator(self, Jaco):
      self.Jaco_ = Jaco

  # def setXdot(self, xdot, t):
  #     pass
  #
  # def postProcessContinuationStep(self):
  #     pass
  #
  # def preProcessContinuationStep(self):
  #     pass

######################################################################

# Main routine
def main():


# basilisk options
  Nx = 32
  nl = 16
  x = np.linspace(0, 1, Nx)
  y = np.linspace(0, 1, Nx)
  s = np.linspace(-1, 0, nl)
  X,Y = np.meshgrid(x,y)
  ne = Nx*Nx*nl
  b0 = np.zeros((nl,Nx,Nx))
  db0 = np.zeros((nl,Nx,Nx))

  for l in range(0,nl):
    b0[l,:,:] = np.sin(2*np.pi*X)*np.cos(2*np.pi*X)


  bas.init_grid(Nx)
  bas.pyinit_vertgrid(nl)
  bas.pyset_vars()
  contpar = 1
  bas.pyset_contpar(contpar)
  bas.pystep(b0.reshape(ne),db0.reshape(ne))

  # prepare output file
  fout = open("sol.dat","w+")
  fout.close()
  fout = open("cont_diag.dat","w+")
  fout.close()

  params = {}
  params["N"] = Nx
  params["nl"] = nl
  with open("params.dat",'wb') as fout:
    pickle.dump(params,fout)

  flagjac = 0

  sigma = options.sigma

  # Communicator
  comm    = Epetra.PyComm()
  myPID   = comm.MyPID()
  numProc = comm.NumProc()

  # Suppress 'Aztec status AZ_loss: loss of precision' messages
  comm.SetTracebackMode(0)

  # Get parameters
  scale = options.scale
  maxNewtonIters = options.maxNewtonIters
  N = ne

  #####################################################################
  paramList = \
    LOCA.Epetra.defaultContinuationParameters(comm=comm,
                                              verbosity=options.verbosity)
  # from options_trilinos import paramList

  lsParams    = paramList["NOX"]["Direction"]["Newton"]["Linear Solver" ]
  printParams = paramList["NOX"]["Printing"]
  stepper     = paramList["LOCA"]["Stepper"]
  printParams["Output Information"] = NOX.Utils.StepperIteration + \
                                      NOX.Utils.StepperDetails
  stepper["Max Steps"] = 200
  stepper["Min Value"] = -10.0
  stepper["Max Value"] = 10.0
  stepper["Initial Value"] = 1e-7
  paramList["LOCA"]["Stepper"]["Continuation Parameter"] = "sigma"
  paramList["LOCA"]["Step Size"]["Initial Step Size"] = 1e-7
  paramList["LOCA"]["Step Size"]["Max Step Size"] = 1e-2

  stepper["Compute Eigenvalues"] = False
  stepper["Eigensolver"]["Shift"] = 0.0
  stepper["Eigensolver"]["Num Blocks"] = 20
  stepper["Eigensolver"]["Block Size"] = 4
  stepper["Eigensolver"]["Operator"] = "Jacobian Inverse"
  stepper["Eigensolver"]["Sorting Order"] = "LM"
  p = LOCA.ParameterVector()
  p.addParameter("sigma", 0.0)

  # create Jacobian Operator
  Jaco = ChanJac(comm, N)

  interface = ChanProblemInterface(comm, N, Jacob_ = Jaco)
  interface.setParameters(p)

  soln    = interface.InitialGuess()
  noxSoln = NOX.Epetra.Vector(soln, NOX.Epetra.Vector.CreateView)


  print('--> initial guess, done!')
  iReq = interface
  iJac = interface
  if flagjac == 0:
    iJac = NOX.Epetra.MatrixFree(printParams, interface, noxSoln)
    Jaco = iJac
      
  linSys = NOX.Epetra.LinearSystemAztecOO(printParams,
                                          lsParams,
                                          iReq,
                                          iJac,
                                          Jaco,
                                          soln)

  shiftedLinSys = NOX.Epetra.LinearSystemAztecOO(printParams,
                                                 lsParams,
                                                 iReq,
                                                 iJac,
                                                 Jaco,
                                                 soln)
  # Create the Group
  locaSoln = NOX.Epetra.Vector(soln, NOX.Epetra.Vector.CreateView)

  globalData = LOCA.createGlobalData(paramList)
  print('--> group', end="")
  group = LOCA.Epetra.Group(globalData, printParams, iReq, locaSoln, linSys, p)

  print('done!')

  # Create the convergence tests
  normF = NOX.StatusTest.NormF(1.0e-8)
  maxIters = NOX.StatusTest.MaxIters(maxNewtonIters)
  converged = NOX.StatusTest.Combo(NOX.StatusTest.Combo.OR, normF, maxIters)

  # Create the stepper
  print("- Making stepper")
  stepper = LOCA.Stepper(globalData, group, converged, paramList)

  # Perform continuation run
  print("- About to run stepper")

  status = stepper.run()
  print("- Completed stepper run")

  if 1 == 0:
      import pylab as plt

      lastpt = np.array(interface.history)[-1,:]

      b0 = np.reshape(lastpt,(nl,Nx,Nx))

      print(np.array(interface.history).shape)

      plt.subplot(121)
      plt.contourf(X,Y,b0[-1,:,:].T,20); plt.colorbar()
      plt.subplot(122)
      plt.contourf(x,s,b0[:,2,:],20); plt.colorbar()
      plt.show()
      # plt.subplot(121)
      # plt.plot(array(interface.history)[:,0],array(interface.history)[:,1],'.-')
      # plt.subplot(122)
      # plt.plot(array(interface.history)[:,0],'.-')
      # plt.show()
      print("- plot done!")

  bas.pytrash_vars()

  return status

# ######################################################################

if __name__ == "__main__":

    status = main()

    # if status == LOCA.Abstract.Iterator.NotFinished:
    #   print("Iterator status = Not Finished")
    #   print("End Result: TEST FAILED")
    # elif status == LOCA.Abstract.Iterator.Failed:
    #   print("Iterator status = Failed")
    #   print("End Result: TEST FAILED")
    # elif status == LOCA.Abstract.Iterator.Finished:
    #   print("Iterator status = Finished")
    #   print("End Result: TEST PASSED")
    # elif status == LOCA.Abstract.Iterator.LastIteration:
    #   print("Iterator status = Last Iteration")
    #   print("End Result: TEST PASSED")

sys.exit(status)
