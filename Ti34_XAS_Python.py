import time

from scipy import sparse
import numpy as np
from itertools import combinations_with_replacement
from math import factorial
from scipy.sparse.linalg import eigsh
from numpy.linalg import eigh
from time import perf_counter
import pickle
import os
import glob
from scipy.sparse.linalg import LinearOperator
from Broaden_Spectrum import GetBroadeningList, BroadenGamma, BroadenSigma
import matplotlib.pyplot as plt
from time import perf_counter
from KK_And_Merge import *
import copy

global OpsTi

with open("Ti34OpsPython" + '.pkl', 'rb') as f:
    OpsTi = pickle.load(f)


def Lanczos(HS, v=None, m=100):
    """
    Purpose: Perform the Lanczos Algorithm
        :param HS: A matrix in sparse format (csr, csc, bsr, COO, lil, diag)
        :param v: An arbitrary vector of dimensions 1xn where it's modulus is 1
        :param m: The number of iterations
    :return:
    """
    # a parameter used to determine if all imaginary values are small enough to consider zero
    successful = True

    # a check to make sure that the input Hamiltonian is square
    if (HS.get_shape()[0] != HS.get_shape()[1]):
        raise Exception('Must be a square matrix.')

    # get the dimensions of HS
    n = HS.get_shape()[0]

    # creates a random vector 'v' if one is not inputted into the function
    if v is None:
        v = np.random.random(n)
        v = v/float(np.linalg.norm(v))

    # Checks if the number of iterations is smaller or equal to the dimensions of the matrix
    if m > n:
        m = n

    # this performs the dot product for the matrix HS and an arbitrary vector 'v'
    def MV(v):
        return HS.dot(v)

    # creates a linear operator to perform the dot product
    MyOp = LinearOperator((n, n), matvec=MV)

    # This section performs the Lanczos Algorithm for 'm' iterations
    # - 'V' stores the vectors 'v' at each iteration
    # - 'T' stores the elements of the tridiagonal matrix created by the Lanczos Algorithm
    # - 'beta' are the off diagonal elements of the tridiagonal matrix
    # - 'vo' keeps track of the previous 'v' for computation
    # - 'alpha' are the diagonal elements of the tridiagonal matrix
    # - 'w' a intermediate value used to calculate 'beta' and 'alpha'
    # initialization to start the algorithm
    V = np.zeros((m,n), dtype=complex)
    T = np.zeros((m,m), dtype=complex)
    beta = 0
    vo = np.zeros(n)

    # for loop only performs all the calculation, except for the last iteration 'm'
    for j in range(m-1):
        # gets all the elements for the jth component
        w = MyOp.matvec(v)
        alpha = np.dot(np.conjugate(w),v)
        w = w - alpha*v - beta*vo
        beta = np.linalg.norm(w)

        vo = v
        v = w/beta

        # inputs the elements into the proper position in the tridiagonal matrix
        T[j,j]  = alpha
        T[j,j+1] = beta
        T[j+1, j] = beta
        V[j,:] = v

        # A check to determine if the imaginary component is to large and cannot be considered zero
        if abs(beta.real - beta.imag)<beta.real - 1e-7 or abs(alpha.real - alpha.imag) < alpha.real - 1e-7:
            successful = False

    # performs the last iteration of the algortihm (avoids having to add an if statement in the for loop)
    w = MyOp.matvec(v)
    alpha = np.dot(np.conjugate(w),v)
    w = w - alpha*v - beta*vo
    T[m-1, m-1] = alpha
    V[m-1, :] = w/np.linalg.norm(w)

    # A check to determine if the imaginary component is too large and cannot be considered zero
    if abs(beta.real - beta.imag) < beta.real - 1e-7 or abs(alpha.real - alpha.imag) < alpha.real - 1e-7:
        successful = False

    # If an imaginary component was not small enough then let the user know
    if not(successful):
        raise Exception('The imaginary components are too large compared to the real components.')

    return np.real(T), V

def CreateXAS(v, E, Hf, Tmat, NIter=100, Sticks=True, Gamma=0.2, Sigma=0.0):
    """
    Purpose: Return a range of spectrums
        :param v: A list of initial state eigenvectors
        :param E: A list of energies corresponding to those eigenvectors 'v'
        :param Hf: Final state Hamiltonian in sparse format
        :param Tmat: A list of transition matrices
        :param NIter: Number of lanczos iterations
        :param Sticks: Do not broaden if True, broaden if False
        :param Gamma: Gamma value
        :param Sigma: Sigma value, not used if zero is given
    :return Spec: A numpy list of form Spec[m,n]
                - the index m refers to a spectrum which is made up of two elements [spectrum, [i, j]]
                - the index n refers to one of the two elements
                    - n=0 gives a list of the energies and Intensities
                        - Spec[0,0] = spectrum
                        -spectrum[:, 0] = Energy
                        -spectrum[:, 1] = Intensity
                    - n=1 gives an array of two numbers [i, j]
                        -the first number corresponds to the eigenvector, v
                        -the second number refers to the transition matrix, T
                        - the numbers corresponds to the v and T used as inputted in their respective lists
                        -Example: [v1, v2, v3] and [T1, T2, T3]
                        -Spec[0, 1] = [0, 0] corresponds to the spectrum using v1 and T1
    """
    # this is under the assumption that all vectors v are of the same length, as they must for the Lanczos Algorithm
    #tmp = [v[0,i] for i in range(len(v.T))]
    #[evec[0,i] for i in range(len(evec))]
    #print(tmp)
    #print(len(v))

    n = len(v[0])
    #print("n", n)
    row = np.arange(n, dtype=int)
    col = np.zeros(n, dtype=int)
    if not(Sticks):
        gamma = []
        gamma.append([])  # energy values
        gamma.append([])  # fwhm's
        gamma.append([])  # artcan widths interpolating between fwhm's
        gamma[0] = [-9999, 9999]
        gamma[1] = [
            Gamma]  # gamma[1][j] is fwhm to use up to gamma[0][j], where gamma changes to gamma[1][j+1] by interpolating using arctan of width gamma[2][j]
        gamma[2] = [Gamma / 10, Gamma / 10]  # first one is for the -9999 "interface"
        if Sigma != 0:
            # Use a constant Gaussian Broadening
            sigma = []
            sigma.append([])  # energy values
            sigma.append([])  # fwhm's
            sigma.append([])  # artcan widths interpolating between fwhm's
            sigma[0] = [-9999, 9999]
            sigma[1] = [
                Sigma]  # gamma[1][j] is fwhm to use up to gamma[0][j], where gamma changes to gamma[1][j+1] by interpolating using arctan of width gamma[2][j]
            sigma[2] = [Sigma / 10, Sigma / 10]  # first one is for the -9999 "interface"

    # iterate through each eigenvectors and transition matrix
    Spec = []
    
    idx = 0
    for i in range(len(v)):
        Spec.append([])
        # transform vector into sparse form and determine the complex conjugate
        #print(v[i], v[i,:], row, col)
        
        psi = sparse.csr_matrix((v[i].ravel(),(row, col)), shape=(n, 1))
        psiT = complex_conjugate(psi)
        
        #psi = v[i]
        #psiT = np.conjugate(psi.T)
        
        # select the corresponding energy
        Ei = E[i]

        for j in range(len(Tmat)):
            Spec[i].append([])
            # selects transition matrix and determines it's complex conjugate
            T = Tmat[j]
            TT = complex_conjugate(T)

            # constant required in multiple calculations
            x = psiT.dot(T.dot(TT.dot(psi)))
            x = x.data[0]

            # determines the starting vector for the Lanczos algorithm
            phi = (TT.dot(psi).toarray())/np.sqrt(x)
            phi = phi[:, 0]

            #print("Phi")
            #print(phi)
            # performs the Lanzcros Algorithm
            H, V = Lanczos(Hf, v=phi, m=NIter)

            # determines the eigenvalues and eigenvectors of the tridiagonal matrix
            eval, evec = np.linalg.eigh(H)

            # calculate the intensity
            Intensity = np.array([evec[0, i]**2 for i in range(len(evec))])*x
            # determine the energy relative to the intial energy
            energy = eval - Ei

            stick = np.zeros((len(energy), 2))
            stick[:, 0] = energy
            stick[:, 1] = Intensity

            # determines if broadening should be used

            if not(Sticks):

                gammaVals = GetBroadeningList(gamma, stick[:,0])
                EnergyRange = np.linspace(energy[0]-5, energy[-1]+5, 500)
                specBroadG = BroadenGamma(EnergyRange, stick, gammaVals)

                # determines if gaussian broadening needs to be used
                if Sigma == 0:
                    spectrum_i = specBroadG
                else:
                    sigmaVals = GetBroadeningList(sigma, stick[:,0])
                    specBroadS = BroadenSigma(EnergyRange, stick, sigmaVals)
                    sigmaVals2 = GetBroadeningList(sigma, specBroadS[:,0])
                    spectrum_i = BroadenSigma(EnergyRange, specBroadG, sigmaVals2)
                Spec[i][j] = spectrum_i
                #Spec[idx] = [spectrum_i, [i, j]]
            else:
                Spec[i][j] = stick
                #Spec[idx] = [stick, [i,j]]
                # append each spectrum to a list
            idx += 1


    return Spec
    
    

def CreateXASExact(v, E, Hf, Tmat, Sticks=True, Gamma=0.2, Sigma=0):
    """
    Purpose: Return a range of spectrums
        :param v: A list of initial state eigenvectors
        :param E: A list of energies corresponding to those eigenvectors 'v'
        :param Hf: Final state Hamiltonian in sparse format
        :param Tmat: A list of transition matrices
        :param NIter: Number of lanczos iterations
        :param Sticks: Do not broaden if True, broaden if False
        :param Gamma: Gamma value
        :param Sigma: Sigma value, not used if zero is given
    :return Spec: A numpy list of form Spec[m,n]
                - the index m refers to a spectrum which is made up of two elements [spectrum, [i, j]]
                - the index n refers to one of the two elements
                    - n=0 gives a list of the energies and Intensities
                        - Spec[0,0] = spectrum
                        -spectrum[:, 0] = Energy
                        -spectrum[:, 1] = Intensity
                    - n=1 gives an array of two numbers [i, j]
                        -the first number corresponds to the eigenvector, v
                        -the second number refers to the transition matrix, T
                        - the numbers corresponds to the v and T used as inputted in their respective lists
                        -Example: [v1, v2, v3] and [T1, T2, T3]
                        -Spec[0, 1] = [0, 0] corresponds to the spectrum using v1 and T1
    """
    # this is under the assumption that all vectors v are of the same length, as they must for the Lanczos Algorithm
    n = len(v[0])
    row = np.arange(n, dtype=int)
    col = np.zeros(n, dtype=int)
    if not(Sticks):
        gamma = []
        gamma.append([])  # energy values
        gamma.append([])  # fwhm's
        gamma.append([])  # artcan widths interpolating between fwhm's
        gamma[0] = [-9999, 9999]
        gamma[1] = [
            Gamma]  # gamma[1][j] is fwhm to use up to gamma[0][j], where gamma changes to gamma[1][j+1] by interpolating using arctan of width gamma[2][j]
        gamma[2] = [Gamma / 10, Gamma / 10]  # first one is for the -9999 "interface"
        if Sigma != 0:
            # Use a constant Gaussian Broadening
            sigma = []
            sigma.append([])  # energy values
            sigma.append([])  # fwhm's
            sigma.append([])  # artcan widths interpolating between fwhm's
            sigma[0] = [-9999, 9999]
            sigma[1] = [
                Sigma]  # gamma[1][j] is fwhm to use up to gamma[0][j], where gamma changes to gamma[1][j+1] by interpolating using arctan of width gamma[2][j]
            sigma[2] = [Sigma / 10, Sigma / 10]  # first one is for the -9999 "interface"

    # iterate through each eigenvectors and transition matrix
    #Spec = [[] for i in range(len(v)*len(Tmat))]
    idx = 0
    #print("length v: ", len(v))
    Spec = []
    
    idx = 0
    for i in range(len(v)):
        Spec.append([])
        # transform vector into sparse form and determine the complex conjugate
        psi = sparse.csr_matrix((v[i],(row, col)), shape=(n, 1))

        psiT = complex_conjugate(psi)
        # select the corresponding energy
        Ei = E[i]

        for j in range(len(Tmat)):
            Spec[i].append([])
            # selects transition matrix and determines it's complex conjugate
            T = Tmat[j]
            TT = complex_conjugate(T)

            # constant required in multiple calculations
            #x = psiT.dot(T.dot(TT.dot(psi)))
            #x = x.data[0]

            
            # determines the starting vector for the Lanczos algorithm
            phiT = (TT.dot(psi).toarray())#/np.sqrt(x)
            #phiT = phiT[:, 0]

            
            valsf,vecsf = np.linalg.eigh(Hf.todense())

            
            # calculate the intensity
            
            Intensity = np.square(np.dot(vecsf.T,phiT))
            
            # determine the energy relative to the intial energy
            energy = valsf - Ei

            stick = np.zeros((len(energy), 2))
            stick[:, 0] = energy
            stick[:, 1] = Intensity.flatten()

            # determines if broadening should be used

            if not(Sticks):

                gammaVals = GetBroadeningList(gamma, stick[:,0])
                EnergyRange = np.linspace(np.min(energy)-5, np.max(energy)+5, 500)
                specBroadG = BroadenGamma(EnergyRange, stick, gammaVals)

                # determines if gaussian broadening needs to be used
                if Sigma == 0:
                    spectrum_i = specBroadG
                else:
                    sigmaVals = GetBroadeningList(sigma, stick[:,0])
                    specBroadS = BroadenSigma(EnergyRange, stick, sigmaVals)
                    sigmaVals2 = GetBroadeningList(sigma, specBroadS[:,0])
                    spectrum_i = BroadenSigma(EnergyRange, specBroadG, sigmaVals2)

                Spec[i][j] = spectrum_i
            else:
                Spec[i][j] = stick
                # append each spectrum to a list
            idx += 1


    return np.array(Spec)
    
    
def complex_conjugate(M):
    """
    Purpose: Take the transpose and complex conjugate of the sparse matrix
        :param M: Input sparse matrix in csr_matrix format
        :return M_ct: The transposed and complex conjugate in csr_matrix form
    """
    M.check_format(full_check=True)# checks to make sure the format is correct
    M_t = M.transpose(copy=True)
    M_ct = M_t.conjugate(copy=True)
    return M_ct



def save_obj(obj, name ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open( name + '.pkl', 'rb') as f:
        return pickle.load(f)
        
def txtimport(textfile):
    """
    Purpose: To read in files where the first line states the matrix dimensions and convert into sparse matrix
    :param textfile: The file to be converted
    :return: The sparse matrix
    """
    M = np.loadtxt(textfile, skiprows=1) # reads in the file
    
    if(len(M.shape) < 2 and len(M) > 0):
      #print(M)
      M = M.reshape(1,3)
      #print(M)
    dim = np.loadtxt(textfile, max_rows=1)
    row = int(dim[0])
    col = int(dim[1])
    if(len(M) == 0):
      M = np.array([[1,1,0]])
#      return sparse.csr_matrix((0, (0, 0)), [row,col])
#    else:
    
    return sparse.csr_matrix((M[:,2], (M[:,0]-1, M[:,1]-1)), [row,col])

def WriteOpsToPickle():
  Ops = dict()
  for files in glob.glob("*.dat"):
    p1 = os.path.splitext(files)[0]
    print(p1)
    Ops[p1] = txtimport(files)

  save_obj(Ops,"Ti34OpsPython")



def MergeWithOffRes(inSpec,EShift,S1,S2,m,b,E1,w1,E2,w2,offresFile,element,c1,c2,fi1,fi2,ff1,ff2):

  mySpec = {}
  edgeJump = {}
  

  inSpec[:,0] = inSpec[:,0] + EShift

  edgeJump = S2*((np.arctan(w1*(inSpec[:,0]-E1))/(np.pi)+1/2)+ 0.5*(arctan(w2*(inSpec[:,0]-E2))/(math.pi)+1/2))

  inSpec[:,1] = S1*inSpec[:,1]*inSpec[:,0] + m*(inSpec[:,0]-inSpec[0,0]) + b + edgeJump


  #Subtract gaussian in postedge for exafs-like feature
  gw = 5
  ga = 25
  xc = 468
  ep = 4*math.log(2)/gw/gw
  den = gw * math.sqrt(math.pi/4/math.log(2))

  inSpec[:,1] = inSpec[:,1] - ga/den*np.exp(-ep*(inSpec[:,0]-xc)**2)


  
  
  #plt.plot(inSpec[:,0],inSpec[:,1])
  #plt.show()
  
  
  
  
  
 #Example call:  python KK_And_Merge.py spec.dat Ti.ff Ti 432 437
    
  offres = GetSpecFromFile(offresFile)    
  offres = offres.T
  
  
  mergedE = np.hstack((offres[np.where(offres[:,0] < c1)[0],0],inSpec[np.where((inSpec[:,0]>c1) & (inSpec[:,0]<c2))[0],0],offres[np.where(offres[:,0]>c2)[0],0]))
  
  mergedf2 = np.hstack((offres[np.where(offres[:,0] < c1)[0],2],inSpec[np.where((inSpec[:,0]>c1) & (inSpec[:,0]<c2))[0],1],offres[np.where(offres[:,0]>c2)[0],2]))
   
  
  
  #print(offres.shape)
  
  #spec = MergeOffRes(inSpec,offres,c1,c2)

  offset = GetOffset(element) #'Ti'

  mergedf1 = KK_Robert(mergedE, mergedf2, offset)

  return np.transpose(np.vstack((mergedE,mergedf1,mergedf2)))
'''
  f = open('merged.ff','w')
  for i in range(len(mergedE)):
      f.write(str(mergedE[i]))
      f.write(' ')
      f.write(str(mergedf1[i]))
      f.write(' ')
      f.write(str(mergedf2[i]))
      f.write(' ')
      f.write('\n')
  f.close()
 '''
  




#Only run once to get an operator file
#WriteOpsToPickle()




def GetTiFormFactor(dExy,dExzyz,dEx2y2,dEz2, nd=1,T=300,tenDq=2.12):

#prepath = "ff/"
#OrbE = np.loadtxt(prepath + "OrbitalEnergies.txt")
  dExy = float(dExy)
  dExzyz = float(dExzyz)
  dEx2y2 = float(dEx2y2)
  dEz2 = float(dEz2)
  #Ops = load_obj("Ti34OpsPython")
  Ops = copy.deepcopy(OpsTi)
  
  KelvinToeV = 8.61735E-5
  T = T * KelvinToeV
  
  if(nd == 0):
    F2dd = 0
    F4dd = 0
    F2ddX = 0
    F4ddX = 0
    F2pdX = 3.628
    G1pdX = 3.153*1.0
    G3pdX = 1.792*1.0
    zeta_2p = 3.820#--3.738
    zeta_3d = 0.019
    zeta_3dX = 0.032
  
  else:
    F2dd = 0
    F4dd = 0
    F2ddX = 5.896*1.0018
    F4ddX = 3.704*1.0018
    F2pdX = 3.181*1.0018
    G1pdX = 3.193*0.85*1.25
    G3pdX = 1.814*0.85*1.25
    zeta_2p = 3.71
    zeta_3d = 0.019
    zeta_3dX = 0.026


  F0dd    = (F2dd+F4dd)*2/63
  F0ddX    = (F2ddX+F4ddX)*2/63
  F0pdX    =  G1pdX*1/15 + G3pdX*3/70

  Bz      = 0.000001


  H = (0.6*tenDq + dEz2) * Ops["p6d" + str(nd) + "_Oppz2"] + \
      (0.6*tenDq + dEx2y2) * Ops["p6d" + str(nd) + "_Oppx2y2"] + \
      (-0.4*tenDq + dExy) * Ops["p6d" + str(nd) + "_Oppxy"] + \
      (-0.4*tenDq + dExzyz) * Ops["p6d" + str(nd) + "_Oppxzyz"] + \
      Bz * (2*Ops["p6d" + str(nd) + "_OppSz"] + Ops["p6d" + str(nd) + "_OppLz"])

  HXAS = (0.6*tenDq + dEz2) * Ops["p5d" + str(nd+1) + "_Oppz2"] + \
        (0.6*tenDq + dEx2y2) * Ops["p5d" + str(nd+1) + "_Oppx2y2"] + \
        (-0.4*tenDq + dExy) * Ops["p5d" + str(nd+1) + "_Oppxy"] + \
        (-0.4*tenDq + dExzyz) * Ops["p5d" + str(nd+1) + "_Oppxzyz"] + \
        Bz * (2*Ops["p5d" + str(nd+1) + "_OppSz"] + Ops["p5d" + str(nd+1) + "_OppLz"]) + \
        F0ddX * Ops["p5d" + str(nd+1) + "_OppF0"] + F2ddX * Ops["p5d" + str(nd+1) + "_OppF2"] + F4ddX * Ops["p5d" + str(nd+1) + "_OppF4"] +\
         F0pdX * Ops["p5d" + str(nd+1) + "_OppUpdF0"] + F2pdX * Ops["p5d" + str(nd+1) + "_OppUpdF2"] + G1pdX * Ops["p5d" + str(nd+1) + "_OppUpdG1"] + G3pdX * Ops["p5d" + str(nd+1) + "_OppUpdG3"] +\
         + zeta_3dX * Ops["p5d" + str(nd+1) + "_Oppzeta3d"] + zeta_2p * Ops["p5d" + str(nd+1) + "_Oppzeta2p"]  


  if(nd == 0):
    row, column = HXAS.get_shape()
    I1 = sparse.identity(row, format='csr')  # creates the identity matrix for H1
    t2gShift = -0.25
    HXAS = HXAS +  (3.75*I1 - Ops["p5d" + str(nd+1) + "_OppJsqr_2p"])/3 * Ops["p5d" + str(nd+1) + "_OppNt2g"] * t2gShift



  eval, evec = eigh(H.todense())



  #this seems to do nothing?
  tmp = np.zeros((len(eval),len(eval)))
  for i in range(len(eval)):
    for j in range(len(eval)):
      tmp[i,j] = evec[i,j] 
  evec = tmp
 

  Egrnd = eval[0]
  dEList = eval - Egrnd
  myEps = 1e-6
  
  dE = -np.log(myEps)*T

  v = evec[:,np.where(dEList < dE)[0]].T
  
  E = eval[np.where(dEList < dE)]


  dz = np.exp(-(E-eval[0])/T)
  Z = np.sum(dz)
  dz = dz / Z
  
  Tmat = [Ops["p5d" + str(nd+1) + "_TXASx"],Ops["p5d" + str(nd+1) + "_TXASz"]]#, T, T]



  
  myNIter = 25
  if(nd == 1):
    myNIter = 55
    
  spec = CreateXAS(v, E, HXAS, Tmat, NIter=myNIter, Sticks=True, Gamma=0.15, Sigma=0.15)
  #spec = CreateXASExact(v, E, HXAS, Tmat, Sticks=False, Gamma=0.1, Sigma=0.1)
  
  
  
  #[[spec[i][j][:,1] = dz[i]*spec[i][j][:,1] for i in range(len(dz))] for j in range(len(Tmat))]
  for i in range(len(dz)):
    for j in range(len(Tmat)):
      spec[i][j][:,1] = dz[i]*spec[i][j][:,1]
  
  
  
  #plt.plot(spec[0][0][:,0],spec[0][0][:,1])
  #plt.plot(spec[0][1][:,0],spec[0][1][:,1])
  #plt.show()
  eners_xx = np.hstack([spec[i][0][:,0] for i in range(len(spec))])
  eners_zz = np.hstack([spec[i][1][:,0] for i in range(len(spec))])
  sticks_xx = np.hstack([spec[i][0][:,1] for i in range(len(spec))])
  sticks_zz = np.hstack([spec[i][1][:,1] for i in range(len(spec))])
  
  #print(len(spec[0][0][:,0]))
  #print(len(eners_xx))
  
  spec_xx = np.transpose(np.vstack((eners_xx,sticks_xx)))
  spec_zz = np.transpose(np.vstack((eners_zz,sticks_zz)))

  
    
    
  gamma = []
  gamma.append([])  # energy values
  gamma.append([])  # fwhm's
  gamma.append([])  # artcan widths interpolating between fwhm's
  
  if(nd == 0):
    gamma[0] = [ -9999,    -0.900,  2.000,  4.750,  6.500,  9999]
    gamma[1] = [            0.064,  0.560,  1.000,  1.360,  0.888]  # gamma[1][j] is fwhm to use up to gamma[0][j], where gamma changes to gamma[1][j+1] by interpolating using arctan of width gamma[2][j]
    gamma[2] = [0.3,        0.300,  0.300,  0.300,  0.300,  0.300]  # first one is for the -9999 "interface"
  else:
    gamma[0] = [ -9999,     0.500,  2.500,  9999]
    gamma[1] = [            0.468,  0.680,  1.100]  # gamma[1][j] is fwhm to use up to gamma[0][j], where gamma changes to gamma[1][j+1] by interpolating using arctan of width gamma[2][j]
    gamma[2] = [0.3,        0.300,  0.300,  0.300]  # first one is for the -9999 "interface"
  
  
  sigma = []
  sigma.append([])  # energy values
  sigma.append([])  # fwhm's
  sigma.append([])  # artcan widths interpolating between fwhm's
  if(nd == 0):
    sigma[0] = [  -9999,     -0.750,  9999]
    sigma[1] = [        0.200, 0.350]  # gamma[1][j] is fwhm to use up to gamma[0][j], where gamma changes to gamma[1][j+1] by interpolating using arctan of width gamma[2][j]
    sigma[2] = [0.300,  0.300, 0.300]  # first one is for the -9999 "interface"
  else:
    sigma[0] = [  -9999,     9999]
    sigma[1] = [        0.25  ]  # gamma[1][j] is fwhm to use up to gamma[0][j], where gamma changes to gamma[1][j+1] by interpolating using arctan of width gamma[2][j]
    sigma[2] = [0.300,  0.300]  # first one is for the -9999 "interface"


  
  epad = 12
  estep = 0.04
  ehigh = np.max(eners_xx) + epad
  elow  = np.min(eners_xx) - epad
  nE = int((ehigh - elow)/estep)
  erange = np.linspace(np.min(eners_xx) - epad,np.max(eners_xx)+epad, nE) 
  gammaVals = GetBroadeningList(gamma, erange)
  sigmaVals = GetBroadeningList(sigma, spec_xx[:,0])
  sigmaValszz = GetBroadeningList(sigma, spec_zz[:,0])
  


  spec_xx = BroadenGamma(erange, BroadenSigma(erange, spec_xx, sigmaVals), gammaVals)
  spec_zz = BroadenGamma(erange, BroadenSigma(erange, spec_zz, sigmaValszz), gammaVals)
  #Spectra_xx.Broaden(,)
  
  #print(np.sum(spec_xx[:,1]), np.sum(spec_zz[:,1]))
  
  #plt.plot(spec_xx[:,0],spec_xx[:,1],label="xx")
  #plt.plot(spec_zz[:,0],spec_zz[:,1],label="zz")
  #plt.legend()
  #plt.show()
  
  
  EShift = 460.35-0.49
  scale = 18.0
  edgeShift = 4.0+2.2-2.2
  E1 = 453.8 + edgeShift
  w1 = 2*0.5/0.8
  E2 = 460.2 + edgeShift
  w2 = 2*0.5/0.8
  offResFile = "Ti.ff"
  element = "Ti"
  c1 = 443 #432
  c2 = 479#477
  fi1 = 435
  fi2 = 452
  ff1 = 474
  ff2 = 487

  S1 = (0.27-0.015)*152*.55
  S2 = 8.95
  m = -0.014
  b = 1.6155-0.005
  
  if(nd == 1):
    EShift = 459.45-0.49
    m = -0.013
    S2 = 8.95

  merged_xx = MergeWithOffRes(spec_xx,EShift,S1,S2,m,b,E1,w1,E2,w2,offResFile,element,c1,c2,fi1,fi2,ff1,ff2)
  merged_zz = MergeWithOffRes(spec_zz,EShift,S1,S2,m,b,E1,w1,E2,w2,offResFile,element,c1,c2,fi1,fi2,ff1,ff2)

  final = np.hstack((merged_xx,merged_zz[:,1:]))
  
  
  #write to file

  #np.savetxt(saveName,final[np.where((final[:,0] > 300) & (final[:,0] < 1000))])

  return final[np.where((final[:,0] > 300) & (final[:,0] < 1000))]
  
  
  
  
  
  
  
#Example call

#GetTiFormFactor(1,300,2.1,0.1,-0.1,0.1,-0.14)

