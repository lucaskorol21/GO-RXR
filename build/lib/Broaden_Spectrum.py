import numpy as np
import matplotlib.pyplot as plt


 
def GetBroadeningList(regions,inSpec):
    #determine point by point fwhm by interpolating using arctan
    gVals = np.zeros((len(inSpec),2))
    gVals[:,0] = inSpec
    
    for j in range(len(regions[1])):
      gVals[:,1] = gVals[:,1] + (np.arctan((gVals[:,0]-regions[0][j])/regions[2][j])/np.pi+0.5) * (np.arctan((regions[0][j+1]-gVals[:,0])/regions[2][j+1])/np.pi+0.5) * regions[1][j]
    
    return gVals

def BroadenGamma(eVals,inSpec,GammaList):
    #Lorentzian broadening of inSpec with FWHM as per the GammaList
    specBroad = np.zeros((len(eVals),2))
    specBroad[:,0] = eVals #energies
    for j in range(len(inSpec[:,0])):
        specBroad[:,1] = specBroad[:,1] + inSpec[j,1]/((eVals-inSpec[j,0])*(eVals-inSpec[j,0])+GammaList[j,1]*GammaList[j,1]/4)*(GammaList[j,1]/np.pi/2)*(eVals[1]-eVals[0])
        
    return specBroad

    
def BroadenSigma(eVals,inSpec,SigmaList):
    #Gaussian broadening of inSpec with FWHM as per the SigmaList
    specBroad = np.zeros((len(eVals),2))
    specBroad[:,0] = eVals #energies
    cfact = 2*np.sqrt(2*np.log(2))

    for j in range(len(inSpec[:,0])):
        specBroad[:,1] = specBroad[:,1] + inSpec[j,1]/SigmaList[j,1]*cfact/np.sqrt(2*np.pi)*np.exp(-(eVals-inSpec[j,0])**2/(2*SigmaList[j,1]*SigmaList[j,1]/cfact/cfact)) * (eVals[1] - eVals[0])
        
    return specBroad
    





 
#EXAMPLE OF USAGE

'''
#Use a Lorentzian broadening that varies with energy
gamma = []
gamma.append([]) #energy values
gamma.append([]) #fwhm's
gamma.append([]) #artcan widths interpolating between fwhm's
gamma[0] = [-9999,-5.0,0.0,5.0,9999]
gamma[1] = [0.2,0.6,0.2,0.8] #gamma[1][j] is fwhm to use up to gamma[0][j], where gamma changes to gamma[1][j+1] by interpolating using arctan of width gamma[2][j]
gamma[2] = [0.2,0.2,0.2,0.2,0.2] #first one is for the -9999 "interface"


#Use a constant Gaussian Broadening
sigma = []
sigma.append([]) #energy values
sigma.append([]) #fwhm's
sigma.append([]) #artcan widths interpolating between fwhm's
sigma[0] = [-9999,9999]
sigma[1] = [0.55] #gamma[1][j] is fwhm to use up to gamma[0][j], where gamma changes to gamma[1][j+1] by interpolating using arctan of width gamma[2][j]
sigma[2] = [0.1,0.1] #first one is for the -9999 "interface"



#Generate a test spectrum
#Just has four peaks
testSpec = np.zeros((4,2))
testSpec[0,0] = -3.0; testSpec[0,1] = 0.50
testSpec[1,0] = -1.5; testSpec[1,1] = 1.00
testSpec[2,0] =  1.5; testSpec[2,1] = 0.75
testSpec[3,0] =  3.0; testSpec[3,1] = 0.25

#The energy mesh we want to broaden our spectrum onto
energies = np.linspace(-5,5,500)


#Generate Lor and Gau broadening lists (a FWHM for each data point)
gammaVals = GetBroadeningList(gamma,testSpec)
sigmaVals = GetBroadeningList(sigma,testSpec)

#Ex 1: broadening with lorentzian
specBroadG = BroadenGamma(energies,testSpec,gammaVals)

#Ex 2: broadening with Gaussian
specBroadS = BroadenSigma(energies,testSpec,sigmaVals)

#Ex 3: broadening with both. Need to create a new sigma list, and then
#apply to already gamma-broadened
sigmaVals2 = GetBroadeningList(sigma,specBroadS)
specBroadGS = BroadenSigma(energies,specBroadG,sigmaVals2)


#Plot the broadened spectra
plt.plot(specBroadG[:,0],specBroadG[:,1], label="Lorentzian")
plt.plot(specBroadS[:,0],specBroadS[:,1], label="Gaussian")
plt.plot(specBroadGS[:,0],specBroadGS[:,1], label="Lor + Gau")
plt.show()
'''
