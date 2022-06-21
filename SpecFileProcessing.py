import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, least_squares
import math
import scipy.signal
import pickle
import os
import glob

#read spec file
#returns a list of the scans, where the scans are dictionaries
#At start of file, #O0 and #01 store the motor names for #P0 and #P1 for scan info
#Then #L stores motor names for scan data
#The linear inclination angle is given as one of the p values, so this can distinguish
#LH from LV, but cannot distinguish circular.  So we need to get polarization from pre-scan header
def GetPol(line,inpol):
  if(line.upper().find("CIRCULAR LEFT") > -1):
    return "CL"
  if(line.upper().find("CIRCULAR RIGHT") > -1):
    return "CR"
  if(line.upper().find("LINEAR VERT") > -1):
    return "LV"
  if(line.upper().find("LINEAR HORIZ") > -1):
    return "LH"
  return inpol

def ReadSpecFile(fname):
  specInfo = []
  specInfoNames = ["Polarization"]
  specData = []
  specNum = 0
  specIdx = 0
  curPol = "LV"
  with open(fname) as f:
    fdata = f.readlines()
    iter = 0
    while(iter < len(fdata)):
      iterold = iter
      curPol = GetPol(fdata[iter],curPol)

      fdata[iter] = fdata[iter].replace("Sample Horz","Sample_Horz")
      fdata[iter] = fdata[iter].replace("Sample Vert","Sample_Vert")
      fdata[iter] = fdata[iter].replace("Sample Dept","Sample_Dept")
      fdata[iter] = fdata[iter].replace("MCP A","MCP_A")
      fdata[iter] = fdata[iter].replace("MCP B","MCP_B")
      fdata[iter] = fdata[iter].replace("SDD ROI","SDD_ROI")
      fdata[iter] = fdata[iter].replace("Mech Curr","Mesh_Curr")
      fdata[iter] = fdata[iter].replace("Samp Curr","Samp_Curr")
      fdata[iter] = fdata[iter].replace("#L Energy","#L Energy1")



      #print(fdata[iter])
      if(fdata[iter].startswith("#O")):
        #the following removes the first word (#O0) and turns the rest into a list
        [ specInfoNames.append(x) for x in fdata[iter].split()[1:] ]

      if(fdata[iter].startswith("#S")):
        #print(fdata[iter])
        specNum = int(fdata[iter].split()[1])
        specIdx = specNum - 1
        
        
      if(fdata[iter].startswith("#P")):
        if(fdata[iter].startswith("#P0")):
          specInfo.append(dict())
          specInfo[specIdx][specInfoNames[0]] = curPol
        pdata = []
        [ pdata.append(float(x)) for x in fdata[iter].split()[1:] ]
        for j in range(len(pdata)):
          specInfo[specIdx][specInfoNames[len(specInfo[specIdx])]] = pdata[j]
        #print(specInfo[specIdx])

      if(fdata[iter].startswith("#L")):
        mnames = []
        [ mnames.append(x) for x in fdata[iter].split()[1:] ]
        #now read scan, until blank line or "#C" line
        tdata = []
        iter = iter + 1
        if(iter < len(fdata)):
          curPol = GetPol(fdata[iter],curPol)
        
        while(iter < len(fdata) and len(fdata[iter]) > 3 and ((not fdata[iter].startswith("#")) or ("injection" in fdata[iter]) or ("communication" in fdata[iter]))):
          if(("injection" not in fdata[iter]) and ("communication" not in fdata[iter])):
            tlist = fdata[iter].split()
            tlist = [float(x) for x in tlist]
            tdata.append(tlist)
          iter = iter+1
          if(iter < len(fdata)):
            curPol = GetPol(fdata[iter],curPol)
        tdata = np.array(tdata)

        #now map the keys to the data lists and store in specData[scanIdx]
        specData.append(dict())
        if(len(tdata) > 0): #otherwise scan was aborted before starting
          for i in range(len(mnames)):
            specData[specIdx][mnames[i]] = tdata[:,i]
        else:
          for i in range(len(mnames)):
            specData[specIdx][mnames[i]] = [0]
        #if(iter - iterold > 1):
        #  iter = iter - 1   
          
      iter = iter+1
      
      
  return specInfo,specData
      

def ReadSpecFileRIXS(fname):
  specInfo = []
  specInfoNames = ["Polarization"]
  specData = []
  specNum = 0
  specIdx = 0
  curPol = "LV"
  with open(fname) as f:
    fdata = f.readlines()
    iter = 0
    while(iter < len(fdata)):
      iterold = iter
      curPol = GetPol(fdata[iter],curPol)

      fdata[iter] = fdata[iter].replace("Sample Horz","Sample_Horz")
      fdata[iter] = fdata[iter].replace("Sample Vert","Sample_Vert")
      fdata[iter] = fdata[iter].replace("Sample Dept","Sample_Dept")
      fdata[iter] = fdata[iter].replace("MCP A","MCP_A")
      fdata[iter] = fdata[iter].replace("MCP B","MCP_B")
      fdata[iter] = fdata[iter].replace("SDD ROI","SDD_ROI")
      fdata[iter] = fdata[iter].replace("Mech Curr","Mesh_Curr")
      fdata[iter] = fdata[iter].replace("Samp Curr","Samp_Curr")
      fdata[iter] = fdata[iter].replace("#L Energy","#L Energy1")
      #print(fdata[iter])
      
      if(fdata[iter].startswith("#O")):
        #the following removes the first word (#O0) and turns the rest into a list
        [ specInfoNames.append(x) for x in fdata[iter].split()[1:] ]

      if(fdata[iter].startswith("#S")):
        #print(fdata[iter])
        specNum = int(fdata[iter].split()[1])
        specIdx = specNum - 1
        
        
      if(fdata[iter].startswith("#P")):
        if(fdata[iter].startswith("#P0")):
          specInfo.append(dict())
          specInfo[specIdx][specInfoNames[0]] = curPol
        pdata = []
        [ pdata.append(float(x)) for x in fdata[iter].split()[1:] ]
        for j in range(len(pdata)):
          specInfo[specIdx][specInfoNames[len(specInfo[specIdx])]] = pdata[j]
        #print(specInfo[specIdx])

      if(fdata[iter].startswith("#L")):
        mnames = []
        [ mnames.append(x) for x in fdata[iter].split()[1:] ]
        #now read scan, until blank line or "#C" line
        tdata = []
        #print(fdata[iter])
        #print(mnames)
        iter = iter + 1
        if(iter < len(fdata)):
          curPol = GetPol(fdata[iter],curPol)
        
        while(iter < len(fdata) and len(fdata[iter]) > 3 and ((not fdata[iter].startswith("#")) or ("injection" in fdata[iter]) or ("communication" in fdata[iter]))):
          if(("injection" not in fdata[iter]) and ("communication" not in fdata[iter])):
            tlist = fdata[iter].split()
            tlist = [float(x) for x in tlist]
            tdata.append(tlist)
          iter = iter+1
          if(iter < len(fdata)):
            curPol = GetPol(fdata[iter],curPol)
        tdata = np.array(tdata)

        #now map the keys to the data lists and store in specData[scanIdx]
        specData.append(dict())
        if(len(tdata) > 0): #otherwise scan was aborted before starting
          for i in range(len(mnames)):
            specData[specIdx][mnames[i]] = tdata[:,i]
        else:
          for i in range(len(mnames)):
            specData[specIdx][mnames[i]] = [0]
        #if(iter - iterold > 1):
        #  iter = iter - 1   
          
      iter = iter+1
      
      
  return specInfo,specData



def ReadSDDFile(fname):
  sddData = []
  specNum = 0
  specIdx = 0
  with open(fname) as f:
    fdata = f.readlines()
    iter = 0
    while(iter < len(fdata)):
        
      if(fdata[iter].startswith("#S")):
        #print(fdata[iter])
        specNum = int(fdata[iter].split()[1])
        specIdx = specNum - 1
        
        
      if(fdata[iter].startswith("#@CALIB")):
        sddData.append(dict())
        tdata = []
        iter = iter + 1
        while(iter < len(fdata) and len(fdata[iter]) > 3 and ((not fdata[iter].startswith("#")) or ("injection" in fdata[iter]) or ("communication" in fdata[iter]))):
          if(("injection" not in fdata[iter]) and ("communication" not in fdata[iter])):
            tdata.append(float(fdata[iter]))
            iter = iter + 1
        sddData[specIdx]["ESDD"] = np.array(tdata)
        

      if(fdata[iter].startswith("#@MCA")):
        tdata = []
        iter = iter + 1
        while(iter < len(fdata) and len(fdata[iter]) > 3 and ((not fdata[iter].startswith("#")) or ("injection" in fdata[iter]) or ("communication" in fdata[iter]))):
          if(("injection" not in fdata[iter]) and ("communication" not in fdata[iter])):
            tlist = fdata[iter].split()
            tlist = [float(x) for x in tlist]
            tdata.append(tlist)
          iter = iter+1
        sddData[specIdx]["Data"] = np.transpose(np.array(tdata))
          
      iter = iter+1
      
      
  return sddData
      



def ReadXESFile(fname):
  xesData = []
  specNum = 0
  specIdx = 0
  with open(fname) as f:
    fdata = f.readlines()
    iter = 0
    while(iter < len(fdata)):
        
      if(fdata[iter].startswith("#S")):
        #print(fdata[iter])
        specNum = int(fdata[iter].split()[1])
        specIdx = specNum - 1
        
        
      if(fdata[iter].startswith("#C MCP Energy Scale")):
        xesData.append(dict())
        tdata = []
        iter = iter + 1
        while(iter < len(fdata) and len(fdata[iter]) > 3 and ((not fdata[iter].startswith("#")) or ("injection" in fdata[iter]) or ("communication" in fdata[iter]))):
          if(("injection" not in fdata[iter]) and ("communication" not in fdata[iter])):
            tdata.append(float(fdata[iter]))
            iter = iter + 1
        xesData[specIdx]["EXES"] = np.array(tdata)
        

      if(fdata[iter].startswith("#@MCA")):
        tdata = []
        iter = iter + 1
        while(iter < len(fdata) and len(fdata[iter]) > 3 and ((not fdata[iter].startswith("#")) or ("injection" in fdata[iter]) or ("communication" in fdata[iter]))):
          if(("injection" not in fdata[iter]) and ("communication" not in fdata[iter])):
            tlist = fdata[iter].split()
            tlist = [float(x) for x in tlist]
            tdata.append(tlist)
          iter = iter+1
        xesData[specIdx]["Data"] = np.transpose(np.array(tdata))
          
      iter = iter+1
      
      
  return xesData
      
            

#
#
#       RXR Functions
#
#
#            
def save_obj(obj, name ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open( name + '.pkl', 'rb') as f:
        return pickle.load(f)

def ReadHenkeToDict():
  path = "henke/"
  savepath = ""
  henke = dict()
  for files in glob.glob(path + "*.nff"):
    p1 = os.path.splitext(files)[0].replace("henke\\","")
    p1 = p1.title() #convert first letter to uppercase
    henke[p1] = np.loadtxt(files,skiprows=1)
    
    henke[p1][np.where(henke[p1][:,1] < -8000),1] = 0
    
    #inData[i][m] = inData[i][m] / np.average(inData[i][m][np.where((inData[i]["MonoEngy"] > elow) & (inData[i]["MonoEngy"] < ehigh) )])
    print(p1)
  save_obj(henke,"HenkePython")
   
   
def WriteReMagX(fname,AScans,AInfo,EScans,EInfo,header):
  
  file = open(fname, "w")
  
  startstr = """ # Complete Status of ReMagX
# File Version
fileversion = 1

# layer configuration
layerversion = 3
# layer 0
 layer = 0
  layer_d = 0
  layer_db = SrTiO3
  layer_material = SrTiO3
  layer_density = 5.12
  layer_delta = 0.00496638
  layer_beta = 0.00166528
  layer_ddelta = 0
  layer_dbeta = 0
  layer_maggamma = 90
  layer_magphi = 90
  layer_sigma = 2
  layer_electron_escape_depth = 0

# layer 1
 layer = 1
  layer_d = 40
  layer_db = LaMnO3
  layer_material = LaMnO3
  layer_density = 6.5
  layer_delta = 0.00423605
  layer_beta = 0.000617009
  layer_ddelta = 0
  layer_dbeta = 0
  layer_maggamma = 90
  layer_magphi = 90
  layer_sigma = 2
  layer_electron_escape_depth = 0

# layer 2
 layer = 2
  layer_d = 40
  layer_db = LaFeO3
  layer_material = LaFeO3
  layer_density = 6.5
  layer_delta = 0.00423605
  layer_beta = 0.000617009
  layer_ddelta = 0
  layer_dbeta = 0
  layer_maggamma = 90
  layer_magphi = 90
  layer_sigma = 2
  layer_electron_escape_depth = 0
  
# layer 3
 layer = 3
  layer_d = 0
  layer_db = vacuum
  layer_material = 
  layer_density = 0
  layer_delta = 0
  layer_beta = 0
  layer_ddelta = 0
  layer_dbeta = 0
  layer_maggamma = 90
  layer_magphi = 90
  layer_sigma = 0
  layer_electron_escape_depth = 0


# element configuration
elementlayerversion = 1

# magnetic configuration with gauss
magneticversion = 4
# gauss function 1
 gauss_function = 1
 gauss_delta = 0
 gauss_beta = 0
 gauss_variance= 1
 gauss_z = 0
 gauss_usespline = 0
 gauss_acoupling = -1
 gauss_varcoupling = -1
 gauss_bindtolayer = -1
 gauss_splinepoints= 0


# fit variables
fitvariableversion = 1


# mode (material or element)
mode = 0

# numerical precision
precision = 1

# disable script
disablescript = 0

# activate parallelization
parallelization = 0

# energy (saved for compatibility)
energy = 385


# angular resolution
resolution = 0

# energy resolution
energyresolution = 0

# min qz to calculate
qzmin = 0

# max qz to calculate
qzmax = 0.25

# number of points to show
numberofpoints = 300

# averaging value
averaging = 0

# temperature for fitting
fittemp = 1

# delta temperature for fitting
fitdtemp = -0.01

# minimum qz for fitting
fitqzmin = 0

# maximum qz for fitting
fitqzmax = 9999

# minimum engergy for fitting
fitenergymin = 0

# maximum energy for fitting
fitenergymax = 9999

# fast calculation of resolution
fastresolution = 1

# use measurement data for determine qz
usemeasurementqz = 1

# show qz(false) or angle (true)
showangle = 0

# fit algorithm
fitalgorithm = -1

# fit error function
fiterrorfunction = 0

# fit resolution
fitresolution = 0

# gausscoupling algorithm
gausscoupling = 1

# multi slicing mode (0=roughness approx. ...)
layersegmentationmode = 0

# algorithm (0=Parratt. ...)
algorithm = 0

# scattering type (0=reflection, 1=transmission)
scattering = 0

# fit measurement multiplicator
fitmeasmult = 0

# fit qz shift
fitdeltaqz = 0

# fit energy shift
fitenergyshift = 0

# energy scan qz
engscanqz = 0.1

# energy scan angle
engscanangle = 10

# energy scan min
engscanmin = 300

# energy scan max
engscanmax = 30000

# polarization of sigma part and pi part in complex numbers
 ray1sigmare = 1
 ray1sigmaim = 0
 ray1pire = 0
 ray1piim = 0
 ray2sigmare = 0
 ray2sigmaim = 0
 ray2pire = 0
 ray2piim = 0

# polarization analyzer
 polanalyzer = 0
# minimal layer thickness
minlayerthickness = 2

# maximal layer thickness
maxlayerthickness = 3

# gradient scale
gradientscale = 1

# segmentation error
segmentationerror = 1e-06

# segmentation model
segmentationmodel = 0

# use magnetization
magnetizationmode = 0

# magnetization direction
maggamma = 90

magphi = 90

# random number seed
evo_random = 2080746695

# population
evo_population = 50

# percent elite of population
evo_elite = 10

# percent of parents
evo_parents = 40

# mutation range in real numbers in percent
evo_range = 10

# probability of mixing
evo_mixing = 33

# probability of mutation
evo_probability = 25

# take elite as parents
evo_takeeliteasparents = 1

# restrict fit parameters
evo_evorestrictfitparameters = 1

# optical constant handling
ocautomatic = 1

# multilayer script
enablemultilayer = 0
multilayerscript = 

# chantler database
usechantlerdatabase = 1

# asymmetry definition
asymmetrydef = 0

# script parameters
# script before calculation

# default script

# active data set index
activedatasetindex = 20
# measurement data
measurementdataversion = 2
"""
  
  if(header == ""):
    file.write(startstr)
  else:
    file.write(header)
    
  dsNum = 1
  for i in range(len(AScans)):
    file.write("datasetnumber = %d \n" % dsNum)
    name = str(AInfo[i][0]) + "_A_" + AInfo[i][3] + "_" + AInfo[i][1]
    file.write("datasettitle = %s \n" % name)
    file.write("datasetenergy = %s \n" % AInfo[i][3])
    file.write("datasetresolution = 0 \n")
    file.write("datasetmultiplicator = 1 \n")
    file.write("datasettranslator = 0 \n")
    file.write("datasetqztranslation = 0 \n")
    file.write("datasetenergyshift = 0 \n")
    if(AInfo[i][1] == "S"):
      file.write("datasetray1sigmare = 1 \n")
      file.write("datasetray1sigmaim = 0 \n")
      file.write("datasetray1pire = 0 \n")
      file.write("datasetray1piim = 0 \n")
    elif(AInfo[i][1] == "P"):
      file.write("datasetray1sigmare = 0 \n")
      file.write("datasetray1sigmaim = 0 \n")
      file.write("datasetray1pire = 1 \n")
      file.write("datasetray1piim = 0 \n")
    elif(AInfo[i][1] == "L"):
      file.write("datasetray1sigmare = 0.707107 \n")
      file.write("datasetray1sigmaim = 0 \n")
      file.write("datasetray1pire = 0 \n")
      file.write("datasetray1piim = 0.707107 \n")
    elif(AInfo[i][1] == "R"):
      file.write("datasetray1sigmare = 0.707107 \n")
      file.write("datasetray1sigmaim = 0 \n")
      file.write("datasetray1pire = 0 \n")
      file.write("datasetray1piim = -0.707107 \n")
    
    file.write("datasetray2sigmare = 1 \n")
    file.write("datasetray2sigmaim = 0 \n")
    file.write("datasetray2pire = 0 \n")
    file.write("datasetray2piim = 0 \n")

    file.write("dataseterrorweight = 1 \n")
    file.write("datasetactive = 0 \n")
    file.write("datasetcomment =  \n")
    file.write("datasetfile =  \n")
    file.write("datasetpoints = %d \n" % len(AScans[i][:,0]))
    for j in range(len(AScans[i][:,0])):
      file.write("dataset_qz = %f \n" % AScans[i][j][2])
      file.write("dataset_R0 = %e \n" % AScans[i][j][3])
      #file.write("dataset_eng = %f \n" % AScans[i][j][0])
    file.write("\n\n")
    dsNum = dsNum + 1
    
    #write asymmetry if possible
    if i>0:
      if(AInfo[i-1][3] == AInfo[i][3]):
        file.write("datasetnumber = %d \n" % dsNum)
        name = str(AInfo[i-1][0]) + "-" + str(AInfo[i][0]) + "_A_" + AInfo[i][3] + "_" + AInfo[i-1][1] + "-" + AInfo[i][1]  + "_Asymm"
        file.write("datasettitle = %s \n" % name)
        file.write("datasetenergy = %s \n" % AInfo[i][3])
        file.write("datasetresolution = 0 \n")
        file.write("datasetmultiplicator = 1 \n")
        file.write("datasettranslator = 0 \n")
        file.write("datasetqztranslation = 0 \n")
        file.write("datasetenergyshift = 0 \n")
        if(AInfo[i-1][1] == "S"):
          file.write("datasetray1sigmare = 1 \n")
          file.write("datasetray1sigmaim = 0 \n")
          file.write("datasetray1pire = 0 \n")
          file.write("datasetray1piim = 0 \n")
        elif(AInfo[i-1][1] == "P"):
          file.write("datasetray1sigmare = 0 \n")
          file.write("datasetray1sigmaim = 0 \n")
          file.write("datasetray1pire = 1 \n")
          file.write("datasetray1piim = 0 \n")
        elif(AInfo[i-1][1] == "L"):
          file.write("datasetray1sigmare = 0.707107 \n")
          file.write("datasetray1sigmaim = 0 \n")
          file.write("datasetray1pire = 0 \n")
          file.write("datasetray1piim = 0.707107 \n")
        elif(AInfo[i-1][1] == "R"):
          file.write("datasetray1sigmare = 0.707107 \n")
          file.write("datasetray1sigmaim = 0 \n")
          file.write("datasetray1pire = 0 \n")
          file.write("datasetray1piim = -0.707107 \n")

        if(AInfo[i][1] == "S"):
          file.write("datasetray2sigmare = 1 \n")
          file.write("datasetray2sigmaim = 0 \n")
          file.write("datasetray2pire = 0 \n")
          file.write("datasetray2piim = 0 \n")
        elif(AInfo[i][1] == "P"):
          file.write("datasetray2sigmare = 0 \n")
          file.write("datasetray2sigmaim = 0 \n")
          file.write("datasetray2pire = 1 \n")
          file.write("datasetray2piim = 0 \n")
        elif(AInfo[i][1] == "L"):
          file.write("datasetray2sigmare = 0.707107 \n")
          file.write("datasetray2sigmaim = 0 \n")
          file.write("datasetray2pire = 0 \n")
          file.write("datasetray2piim = 0.707107 \n")
        elif(AInfo[i][1] == "R"):
          file.write("datasetray2sigmare = 0.707107 \n")
          file.write("datasetray2sigmaim = 0 \n")
          file.write("datasetray2pire = 0 \n")
          file.write("datasetray2piim = -0.707107 \n")



        file.write("dataseterrorweight = 1 \n")
        file.write("datasetactive = 0 \n")
        file.write("datasetcomment =  \n")
        file.write("datasetfile =  \n")
        file.write("datasetpoints = %d \n" % len(AScans[i][:,0]))
        for j in range(len(AScans[i][:,0])):
          file.write("dataset_qz = %f \n" % AScans[i][j][2])
          #print(AScans[i-1][j][3]+AScans[i][j][3])
          file.write("dataset_A = %e \n" % ((AScans[i-1][j][3]-AScans[i][j][3])/(AScans[i-1][j][3]+AScans[i][j][3])))
          #file.write("dataset_eng = %f \n" % AScans[i][j][0])
        file.write("\n\n")
        dsNum = dsNum + 1
    
    
  for i in range(len(EScans)):
    file.write("datasetnumber = %d \n" % dsNum)
    name = str(EInfo[i][0]) + "_E" + str(round(float(EInfo[i][3]),2)) + "_Th" + str(round(float(EInfo[i][4]),2)) + "_" + EInfo[i][1]
    file.write("datasettitle = %s \n" % name)
    file.write("datasetenergy = %s \n" % EInfo[i][3])
    file.write("datasetresolution = 0 \n")
    file.write("datasetmultiplicator = 1 \n")
    file.write("datasettranslator = 0 \n")
    file.write("datasetqztranslation = 0 \n")
    file.write("datasetenergyshift = 0 \n")
    if(EInfo[i][1] == "S"):
      file.write("datasetray1sigmare = 1 \n")
      file.write("datasetray1sigmaim = 0 \n")
      file.write("datasetray1pire = 0 \n")
      file.write("datasetray1piim = 0 \n")
    elif(EInfo[i][1] == "P"):
      file.write("datasetray1sigmare = 0 \n")
      file.write("datasetray1sigmaim = 0 \n")
      file.write("datasetray1pire = 1 \n")
      file.write("datasetray1piim = 0 \n")
    elif(EInfo[i][1] == "L"):
      file.write("datasetray1sigmare = 0.707107 \n")
      file.write("datasetray1sigmaim = 0 \n")
      file.write("datasetray1pire = 0 \n")
      file.write("datasetray1piim = 0.707107 \n")
    elif(EInfo[i][1] == "R"):
      file.write("datasetray1sigmare = 0.707107 \n")
      file.write("datasetray1sigmaim = 0 \n")
      file.write("datasetray1pire = 0 \n")
      file.write("datasetray1piim = -0.707107 \n")
      
    
    file.write("datasetray2sigmare = 1 \n")
    file.write("datasetray2sigmaim = 0 \n")
    file.write("datasetray2pire = 0 \n")
    file.write("datasetray2piim = 0 \n")

    file.write("dataseterrorweight = 1 \n")
    file.write("datasetactive = 0 \n")
    file.write("datasetcomment =  \n")
    file.write("datasetfile =  \n")
    file.write("datasetpoints = %d \n" % len(EScans[i][:,0]))
    for j in range(len(EScans[i][:,0])):
      file.write("dataset_qz = %f \n" % EScans[i][j][2])
      file.write("dataset_R0 = %e \n" % EScans[i][j][3])
      file.write("dataset_eng = %f \n" % EScans[i][j][0])
    file.write("\n\n")
    dsNum = dsNum + 1    
    
    
    
    #write asymmetry if possible
    if i>0:
      if(abs(float(EInfo[i-1][3]) - float(EInfo[i][3])) < 0.015 and abs(float(EInfo[i-1][4]) - float(EInfo[i][4])) < 0.1 ):
        file.write("datasetnumber = %d \n" % dsNum)
        name = str(EInfo[i-1][0]) + "-" + str(EInfo[i][0]) + "_E" + str(round(float(EInfo[i][3]),2)) + "_Th" + str(round(float(EInfo[i][4]),2)) + "_" + EInfo[i-1][1] + "-" + EInfo[i][1] + "_Asymm"
        file.write("datasettitle = %s \n" % name)
        file.write("datasetenergy = %s \n" % EInfo[i][3])
        file.write("datasetresolution = 0 \n")
        file.write("datasetmultiplicator = 1 \n")
        file.write("datasettranslator = 0 \n")
        file.write("datasetqztranslation = 0 \n")
        file.write("datasetenergyshift = 0 \n")
        if(EInfo[i-1][1] == "S"):
          file.write("datasetray1sigmare = 1 \n")
          file.write("datasetray1sigmaim = 0 \n")
          file.write("datasetray1pire = 0 \n")
          file.write("datasetray1piim = 0 \n")
        elif(EInfo[i-1][1] == "P"):
          file.write("datasetray1sigmare = 0 \n")
          file.write("datasetray1sigmaim = 0 \n")
          file.write("datasetray1pire = 1 \n")
          file.write("datasetray1piim = 0 \n")
        elif(EInfo[i-1][1] == "L"):
          file.write("datasetray1sigmare = 0.707107 \n")
          file.write("datasetray1sigmaim = 0 \n")
          file.write("datasetray1pire = 0 \n")
          file.write("datasetray1piim = 0.707107 \n")
        elif(EInfo[i-1][1] == "R"):
          file.write("datasetray1sigmare = 0.707107 \n")
          file.write("datasetray1sigmaim = 0 \n")
          file.write("datasetray1pire = 0 \n")
          file.write("datasetray1piim = -0.707107 \n")

        if(EInfo[i][1] == "S"):
          file.write("datasetray2sigmare = 1 \n")
          file.write("datasetray2sigmaim = 0 \n")
          file.write("datasetray2pire = 0 \n")
          file.write("datasetray2piim = 0 \n")
        elif(EInfo[i][1] == "P"):
          file.write("datasetray2sigmare = 0 \n")
          file.write("datasetray2sigmaim = 0 \n")
          file.write("datasetray2pire = 1 \n")
          file.write("datasetray2piim = 0 \n")
        elif(EInfo[i][1] == "L"):
          file.write("datasetray2sigmare = 0.707107 \n")
          file.write("datasetray2sigmaim = 0 \n")
          file.write("datasetray2pire = 0 \n")
          file.write("datasetray2piim = 0.707107 \n")
        elif(EInfo[i][1] == "R"):
          file.write("datasetray2sigmare = 0.707107 \n")
          file.write("datasetray2sigmaim = 0 \n")
          file.write("datasetray2pire = 0 \n")
          file.write("datasetray2piim = -0.707107 \n")



        file.write("dataseterrorweight = 1 \n")
        file.write("datasetactive = 0 \n")
        file.write("datasetcomment =  \n")
        file.write("datasetfile =  \n")
        file.write("datasetpoints = %d \n" % len(EScans[i][:,0]))
        for j in range(len(EScans[i][:,0])):
          file.write("dataset_qz = %f \n" % EScans[i][j][2])
          file.write("dataset_A = %e \n" % ((EScans[i-1][j][3]-EScans[i][j][3])/(EScans[i-1][j][3]+EScans[i][j][3])))
          file.write("dataset_eng = %f \n" % EScans[i][j][0])
        
        file.write("\n\n")
        dsNum = dsNum + 1
    
    
    
    
    
  file.close()
  

      
      
#
#
#
#      General data processing
#
#
#
#
      
def NormByMonitor(inData,scanNums,sigList,monitor):

  for i in range(len(inData)):
    if(i+1 in scanNums):
      print(monitor, np.mean(inData[i][monitor]))
      for x in sigList:
        #print(i+1)
        inData[i][x] = inData[i][x] / inData[i][monitor]
  return inData
  
def GetDirBeamCorrection(corrInfo):
  cInfo,cData = ReadSpecFile(corrInfo[0])
  corr = dict()
  for i in range(1,len(corrInfo)):
    cData[corrInfo[i][0]-1]["PicoAm3"] = cData[corrInfo[i][0]-1]["PicoAm3"] / (cData[corrInfo[i][0]-1]["I0_BD3"] )
    tcorr = np.transpose(np.array([cData[corrInfo[i][0]-1]["MonoEngy"], cData[corrInfo[i][0]-1]["PicoAm3"]]))
    
    if(corrInfo[i][1] in corr):
      corr[corrInfo[i][1]] = np.vstack((corr[corrInfo[i][1]],tcorr))
    else:
      corr[corrInfo[i][1]] = tcorr     
  return corr
    
  
def GetECalList(inPairs):
  #list of pairs of [measE,calE]
  #automatically gives the number of polynomial coefficients to shift measE to calE
  EMat = np.zeros((len(inPairs),len(inPairs)))
  coeffs = np.zeros((len(inPairs)))
  ECal = np.zeros((len(inPairs)))
  
  for i in range(len(inPairs)):
    ECal[i] = inPairs[i][1]
    for j in range(len(inPairs)):
      EMat[i][j] = inPairs[i][0]**(len(inPairs)-j-1)
      
  coeffs = np.matmul(np.linalg.inv(EMat),ECal)

  return np.flip(coeffs,axis=0)
  
#Calibrate energy using arbitrary number of polynomial coefficients
def ECalPoly(inData,scanNums,coeffs):
  for i in range(len(inData)):
    if(i+1 in scanNums and "MonoEngy" in inData[i]):
      enew = inData[i]["MonoEngy"] * 0
      for j in range(len(coeffs)):
        enew = enew + coeffs[j] * inData[i]["MonoEngy"] ** j
      inData[i]["MonoEngy"] = enew
  return inData

#Get calibrated value of a single energy
def ECalPolySingle(Ein,coeffs):
  Enew = 0
  for j in range(len(coeffs)):
    Enew = Enew + coeffs[j] * Ein ** j
  return Enew  
  
#Polarization name convention conversions
def ConvPolVtoS(inPol):
  if(inPol == "LV"):
    return "S"
  elif(inPol == "LH"):
    return "P"
  elif(inPol == "CL"):
    return "L"
  elif(inPol == "CR"):
    return "R"
  else:
    return "U"
    
def ConvPolStoV(inPol):
  if(inPol == "S"):
    return "LV"
  elif(inPol == "P"):
    return "LH"
  elif(inPol == "L"):
    return "CL"
  elif(inPol == "R"):
    return "CR"
  else:
    return "U"
    

#Return RXR data for given energy scans (Fixed Q or Fixed Angle)   
def GetRXREScanData(inData,inInfo,inScans):
  retList = []
  
  for iscan in range(len(inScans[:,0])):
    snum = inScans[iscan,0]

    if(snum-1 <= len(inData)):
      i = snum-1
      tr1 = inScans[np.nonzero(inScans[:,0] == i+1),1][0][0] #inScans[inScans[:,0].index(i+1),1]
      tr2 = inScans[np.nonzero(inScans[:,0] == i+1),2][0][0] #inScans[inScans[:,0].index(i+1),2]
      if(tr2 == 0):
        retList.append(np.transpose(np.array([inData[i]["MonoEngy"][tr1:],0*inData[i]["MonoEngy"][tr1:],0*inData[i]["MonoEngy"][tr1:],inData[i]["PicoAm3"][tr1:]])))
      else:
        retList.append(np.transpose(np.array([inData[i]["MonoEngy"][tr1:-tr2],0*inData[i]["MonoEngy"][tr1:-tr2],0*inData[i]["MonoEngy"][tr1:-tr2],inData[i]["PicoAm3"][tr1:-tr2]])))
      
      retList[len(retList)-1][:,1] = inInfo[i]["Theta"]
      retList[len(retList)-1][:,2] = 0.001013546247 * retList[len(retList)-1][:,0] * math.sin(inInfo[i]["Theta"]*math.pi/180) 

  return retList

#Return RXR info for given scans (Fixed Q or Fixed Angle)
def GetRXREScanInfo(inData,inInfo,inScans):
  retList = []

  for iscan in range(len(inScans[:,0])):
    snum = inScans[iscan,0]

    if(snum-1 <= len(inData)):
      i = snum-1

      retList.append([i+1,ConvPolVtoS(inInfo[i]["Polarization"]),"E",inData[i]["MonoEngy"][0],inInfo[i]["Theta"]])
  return np.array(retList)


#Return RXR data for given energy scans (Fixed Energy)   
def GetRXRAScanData(inData,inInfo,inScans):
  retList = []
  for iscan in range(len(inScans[:,0])):
    snum = inScans[iscan,0]

    if(snum-1 <= len(inData)):
      i = snum-1
      tr1 = inScans[np.nonzero(inScans[:,0] == i+1),1][0][0] #inScans[inScans[:,0].index(i+1),1]
      tr2 = inScans[np.nonzero(inScans[:,0] == i+1),2][0][0] #inScans[inScans[:,0].index(i+1),2]
      if(tr2 == 0):
        retList.append(np.transpose(np.array([0*inData[i]["Theta"][tr1:],inData[i]["Theta"][tr1:],0*inData[i]["Theta"][tr1:],inData[i]["PicoAm3"][tr1:]])))
      else:
        retList.append(np.transpose(np.array([0*inData[i]["Theta"][tr1:-tr2],inData[i]["Theta"][tr1:-tr2],0*inData[i]["Theta"][tr1:-tr2],inData[i]["PicoAm3"][tr1:-tr2]])))
      
      
      retList[len(retList)-1][:,0] = inInfo[i]["MonoEngy"]
      retList[len(retList)-1][:,2] = 0.001013546247 * inInfo[i]["MonoEngy"] * np.sin(math.pi / 180 * retList[len(retList)-1][:,1] ) 

  return retList

#Return RXR info for given scans (Fixed Energy)
def GetRXRAScanInfo(inData,inInfo,inScans):
  retList = []
  for iscan in range(len(inScans[:,0])):
    snum = inScans[iscan,0]

    if(snum-1 <= len(inData)):
      i = snum-1
      #print(i+1)
      retList.append([i+1,ConvPolVtoS(inInfo[i]["Polarization"]),"A",inInfo[i]["MonoEngy"],inData[i]["Theta"][0]])
  return np.array(retList)

#Numerical integral of a gaussian?  
def GaussInt(fwhm,x0,x1,N):
  c = fwhm/2.35482004503
  c2 = 2*c*c
  sum = 0
  for i in range(N):
    xc = x0 + i/(N-1)*(x1-x0)
    sum = sum + math.exp(-xc*xc/c2)
  return sum / c / math.sqrt(2*math.pi) * 1 /(N-1)*(x1-x0)

def GaussianIntegral(fwhm,xc,x1,x2):
  #100 points per fwhm
  c = fwhm/2/np.sqrt(2*np.log(2))
  c2 = 2*c**2
  a = 1/c/np.sqrt(2*np.pi)
  N = int(np.abs(x2-x1)/fwhm*100)
  evals = np.linspace(x1,x2,N)
  gau = a*np.exp(-(evals-xc)**2/c2)
  return np.sum(gau[:-1])*(evals[1]-evals[0])
  #print(fwhm,xc,x1,x2)
  #print(np.sum(gau)/N, np.sum(gau[:-1])*(evals[1]-evals[0]))
  #exit()
  

#Apply direct beam correction to set of scans  
def ApplyDirBeam(inScans,inInfo,corrData):
  for i in range(len(inScans)):
    inScans[i][:,3] = inScans[i][:,3] / np.interp(inScans[i][:,0],corrData[ConvPolStoV(inInfo[i][1])][:,0],corrData[ConvPolStoV(inInfo[i][1])][:,1])    
  return inScans

#Subtract constant background from a set of RXR scans   
def BackgroundSubRXR(inScans,bgVal):
  for i in range(len(inScans)):
    inScans[i][:,3] = inScans[i][:,3] - bgVal
  return inScans

#Not implemented - take absolute value of scan intensities  
def AbsValData(inScans):
  return inScans
  

def ApplyGeometryCorrection(inScans,geoData):
#ordering of pars is: fwhm,samplesize,detectorsize,detectordistance,relbeampos,samplex,sampley,detectordelta,detectorx

#footprint on sample (i.e. effective fwhm in sample plane) is d=fwhm/sin(theta)
#integrated flux on sample is 2*int(gauss of width d)| from 0 to L/2

#  for i in range(len(inScans)):
#    d = geoData[0] / np.sin(math.pi/180*inScans[i][:,1]) / 2.35482004503
#    eval = geoData[1]/2 / math.sqrt(2) / d
#    inScans[i][:,3] = inScans[i][:,3] / (scipy.special.erf(eval))
#  return inScans

#ordering of pars is: fwhm,samplesize,detectorsize,detectordistance,relbeampos,samplex,sampley,detectordelta,detectorx

#footprint on sample plane (i.e. effective fwhm in sample plane) is d=fwhm/sin(theta)
#x or y offset shifts center of gaussian in sample plane to point p
#integrated flux on sample is int(gauss of width d) from 0 to (width/2-p) and plus integral from 0 to width/2+p

  for i in range(len(inScans)):
    #inScans[i][:,3] = inScans[i][:,3] * GaussInt(geoData[0],0,geoData[1]*np.sin(math.pi/180*inScans[i][:,1])/2,20) / GaussInt() 
    
    d = geoData[0] / np.sin(math.pi/180*inScans[i][:,1]) / 2.35482004503
    
    #GaussianIntegral(1,0,-1/2.35,1/2.35)
    
    p = np.abs(-geoData[5]/np.tan(math.pi/180*inScans[i][:,1]) + geoData[6])

    for j in range(len(d)):
      inScans[i][j,3] = inScans[i][j,3] / GaussianIntegral(d[j],p[j],-geoData[1]/2,geoData[1]/2) #(scipy.special.erf(eval1)+scipy.special.erf(eval2))
    
    
    #eval1 = np.abs(geoData[1]/2-p) / math.sqrt(2) / d
    #eval2 = np.abs(geoData[1]/2+p) / math.sqrt(2) / d
    #if(inScans[i][2,1] < 10):
    #  print(inScans[i][2,1], eval1[2], eval2[2], scipy.special.erf(eval1[2]), scipy.special.erf(eval2[2]))
    #  exit()
    #print(inScans[i][:,1])
    #print((scipy.special.erf(eval)))
    #exit()
    #inScans[i][:,3] = inScans[i][:,3] / GaussianIntegral(d,p,-geoData[1]/2,geoData[1]/2) #(scipy.special.erf(eval1)+scipy.special.erf(eval2))
    #inScans[i][:,3] = inScans[i][:,3] / GaussInt(geoData[0] / (1 - np.sin(math.pi/180*inScans[i][:,1]),0,geoData[2],20))

#with a scaling factor, this gives integral of gaussian (so can calculate the amount of photons hitting sample)



  return inScans  
  

def RemoveBadPoints(inData,scanNums,mons):  
  for i in range(len(inData)):
    if(i+1 in scanNums):
      for m in mons:
        if(m in inData[i]):
          for k in range(0,len(inData[i][m])):
            if(k==0 and len(inData[i][m]) > 2):
              if((inData[i][m][k]/(inData[i][m][k+1]+1e-20) > 2 or inData[i][m][k]/(inData[i][m][k+1]+1e-20) < 0.5) and (inData[i][m][k+1]/(inData[i][m][k+2]+1e-20) > 2 or inData[i][m][k+1]/(inData[i][m][k+2]+1e-20) < 0.5)):   
                inData[i][m][k] = inData[i][m][k+1]
            elif(k==(len(inData[i][m])-1) and len(inData[i][m]) > 2):
              if((inData[i][m][k]/(inData[i][m][k-1]+1e-20) > 2 or inData[i][m][k]/(inData[i][m][k-1]+1e-20) < 0.5) and (inData[i][m][k-1]/(inData[i][m][k-2]+1e-20) > 2 or inData[i][m][k-1]/(inData[i][m][k-2]+1e-20) < 0.5)):   
                inData[i][m][k] = inData[i][m][k-1]
            elif(len(inData[i][m]) > 2):
              if((inData[i][m][k]/(inData[i][m][k-1]+1e-20) > 2 or inData[i][m][k]/(inData[i][m][k-1]+1e-20) < 0.5) and (inData[i][m][k]/(inData[i][m][k+1]+1e-20) > 2 or inData[i][m][k]/(inData[i][m][k+1]+1e-20) < 0.5)):   
                #print(k, len(inData[i][m]))
                inData[i][m][k] = (inData[i][m][k-1] + inData[i][m][k+1])/2
          
  return inData

def XASSubBG(inData,monE,mons,scanNums,elow,ehigh):
  for scan in scanNums:
    i = scan-1
    idx = (inData[i][monE] > elow) & (inData[i][monE] < ehigh)

    xmean = np.average(inData[i][monE][idx])
    for mon in mons:
      ymean = np.average(inData[i][mon][idx])
      m = np.sum((inData[i][monE][idx]-xmean) * (inData[i][mon][idx]-ymean))/ np.sum((inData[i][monE][idx]-xmean)**2)
      b = ymean - m*xmean
      inData[i][mon] = inData[i][mon] - (m*inData[i][monE]+b)
  return 

  
  
def XASSubBGCons(inData,mons,scanNums,elow,ehigh):
  
  #for each scan, fit a linear line to the region, and subtract that line from spectrum
  #currently just substracting a slope=0 line
  for i in range(len(inData)):
    if(i+1 in scanNums):
      for m in mons:
        if(m in inData[i]):
          inData[i][m] = inData[i][m] - np.average(inData[i][m][np.where((inData[i]["MonoEngy"] > elow) & (inData[i]["MonoEngy"] < ehigh) )])

  return inData

#This is for ipfy  
def InvertMonitor(inData,scanNums,mons):
  
  for i in range(len(inData)):
    if(i+1 in scanNums):
      for m in mons:
        if(m in inData[i]):
          inData[i][m] = 1/inData[i][m]

  return inData

  
def XASNormRegion(inData,mons,scanNums,elow,ehigh):

  for i in range(len(inData)):
    if(i+1 in scanNums):
      for m in mons:
        if(m in inData[i]):
          inData[i][m] = inData[i][m] / np.average(inData[i][m][np.where((inData[i]["MonoEngy"] > elow) & (inData[i]["MonoEngy"] < ehigh) )])
     
  return inData
  
def ProcessXAS(fname,scanNums,mon,ECal,eVals,sum):

  sInfo, sData = ReadSpecFile(fname)
  sData = NormByMonitor(sData,scanNums,[mon],"I0_BD3")
  sData = ECalPoly(sData,scanNums,[0,1,0.00])
  
  sData = XASSubBG(sData,[mon],scanNums,eVals[0],eVals[1])
  sData = XASNormRegion(sData,[mon],scanNums,eVals[2],eVals[3])
  
  #PlotSpecData(sData,scanNums,"MonoEngy","MCP_REIXS")

  myener = sData[scanNums[0]-1]["MonoEngy"]
  myspec = sData[scanNums[0]-1][mon] 
  if(sum):
    for i in scanNums[1:]:
      myspec = myspec + sData[i-1][mon]
  else:
    for i in scanNums[1:]:
      myspec = np.vstack((myspec,sData[i-1][mon]))
    myspec = np.transpose(myspec)

  return myener,myspec
  

      
def ProcessRXR(fname,sScans,ECal,Geo,Corr,sType):
  #read in spec file
  
  scanNums = sScans[:,0]
  sInfo, sData = ReadSpecFile(fname)
  #sData = NormByMonitor(sData,scanNums,["I0_BD3_Amp"],"Seconds")
  sData = NormByMonitor(sData,scanNums,["TEY_REIXS","MCP_REIXS","PicoAm3"],"I0_BD3")
  sData = ECalPoly(sData,scanNums,ECal)
  corrData = GetDirBeamCorrection(Corr)
  
  #need to calibrate energy in sInfo
  for i in range(len(sInfo)):
    sInfo[i]["MonoEngy"] = round(ECalPolySingle(sInfo[i]["MonoEngy"],ECal),2)
  
  #SaveScans("MgO-GaN22",sData,sInfo,sScans,["MonoEngy","TEY_REIXS","MCP_REIXS"])

  
  
  
  if(sType == "E"):
    scansInfo = GetRXREScanInfo(sData,sInfo,sScans) #snum, pol(s/p/L/R), E/A, EStart,ThStart
    scans = GetRXREScanData(sData,sInfo,sScans) #Ener,Th,Qz,R
  elif(sType == "A" or sType == "Q"):
    scansInfo = GetRXRAScanInfo(sData,sInfo,sScans) #snum, pol(s/p/L/R), E/A, EStart,ThStart
    scans = GetRXRAScanData(sData,sInfo,sScans) #Ener,Th,Qz,R
  else:
    print("Unknown scan type ", sType)
    exit()
  

  # for i in range(len(scansInfo)):
    # print(scansInfo[i][3])
    # exit()
    # scansInfo[i][3] = ECalPolySingle(scansInfo[i][3],ECal)
  
  scans = ApplyDirBeam(scans,scansInfo,corrData)

  #np.savetxt("temp.dat",scans[0])





  #
  #
  #
  #   TEMP: ADDED JAN 6 2019
  #
  #
  #
  for i in range(len(scans)):
    scans[i][:,3] = np.abs(scans[i][:,3])
  
    #find minimum
  mins = np.zeros((len(scans)))
  for i in range(len(scans)):
    mins[i] = np.amin(scans[i][:,3])
  

  #print(mins)
  #print(np.argpartition(mins,int(len(mins)/10)))
  #print(mins[np.argpartition(mins,int(len(mins)/10))])
  
  if(np.amin(mins) < 0):
    mins = mins[mins < 0]

  
  min = np.mean(mins)
  min = np.amin(mins)
  print("Mins: \n", mins)
  
  print("MIN: ", min)
  
  
  
  #if(min < 1e-7):
  #  scans = BackgroundSubRXR(scans, min - 1e-7)
  
  scans = ApplyGeometryCorrection(scans,Geo)
  
  #print(sInfo[17])

#  plt.figure()
  #plt.plot(sData[50]["MonoEngy"],sData[50]["MCP_REIXS"])#/sData[50]["I0_BD3_Amp"])
#  plt.plot(scans[0][:,0],scans[0][:,3])
#  plt.show()
  
  return scans,scansInfo



def SaveScans(fpre,inData,inInfo,inScans,monitors):

  for i in range(len(inData)):
    if(i+1 in inScans[:,0]):
      file = open(fpre + "_S" + format(i+1,'03') + ".dat", "w")
      for j in monitors:
        file.write("%s " % j)
      file.write("\n")
      
      for j in range(len(inData[i][list(inData[i].keys())[0]][:])):
        for k in monitors:
          if(k in inData[i]):
            file.write("%e " % inData[i][k][j])
        file.write("\n")
        
      file.close()



def WriteQUADData():
  print("hi")

def GetHenke(el):

  return ff
  
  
def MergeOffRes(ener,XAS,element,eVals):#,eL1,eH1,eL2,eH2):
  #get henke data
  henke = load_obj("HenkePython")
  orSpec = henke[element]
  
  orf2 = np.interp(ener,orSpec[:,0],orSpec[:,2])
  y2pre  = np.ndarray.flatten(orf2[np.where((ener > eVals[0]) & (ener < eVals[1]))])
  y2post = np.ndarray.flatten(orf2[np.where((ener > eVals[2]) & (ener < eVals[3]))])
  
  #spec = A*spec + B + C*energy
  #minimize error in pre-post with offres
  
  #fit to preedge
  
  #scale to postedge
  #insert exp into henke
  #return
  
  
  
  mydat = np.transpose(np.vstack((ener,XAS)))
  #evals = [eL1,eH1,eL2,eH2]
  
  
  
  #three pars:
  #Y = AY + B + CX (scale, vertical shift, straight line subtraction)
  
  #rough order of magnitude initialization
  pars = [1/(XAS[-1]+1e-10),XAS[0]/(XAS[-1]+1e-10), 0.1*(XAS[-1]-XAS[0])/(ener[-1]-ener[0])/(XAS[-1]+1e-10)]
  #print(pars)
  
  def XLDFun(inx,inpar,evals):
    eL1 = evals[0]
    eH1 = evals[1]
    eL2 = evals[2]
    eH2 = evals[3]
    #print(inpar)
    
    curpar = np.copy(inpar)
    curpar[:,1] = inx[0]*curpar[:,1] + inx[1] + inx[2] * curpar[:,0]
    
    y1pre  = np.ndarray.flatten(curpar[np.where((curpar[:,0] > eL1) & (curpar[:,0] < eH1)),1])
    y1post = np.ndarray.flatten(curpar[np.where((curpar[:,0] > eL2) & (curpar[:,0] < eH2)),1])
    
    
    #resid = [xld in pre/post] | [pre1] | [pre2] | [post1 - norm(post1)] | [post2 - norm(post2)]
    resid = y1pre - y2pre
    resid = np.hstack((resid,y1post-y2post))
    #print(resid)
    return resid
    

  res_lsq = least_squares(XLDFun,pars,args=(mydat,eVals),bounds=([0,-np.inf,-np.inf],np.inf))
  
  #print(res_lsq)
  print("In MergeOffres: ",res_lsq.message)
  
  
  XAS = XAS * res_lsq.x[0] + res_lsq.x[1] + res_lsq.x[2] * ener
  
  #now merge
  #print(np.shape(np.ndarray.flatten(orSpec[np.where(orSpec[:,0] < eL1),0])),np.shape(ener))
  eL1 = eVals[0]; eH1 = eVals[1]; eL2 = eVals[2]; eH2 = eVals[3];
  
  XAS = np.hstack((np.ndarray.flatten(orSpec[np.where(orSpec[:,0] < eL1),2]),XAS[np.where((ener >= eL1)&(ener <= eH2))]))
  XAS = np.hstack((XAS,np.ndarray.flatten(orSpec[np.where(orSpec[:,0] > eH2),2])))
  ener = np.hstack((np.ndarray.flatten(orSpec[np.where(orSpec[:,0] < eL1),0]),ener[np.where((ener >= eL1)&(ener <= eH2))]))
  ener = np.hstack((ener,np.ndarray.flatten(orSpec[np.where(orSpec[:,0] > eH2),0])))
  
  #now kramers kronig
  oset = Getf1Offset(element)
  f1 = KKf1Fromf2(ener,XAS,oset)
  
  #plt.figure()
  #plt.plot(ener,XAS)
  #plt.plot(ener,f1)
  #plt.show()
  
  return np.transpose(np.vstack((ener,f1,XAS)))
  
  
  
def KKf1Fromf2(e,im,offset):
#The inner loop here is vectorized compared to the slow version
#The algorithm itself is from Sebastian's online (javascript) kramers kronig calculator
#To do the vectorization I had to remove his "dist" condition where he only uses the log expressions near the current energy point
#The speedup from vectorization was much faster than any slowdown due to using the long log expression for all points

  re = np.zeros(len(im)) + offset

  newE = e + 0.00234563
  
  for i in range(1,len(im)-1):
    xm = e[i-1]
    x0 = e[i]
    xp = e[i+1]
    beta = im[i]
    
    m1 = -beta/(xp-x0)
    b1 = beta - m1*x0
    
    m2 = beta / (x0-xm)
    b2 = beta - m2*x0
    
    delta= 0.5 * b1 * np.log(np.fabs((newE-xp)*(newE+xp))/np.fabs((newE-x0)*(newE+x0))) + m1*(xp-x0) + 0.5 * newE * m1 * np.log((newE+x0)/(newE+xp)) + 0.5 * newE * m1 * np.log(np.fabs(newE-xp)/np.fabs(newE-x0))
    delta = delta + 0.5 * b2 * np.log(np.fabs((newE-x0)*(newE+x0))/np.fabs((newE-xm)*(newE+xm))) + m2*(x0-xm) + 0.5 * newE * m2 * np.log((newE+xm)/(newE+x0)) + 0.5 * newE * m2 * np.log(np.fabs(newE-x0)/np.fabs(newE-xm))
    re = re - 2.0/3.14159265359 * delta
      
  return re

def Getf1Offset(element):
  switcher = { "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, 
               "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
               "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
               "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
               "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
               "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103,
               "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111 }
  return switcher.get(element, -1000)
      
  
  
  
  
  
  
  

  
def PlotFixedA(scans):
  plt.figure()
  plt.yscale('log')
  for i in range(len(scans)):
    plt.plot(scans[i][:,2],scans[i][:,3]*(0.1**i))
  plt.show()
 
def PlotSpecData(inData,inScans,monitorx,monitory):
  plt.figure()
  for i in range(len(inData)):
    if(i+1 in inScans):
      if(monitorx in inData[i] and monitory in inData[i]):
        plt.plot(inData[i][monitorx][:],inData[i][monitory][:])

  plt.show()      

  
def SubLinearBackground(x,y,indices):
  #fit a line using
  slope, intercept, r_value, p_value, std_err = stats.linregress(x[indices],y[indices])
  return y - (slope*x + intercept)
  
def SubtractPolyFit(inX,inY,deg,Elow,Ehigh):
  #subtract polynomial of degree deg fitted between regions Elow and Ehigh
  
  coeffs = np.polyfit(inX[np.where((inX > Elow) & (inX < Ehigh))],inY[np.where((inX > Elow) & (inX < Ehigh))],deg)
  #print(coeffs)
  for i in range(len(coeffs)):
    inY = inY - coeffs[i] * inX ** (len(coeffs)-i-1)
  
  return inY
  
#An example of data processing   
def RunCode():  
 
  #print(scipy.special.erf(0),scipy.special.erf(0.5),scipy.special.erf(1),scipy.special.erf(10))
  #exit()


  #this reads henke from subdir and then saves a python henke file
  #ReadHenkeToDict()
  #exit()
  #so from now on only need the python file, which can be loaded using
  #henke = load_obj("HenkePython")


  henke = load_obj("HenkePython")
  plt.figure()
  plt.plot(henke["mn"][:,0],henke["mn"][:,1])
  plt.plot(henke["mn"][:,0],henke["mn"][:,2])
  plt.show()
  exit()

  ener,lvspec = ProcessXAS("MgO-GaN22-8nm.spec",list(range(19,66,2)),[0,1])
  ener,lhspec = ProcessXAS("MgO-GaN22-8nm.spec",list(range(20,67,2)),[0,1])

  plt.figure()
  plt.plot(ener,lvspec)
  plt.plot(ener,lhspec)
  plt.show()

  exit()

  sInfo, sData = ReadSpecFile("MgO-GaN22-8nm.spec")
  plscans = [24, 25]
  PlotSpecData(sData,plscans,"MonoEngy","MCP_REIXS")


  exit()

  EScan,AScan,ECal,Geo,Corr = GetSampleInfo("MgO-GaN22-8nm.spec")
  AsData,AsInfo = ProcessRXR("MgO-GaN22-8nm.spec", AScan,ECal,Geo,Corr,"A")
  EsData,EsInfo = ProcessRXR("MgO-GaN22-8nm.spec", EScan,ECal,Geo,Corr,"E")
    
    
    
  #WriteReMagX("MgO-GaN22-8nm.all",AsData,AsInfo,EsData,EsInfo)
    
    
  PlotFixedA(AsData)
  
#RunCode()  