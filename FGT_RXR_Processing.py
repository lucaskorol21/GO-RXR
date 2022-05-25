
#the line below should be the directory where SpecFileProcessing.py is stored
#import sys
#sys.path.insert(1, 'C:\\Users\\Robert\\Google Drive\\Beamtimes')


from SpecFileProcessing import *


datadir = "./Data/"



def GetReMagXHeader(names,densities,thicknesses):
  header = """# Complete Status of ReMagX
# File Version
fileversion = 1

# layer configuration
layerversion = 3\n"""

  for i in range(len(names)):
    header = header + "# layer " + str(i) + "\n"
    header = header + "layer = " + str(i) + "\n"
    header = header + "layer_d = " + str(thicknesses[i]) + "\n"
    header = header + "layer_db = " + names[i] + "\n"
    header = header + "layer_material = " + names[i] + "\n"
    header = header + "layer_density = " + str(densities[i]) + "\n"
    header = header + "layer_delta = 0.001\n" 
    header = header + "layer_beta = 0.001\n" 
    header = header + "layer_ddelta = 0\n" 
    header = header + "layer_dbeta = 0\n" 
    header = header + "layer_maggamma = 90\n" 
    header = header + "layer_magphi = 90\n" 
    header = header + "layer_sigma = 1.5\n"
    header = header + "layer_electron_escape_depth = 0\n\n"
    
    
    
  header = header + "# layer " + str(len(names)) + "\n"
  header = header + " layer = " + str(len(names)) + "\n"
  
  header = header + """
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

 element = Ge
  el_f1 = 0
  el_f2 = 0
  el_f1m = 0
  el_f2m = 0
   el_layer = 0
    el_d = 0
    el_density = 0.0732792
    el_sigma = 0.1
   el_layer = 1
    el_d = 1.05
    el_density = 0
    el_sigma = 0.1
   el_layer = 2
    el_d = 20
    el_density = 0
    el_sigma = 0.1
   el_layer = 3
    el_d = 0
    el_density = 0
    el_sigma = 0

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
energy = 500


# angular resolution
resolution = 0.0001

# energy resolution
energyresolution = 0

# min qz to calculate
qzmin = 0

# max qz to calculate
qzmax = 0.505

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
showangle = 1

# fit algorithm
fitalgorithm = 0

# fit error function
fiterrorfunction = 0

# fit resolution
fitresolution = 0

# gausscoupling algorithm
gausscoupling = 0

# multi slicing mode (0=roughness approx. ...)
layersegmentationmode = 1

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
engscanqz = 0.1821

# energy scan angle
engscanangle = 80

# energy scan min
engscanmin = 680

# energy scan max
engscanmax = 750

# polarization of sigma part and pi part in complex numbers
 ray1sigmare = 1
 ray1sigmaim = 0
 ray1pire = 0
 ray1piim = 0
 ray2sigmare = 0.707107
 ray2sigmaim = 0
 ray2pire = 0
 ray2piim = -0.707107

# polarization analyzer
 polanalyzer = 0
# minimal layer thickness
minlayerthickness = 0.1

# maximal layer thickness
maxlayerthickness = 2

# gradient scale
gradientscale = 1

# segmentation error
segmentationerror = 1e-06

# segmentation model
segmentationmodel = 1

# use magnetization
magnetizationmode = 0

# magnetization direction
maggamma = 90

magphi = 90

# random number seed
evo_random = 100

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



scriptline =  --params: 
scriptline =  --Fe thickness for all
scriptline =  --Fe sigma for all
scriptline =  --Fe density for FGT layer
scriptline =  --Ge density for substrate and top two layers
scriptline =  --O density for top two layers
scriptline =  --C density for top layer
scriptline =  
scriptline =  for i=1,5 do
scriptline =    t = Fe:getthickness(i)
scriptline =    Ge:setthickness(i,t)
scriptline =    Te:setthickness(i,t)
scriptline =    O:setthickness(i,t)
scriptline =    C:setthickness(i,t)
scriptline =  end
scriptline =  
scriptline =  for i=1,3 do
scriptline =    d = Fe:getdensity(i)
scriptline =    Ge:setdensity(i, d/3)
scriptline =    Te:setdensity(i, 2*d/3)
scriptline =  end
scriptline =  
scriptline =  for i=0,5 do
scriptline =    s = Fe:getsigma(i)
scriptline =    Ge:setsigma(i,s)
scriptline =    Te:setsigma(i,s)
scriptline =    O:setsigma(i,s)
scriptline =    C:setsigma(i,s)
scriptline =  end


# default script

# active data set index
activedatasetindex = 0
# measurement data
measurementdataversion = 2

"""


  return header


def PCRDRXR():

  sInfo,sData = ReadSpecFile("./Data/Mn3Ge-PCR133D.spec")
  
  GeometryList = [0.3,5.0,10.0,300.0,0.0,0.0,0.0,0.0,0.0] #parameters for the geometry correction
  CorrectionInfo = [ "./Data/Mn3Ge-PCR133D.spec", [103,"CL"], [104,"CR"], [105,"LV"], [106,"LH"] ]
    


def GetSampleInfo(fnameCorr,fnameSample):
  #Read the data and info from the spec file into lists
  sInfo,sData = ReadSpecFile(fnameSample)
  
  #s1 = sData[97-1]["PicoAm2"]/sData[97-1]["I0_BD3_Amp"]
  #s2 = sData[102-1]["PicoAm2"]/sData[102-1]["I0_BD3_Amp"]
  
  #np.savetxt("temp.dat",np.transpose(np.vstack((sData[97-1]["TwoTheta"],s1,s2))))
  
  if(fnameSample == datadir + "FGT-2L"):
    EScansList = [ [16,0,0], [17,0,0], [20,0,0], [21,0,0],
                   [23,0,0], [24,0,0], [27,0,0], [28,0,0],
                   [31,0,0], [32,0,0], [35,0,0], [36,0,0],
                   [39,0,0], [40,0,0], [41,0,0], [42,0,0],
                   [45,0,0], [46,0,0], [47,0,0], [48,0,0],
                   [51,0,0], [52,0,0], [53,0,0], [54,0,0],
                   [57,0,0], [58,0,0], [59,0,0], [60,0,0],
                   [63,0,0], [64,0,0], [65,0,0], [66,0,0],
                   [75,0,0], [76,0,0], [77,0,0], [78,0,0],
                   [81,0,0], [82,0,0], [83,0,0], [84,0,0],
                   [87,0,0], [88,0,0], [89,0,0], [90,0,0],
                   [93,0,0], [94,0,0], [95,0,0], [96,0,0],
                   [99,0,0], [100,0,0], [101,0,0], [102,0,0],
                   [105,0,0], [106,0,0], [107,0,0], [108,0,0],
                   [67,0,0], [69,0,0], [71,0,0] ]
    
    AScansList = [ [14,9,0], [15,9,0], [18,9,0], [19,9,0], [22,9,0],
                   [25,9,0], [26,9,0], [29,9,0], [30,9,0], [33,9,0], [34,9,0],
                   [37,9,0], [38,9,0], [43,9,0], [44,9,0], [49,9,0], [50,9,0],
                   [55,9,0], [56,9,0], [61,9,0], [62,9,0], [73,9,0], [74,9,0],
                   [79,9,0], [80,9,0], [85,9,0], [86,9,0], [91,9,0], [92,9,0],
                   [97,9,0], [98,9,0], [103,9,0], [104,9,0], [68,9,0], [70,9,0],
                   [72,9,0] ]
    
    #ECalList = GetECalList([[706.1, 707.25], [706.1+100, 707.25+100]]) #list of pairs of [measE,calE]
    #ECalList = GetECalList([[706.1, 706.1001], [706.1+100, 706.1001+100]]) #list of pairs of [measE,calE]
    ECalList = GetECalList([[707.8,710],[777.8,780]])

    #ECalList = [0, 1, 0, 0, 0 ] #coefficients for polynomial scaling applied to energy (offset, linear, quadratic, etc)
    GeometryList = [0.3,5.0,10.0,300.0,0.0,0.0,0.0,0.0,0.0] #parameters for the geometry correction
    CorrectionInfo = [ fnameCorr, [111,"CL"], [112,"CR"], [109,"LV"], [110,"LH"] ]
    
 
  
  elif(fnameSample == datadir + "FGT-1L"):
    EScansList = [ [15,0,0], [16,0,0], [19,0,0], [20,0,0],
                   [22,0,0], [23,0,0], [26,0,0], [27,0,0],
                   [30,0,0], [31,0,0], [34,0,0], [35,0,0],
                   [36,0,0], [37,0,0], [40,0,0], [41,0,0],
                   [42,0,0], [43,0,0], [46,0,0], [47,0,0],
                   [48,0,0], [49,0,0], [52,0,0], [53,0,0],
                   [54,0,0], [55,0,0], [58,0,0], [59,0,0],
                   [60,0,0], [61,0,0], [67,0,0], [68,0,0],
                   [69,0,0], [70,0,0], [73,0,0], [74,0,0],
                   [75,0,0], [76,0,0], [79,0,0], [80,0,0],
                   [81,0,0], [82,0,0], [62,0,0], [64,0,0] ]
    
    AScansList = [ [13,9,0], [14,9,0], [17,9,0], [18,9,0], [21,9,0],
                   [24,9,0], [25,9,0], [28,9,0], [29,9,0], [32,9,0], [33,9,0],
                   [38,9,0], [39,9,0], [44,9,0], [45,9,0], [50,9,0], [51,9,0],
                   [56,9,0], [57,9,0], [65,9,0], [66,9,0], [71,9,0], [72,9,0],
                   [77,9,0], [78,9,0], [63,9,0] ]
    
    #ECalList = GetECalList([[706.1, 707.25], [706.1+100, 707.25+100]]) #list of pairs of [measE,calE]
    #ECalList = GetECalList([[706.1, 706.1001], [706.1+100, 706.1001+100]]) #list of pairs of [measE,calE]
    ECalList = GetECalList([[707.8,710],[777.8,780]])
    
    #ECalList = [0, 1, 0, 0, 0 ] #coefficients for polynomial scaling applied to energy (offset, linear, quadratic, etc)
    GeometryList = [0.3,5.0,10.0,300.0,0.0,0.0,0.0,0.0,0.0] #parameters for the geometry correction
    CorrectionInfo = [ fnameCorr, [111,"CL"], [112,"CR"], [109,"LV"], [110,"LH"] ]
    
 
  

        
  else:
    print("Unkown filename ", fnameSample)
    exit()
    
  return np.array(EScansList),np.array(AScansList),np.array(ECalList),np.array(GeometryList),CorrectionInfo
    
  
if __name__ == "__main__":
  
  
    fnameCorr = "FGT-2L"
    samples = ["FGT-2L", "FGT-1L" ]

    names = [ ["Ge", "Fe300Ge100Te200Co", "Fe300Ge100Te200Co", "Fe300Ge100Te200Co", "Ge5O", "Ge5O3C"],
              ["Ge", "Fe300Ge100Te200Co", "Fe300Ge100Te200Co", "Fe300Ge100Te200Co", "Ge5O", "Ge5O3C"] ]

    densities = [ [5.323, 7.3, 7.3, 7.3, 5.323, 5.323],
                  [5.323, 7.3, 7.3, 7.3, 5.323, 5.323] ]

    thicknesses = [ [0,  6,  6,  6, 25, 25],
                    [0,  3,  3,  3, 25, 25] ]
                 
    for sam in range(len(samples)):
        EScan,AScan,ECal,Geo,Corr = GetSampleInfo(datadir + fnameCorr, datadir + samples[sam])
        AsData,AsInfo = ProcessRXR(datadir + samples[sam], AScan,ECal,Geo,Corr,"A")
        EsData,EsInfo = ProcessRXR(datadir + samples[sam], EScan,ECal,Geo,Corr,"E")
  
        remagxHeader = GetReMagXHeader(names[sam],densities[sam],thicknesses[sam])
  
    # print(remagxHeader)
    # exit()
    WriteReMagX(samples[sam] + ".all",AsData,AsInfo,EsData,EsInfo,remagxHeader)




    #example of just plotting some data
    #sInfo, sData = ReadSpecFile("Fe3O4-4ML")
    #plscans = [21, 22]
    #PlotSpecData(sData,plscans,"MonoEngy","MCP")

