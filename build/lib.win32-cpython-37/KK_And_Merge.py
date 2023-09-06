#!/usr/bin/python
import numpy as np
from pylab import *


#read from input file
def GetSpecFromFile(filename):
    inFile=open(filename,"r")
    sticks = []
    sticks.append([]) #energy
    sticks.append([]) #f1
    sticks.append([]) #f2
    
    line = inFile.readline()
    while not line == '':
        splitLine = list(filter(None,line.split(' ')))
        sticks[0].append(float(splitLine[0]))
        
        #sticks[1].append(0.0)
        sticks[1].append(float(splitLine[1]))
        sticks[2].append(float(splitLine[2]))
        line = inFile.readline()
    return np.array(sticks)


#merge with offres
def MergeOffRes(spec,offres,ELow,EHigh):
    merged = []
    merged.append([]) #energy
    merged.append([]) #f1
    merged.append([]) #f2

    for i in range(len(offres[0])):
      if(offres[0,i] < ELow):
        merged[0].append(offres[0,i])
        merged[1].append(offres[1,i])
        merged[2].append(offres[2,i])

    for i in range(len(spec[0])):
      if(spec[0,i] >= ELow and spec[0,i] <= EHigh):
        merged[0].append(spec[0,i])
        merged[1].append(spec[1,i])
        merged[2].append(spec[2,i])

    for i in range(len(offres[0])):
      if(offres[0,i] > EHigh):
        merged[0].append(offres[0,i])
        merged[1].append(offres[1,i])
        merged[2].append(offres[2,i])
        
    return np.array(merged)
    
def KK_RobertSlow(e,im,offset):
  re = zeros(len(im)) + offset
  
  for i in range(1,len(im)-1):
    xm = e[i-1]
    x0 = e[i]
    xp = e[i+1]
    beta = im[i]
    
    m1 = -beta/(xp-x0)
    b1 = beta - m1*x0
    
    m2 = beta / (x0-xm)
    b2 = beta - m2*x0
    
    dEnergy = 0.5*(xp-xm)
    
    dist = 10
    if((xp-xm)*3 > dist):
      dist = (xp-xm)*3
      
    for j in range(len(im)):
      delta=0
      newE = e[j]+0.00234563
      if(math.fabs(x0-newE) < dist):
        delta = delta + 0.5 * b1 * math.log(math.fabs((newE-xp)*(newE+xp))/math.fabs((newE-x0)*(newE+x0))) + m1*(xp-x0) + 0.5 * newE * m1 * math.log((newE+x0)/(newE+xp)) + 0.5 * newE * m1 * math.log(math.fabs(newE-xp)/math.fabs(newE-x0))
        delta = delta + 0.5 * b2 * math.log(math.fabs((newE-x0)*(newE+x0))/math.fabs((newE-xm)*(newE+xm))) + m2*(x0-xm) + 0.5 * newE * m2 * math.log((newE+xm)/(newE+x0)) + 0.5 * newE * m2 * math.log(math.fabs(newE-x0)/math.fabs(newE-xm))
      else:
        delta = delta - dEnergy * (beta*x0)/((newE-x0)*(newE+x0))
      re[j] = re[j] - 2.0/3.14159265359 * delta
      
  return re




def KK_Robert(e,im,offset):
#The inner loop here is vectorized compared to the slow version
#The algorithm itself is from Sebastian's online (javascript) kramers kronig calculator
#To do the vectorization I had to remove his "dist" condition where he only uses the log expressions near the current energy point
#The speedup from vectorization was much faster than any slowdown due to using the long log expression for all points

  re = zeros(len(im)) + offset

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

def GetOffset(element):
  switcher = { "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, 
               "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
               "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
               "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
               "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
               "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103,
               "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111 }
  return switcher.get(element, -1000)
    
    
    



