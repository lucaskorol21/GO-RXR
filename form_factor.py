import time
import numpy as np
import os
import pickle

ff = {'H': 0,
      'He': 0,
      'Li': 0,
      'Be': 0,
      'B': 0,
      'C': 0,
      'N': 0,
      'O': 0,
      'F': 0,
      'Ne': 0,
      'Na': 0,
      'Mg': 0,
      'Al': 0,
      'Si': 0,
      'P': 0,
      'S' : 0,
      'Cl': 0,
      'Ar': 0,
      'K': 0,
      'Ca': 0,
      'Sc': 0,
      'Ti': 0,
      'V': 0,
      'Cr': 0,
      'Mn': 0,
      'Fe': 0,
      'Co': 0,
      'Ni': 0,
      'Cu': 0,
      'Zn': 0,
      'Ga': 0,
      'Ge': 0,
      'As': 0,
      'Se': 0,
      'Br': 0,
      'Kr': 0,
      'Rb': 0,
      'Sr': 0,
      'Y': 0,
      'Zr': 0,
      'Nb': 0,
      'Mo': 0,
      'Tc': 0,
      'Ru': 0,
      'Rh': 0,
      'Pd': 0,
      'Ag': 0,
      'Cd': 0,
      'In': 0,
      'Sn': 0,
      'Sb': 0,
      'Te': 0,
      'I': 0,
      'Xe': 0,
      'Cs': 0,
      'Ba': 0,
      'La': 0,
      'Hf': 0,
      'Ta': 0,
      'W': 0,
      'Re': 0,
      'Os': 0,
      'Ir': 0,
      'Pt': 0,
      'Au': 0,
      'Hg': 0,
      'Tl': 0,
      'Pb': 0,
      'Bi': 0,
      'Po': 0,
      'At': 0,
      'Rn': 0,
      'Fr': 0,
      'Ra': 0,
      'Ac': 0,
      'Ce': 0,
      'Pr': 0,
      'Nd': 0,
      'Pm': 0,
      'Sm': 0,
      'Eu': 0,
      'Gd': 0,
      'Tb': 0,
      'Dy': 0,
      'Ho': 0,
      'Er': 0,
      'Tm': 0,
      'Yb': 0,
      'Lu': 0,
      'Th': 0,
      'Pa': 0,
      'U': 0}

F = 0  # pre-initialized form factor

my_dict = {}
my_dir = os.getcwd() + r'\Magnetic_Scattering_Factor'

elements = ff.keys()
for element in elements:
      for ifile in os.listdir(my_dir):

            if ifile.endswith(element + '.txt'):
                  F = np.loadtxt(my_dir +  "\\" + ifile)
                  my_dict[element] = F


with open('form_factor_magnetic.pkl', 'wb') as f:
      pickle.dump(my_dict, f)

