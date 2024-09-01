import os
import sys

# Get the parent directory of the current script's directory
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# Add the parent directory to the system path
sys.path.append(parent_dir)

import numpy as np
import UTILS.material_structure as ms
from UTILS import ROOT_DIR
import unittest

# This test script can be executed by inputting
#  ->  python -m unittest -v test_reflectivity.py
# into the terminal

def load_reflections(my_path):

    data = np.loadtxt(my_path)
    theta = data[:,0]
    R = data[:,1]
    return theta,R

class TestReflectivity(unittest.TestCase):

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)

        self.root_dir = ROOT_DIR

        self.filename1 = 'Si-Al-10A-550.txt'
        self.my_path1 = self.root_dir + '/TESTS/test_data/' + self.filename1 

        self.filename2 = 'Si-Al-50A-550.txt'
        self.my_path2 = self.root_dir + '/TESTS/test_data/' + self.filename2

        self.filename3 = 'AgBr-AlAs-50A-550.txt'
        self.my_path3 = self.root_dir + '/TESTS/test_data/' + self.filename3

        self.filename4 = 'AgBr-AlAs-50A-550_pi.txt'
        self.my_path4 = self.root_dir + '/TESTS/test_data/' + self.filename4


    def test_single_element_sigma(self):

        # Reflectivity retrieved from the Center of X-ray Optics
        #  - calculated for sigma-polarized light

        theta, Rsol = load_reflections(self.my_path1)  # calculated reflectivity

        sample = ms.slab(2)
        sample.addlayer(0, 'Si', 50, density=0.028)
        sample.addlayer(1,'Al',10, density=0.028)


        E = 550

        qz = np.sin(theta * np.pi / 180) * (E * 0.001013546143)

        qz, Rtest = sample.reflectivity(E,qz,precision=1e-20)

        Rtest = Rtest['S']

        total = sum(abs(np.log(Rsol)-np.log(Rtest)))/len(theta)

        self.assertTrue(total < 0.002)

    def test_single_element_larger_film_sigma(self):
        # Reflectivity retrieved from the Center of X-ray Optics
        #  - calculated for sigma-polarized light

        theta, Rsol = load_reflections(self.my_path2)  # calculated reflectivity

        sample = ms.slab(2)
        sample.addlayer(0, 'Si', 50, density=0.028)
        sample.addlayer(1,'Al',50, density=0.028)

        E = 550

        qz = np.sin(theta * np.pi / 180) * (E * 0.001013546143)

        qz, Rtest = sample.reflectivity(E,qz,precision=1e-20)

        Rtest = Rtest['S']

        total = sum(abs(np.log(Rsol)-np.log(Rtest)))/len(theta)

        self.assertTrue(total < 0.002)

    def test_compound_film_sigma(self):

        # Reflectivity retrieved from the Center of X-ray Optics
        #  - calculated for sigma-polarized light

        theta, Rsol = load_reflections(self.my_path3)  # calculated reflectivity

        sample = ms.slab(2)
        sample.addlayer(0, 'AlAs', 50, density=3.81)
        sample.addlayer(1, 'AgBr', 50, density=6.473)


        E = 550

        qz = np.sin(theta * np.pi / 180) * (E * 0.001013546143)

        qz, Rtest = sample.reflectivity(E,qz,precision=1e-20)

        Rtest = Rtest['S']

        total = sum(abs(np.log(Rsol)-np.log(Rtest)))/len(theta)

        self.assertTrue(total < 0.002)

    def test_compound_film_pi(self):
        # Reflectivity retrieved from the Center of X-ray Optics
        #  - calculated for sigma-polarized light

        theta, Rsol = load_reflections(self.my_path4)  # calculated reflectivity

        sample = ms.slab(2)
        sample.addlayer(0, 'AlAs', 50, density=3.81)
        sample.addlayer(1, 'AgBr', 50, density=6.473)


        E = 550

        qz = np.sin(theta * np.pi / 180) * (E * 0.001013546143)

        qz, Rtest = sample.reflectivity(E,qz,precision=1e-20)

        Rtest = Rtest['P']

        total = sum(abs(np.log(Rsol)-np.log(Rtest)))/len(theta)

        self.assertTrue(total < 0.002)


if __name__ == "__main__":

    unittest.main()