import numpy as np
import UTILS.material_structure as ms
import os
import unittest

# This test script can be executed by inputting
#  ->  python -m unittest -v Testing/reflectivity_test.py
# into the terminal

def load_reflections(fname):

    if os.getcwd().split('\\')[-1] == 'Testing':
        my_path = os.getcwd() + '/test_data/' + fname
    else:
        my_path = os.getcwd() + '/Testing/test_data/' + fname

    data = np.loadtxt(my_path)
    theta = data[:,0]
    R = data[:,1]
    return theta,R

class TestReflectivity(unittest.TestCase):
    def test_single_element_sigma(self):
        # Reflectivity retrieved from the Center of X-ray Optics
        #  - calculated for sigma-polarized light

        theta, Rsol = load_reflections('Si-Al-10A-550.txt')  # calculated reflectivity
        #theta, Rsol = load_reflections('AgBr-AlAs-50A-550.txt')  # calculated reflectivity

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

        theta, Rsol = load_reflections('Si-Al-50A-550.txt')  # calculated reflectivity
        #theta, Rsol = load_reflections('AgBr-AlAs-50A-550.txt')  # calculated reflectivity

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

        theta, Rsol = load_reflections('AgBr-AlAs-50A-550.txt')  # calculated reflectivity

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

        theta, Rsol = load_reflections('AgBr-AlAs-50A-550_pi.txt')  # calculated reflectivity

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