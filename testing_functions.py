from math import sin
import numpy as np
from scipy import interpolate
from scipy import optimize
import Pythonreflectivity as pr
from material_model import *






if __name__ == "__main__":

    mag_optical_profile = np.loadtxt('mag_optical_profile')
    optical_profile = np.loadtxt('optical_constant')

    t1 = optical_profile[:,0]
    d1 = optical_profile[:, 1]
    b1 = optical_profile[:, 3]
    print(d1)
    plt.figure()
    plt.plot(t1,b1,'.')
    plt.plot(t1,d1,'.')
    plt.show()
