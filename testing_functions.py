from math import sin
import numpy as np
from scipy import interpolate
from scipy import optimize




if __name__ == "__main__":
    REAL = np.array([1,2,3,4,5])
    IMAG = np.array([0.1,0.2,0.3,0.4,0.5])

    NEW = np.vectorize(complex)(REAL, IMAG)
    print(NEW)
