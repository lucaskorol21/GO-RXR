from math import sin
import numpy as np
from scipy import interpolate
from scipy import optimize

import matplotlib.pyplot as plt

def function(step):
    print(step)
def smallest_slab_version1(func, precision):
    #Find the smallest slab by decreasing the size of the slab until the variation reaches a certain precision
    maximum_variation = max(np.abs(np.diff(func)))
    def variation_check(max_var, precision):
        return max_var - precision


if __name__ == "__main__":
    x = np.arange(0,10,0.01)
    fun = 10*np.sin(x) + 100
    sv = interpolate.splrep(x, fun)
    print(interpolate.splev(5,sv))
    plt.plot()

