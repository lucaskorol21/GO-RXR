from math import sin
import numpy as np
from scipy import interpolate

import matplotlib.pyplot as plt





if __name__ == "__main__":
    x = np.arange(0,10,0.01)
    fun = 10*np.sin(x) + 100
    sv = interpolate.splrep(x, fun)
    print(interpolate.splev(5,sv))
    plt.plot()

