from math import sin
import numpy as np
import scipy.optimize as optimize

import matplotlib.pyplot as plt


def zero_to_one(func):
    func_min = min(func)
    func_max = max(func)
    amplitude = func_max - func_min

    return (func - func_min) / amplitude


if __name__ == "__main__":
    x = np.arange(0,10,0.01)
    fun = 10*np.sin(x) + 100

    new_func = zero_to_one(fun)

    plt.figure(1)
    plt.plot(x, new_func)
    plt.show()

