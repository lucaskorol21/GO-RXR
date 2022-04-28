from math import sin
import numpy as np
import scipy.integrate as integrate

import matplotlib.pyplot as plt


def total_variation(fx,fy):
    diff_x = np.diff(np.array(fx))
    diff_y = np.diff(np.array(fy))

    total_var = np.sum(np.abs(diff_x)) + np.sum(np.abs(diff_x))
    return total_var

if __name__ == '__main__':
    x = np.arange(1, 100+0.001, 0.001)
    fx = np.sin(x) + np.cos(x)
    fy = x
    r = x[-1]-x[0]
    tv = total_variation(fx, fy)
    z = np.arange(1,100+0.001,r/tv)
    plt.figure()
    plt.plot(np.sin(x)+ np.cos(x), x)
    plt.show()