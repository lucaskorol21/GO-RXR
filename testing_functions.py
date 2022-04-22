from math import sin
import numpy as np
import scipy.integrate as integrate

import matplotlib.pyplot as plt


def slice_size(f,a,dt,dt_min,precision):
    left = (f(a + dt / 2) + 1)* dt / 2
    right = (f(a + dt) + 1)* dt/2
    whole = (f(a + dt) + 1)* dt

    delta = (left+right-whole)/dt


    while(abs(delta)<precision):
        dt = dt + dt_min
        left = (f(a + dt / 2) + 1) * dt / 2
        right = (f(a + dt) + 1) * dt/2
        whole = (f(a + dt) + 1)* dt

        delta = (left+right-whole)/dt

    print(delta, dt)
    if dt == dt_min:
        return dt
    else:
        return dt-dt_min









if __name__ == '__main__':
    a = 0
    b = 2 * np.pi

    dt = 1e-3
    dt_min = 1e-3
    precision = 1e-1
    f = sin

    slicing = list()
    slicing.append(a)
    while (a < b):
        dt = dt_min
        dt = slice_size(f, a, dt, dt_min, precision)
        a = a + dt
        slicing.append(a)

    val = np.sin(np.array(slicing)) + 1

    plt.figure()
    plt.plot(slicing, val, '.')
    plt.show()
