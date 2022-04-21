import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

def spacedmarks(x, y, Nmarks, data_ratio=None):
    import scipy.integrate

    if data_ratio is None:
        data_ratio = plt.gca().get_data_ratio()

    dydx = np.gradient(y, x[1])
    dxdx = np.gradient(x, x[1])*data_ratio
    arclength = scipy.integrate.cumtrapz(np.sqrt(dydx**2 + dxdx**2), x, initial=0)
    marks = np.linspace(0, max(arclength), Nmarks)
    markx = np.interp(marks, arclength, x)
    marky = np.interp(markx, x, y)

    return markx, marky

x = np.linspace(0, 10*np.pi,1000)
y = np.sin(x*2) + np.sin(x+1)

markx, marky = spacedmarks(x,y, 80)

plt.figure()
plt.plot(markx, marky, 'o', color='blue')
plt.plot(x,y)
plt.show()