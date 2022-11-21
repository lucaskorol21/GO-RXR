import numpy as np




if __name__ == '__main__':
    fname = 'Te_R.ff'

    data = np.loadtxt(fname)

    new_data = data[:,0:3]
    E = data[:,0]
    alpha = data[:,1]
    beta = data[:,2]
    new_data = data[:, 0:3]
    new_data[:,2] = data[:,3]
    np.savetxt(fname,new_data)

