import numpy as np
import material_structure as ms

def test_simple_model():
    sample = ms.slab(2)

    E = 833  # energy

    # no roughness
    sample.addlayer(0,'Sr',50)  # Sr substrate
    sample.addlayer(1,'La',50)  # La film

    qz = np.linspace(0.01,0.8,100)
    Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / np.pi  # initial angle

    R = sample.reflectivity(E, qz,precision=1e-20)

    # calculate model analytically

    # find the refractive index

    # constructing the propagation matrices
    n1 = 1
    n2 = 1
    vy = np.cos(Theta)  #vy
    v_z_s_p = 0
    v_z_s_m = 0
    v_z_p_p = 0
    v_z_p_m = 0
    pass

