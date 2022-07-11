from material_structure import *
import pickle
from numba import *
from scipy import interpolate


# Loads all scattering factors when program imported
with open('ff_Altered.pkl', 'rb') as f:
    global ff
    ff = pickle.load(f)
f.close()

# Loads all scattering factors when program imported
with open('ffm_Altered.pkl','rb') as f:
    global ffm
    ffm = pickle.load(f)
f.close()


def resetAlteredSF():
    """
    Purpose: Reset original form factors to original values
    :return:
    """
    global ff
    with open('form_factor.pkl', 'rb') as f:
        ff = pickle.load(f)

    with open('ff_Altered.pkl', 'wb') as handle:
        pickle.dump(ff, handle)

def resetSF():
    """
        Purpose: Reset original form factors to original values
        :return:
        """
    global ff
    with open('ff_Altered.pkl', 'rb') as f:
        ff = pickle.load(f)



def resetAlteredSFM():
    """
    Purpose: Reset original form factors to original value
    :return:
    """
    global ffm
    with open('form_factor_magnetic.pkl', 'rb') as f:
        ffm = pickle.load(f)

    with open('ffm_Altered.pkl', 'wb') as handle:
        pickle.dump(ffm, handle)


def resetSFM():
    """
    Purpose: Reset original form factors to original values
    :return:
    """
    global ffm
    with open('ffm_Altered.pkl', 'rb') as f:
        ffm = pickle.load(f)

def FfEnergyShift(element, dE, opt=False):
    """
    Purpose: set the energy shift for the form factor of a specified element
    :param element: the element symbol
    :param dE: the energy shift in eV
    :param opt: boolean that determines if you want to optimize the energy shift
                    True - use the optimization capability
                    False - set the form factor to the new shifted value
    :return:
    """
    global ff
    if not(opt):  # no optimization
        ff[element][:,0] = ff[element][:,0] + dE  #energy shift

        # save shifted value to file
        with open('ff_Altered.pkl') as f:
            pickle.dump(ff, f)
        f.close()
    else:  # optimization
        resetSF()
        ff[element][:,0] = ff[element][:,0] + dE



def FfmEnergyShift(element,dE, opt = False):
    """
    Purpose: set the energy shift for the magnetic scattering factor
    :param element: symbol for desired element to shift
    :param dE: desired energy shift
    :param opt: boolean that specifies if energy shift will be optimized
                    True - optimization
                    False - no-optimization and save to altered form factor file
    :return:
    """

    global ffm
    if not(opt): # non-optimization
        ffm[element][:,0] = ffm[element][:,0] + dE  # energy shift

        with open("ffm_Altered.pkl") as f:  # save shifted magnetic form factor to altered file
            pickle.dump(ffm, f)
        f.close()
    else:  # optimization
        resetSFM()
        ffm[element][:,0] = ffm[element][:,0] + dE  # energy shift

def form_factor(f,E):

    """
    Purpose: Determines form factors with energy E using linear interpolation
    :param f: List of form factors
    :param E: Desired energy
    :return: Array contains real and imaginary component of form factor at energy E: f=[real, imaginary]
    """
    # Linear interpolation
    fr = interpolate.interp1d(f[:,0],f[:,1])
    fi = interpolate.interp1d(f[:,0],f[:,2])

    if isinstance(E, list) or isinstance(E, np.ndarray):  # handle multiple energy case
        F = np.array([np.array([fr(x), fi(x)]) if x > f[0, 0] and x < f[-1, 0] else np.array([0, 0]) for x in E])
    else:  # handle single energy case
        F = np.array([fr(E), fi(E)]) if E>f[0,0] and E<f[-1,0] else np.array([0,0])
    return F

def find_form_factor(element, E, mag):
    """
    Purpose: Return the magnetic or non-magnetic form factor of a selected element and energy
    :param element: String containing element symbol
    :param E: Energy in electron volts
    :param mag: Boolean specifying if the magnetic form factor is desired
    :return:
    """
    global ffm
    global ff
    if mag:
        mag_keys = list(ffm.keys())
        if element not in mag_keys:
            raise NameError(element + " not found in magnetic form factors")
        F = form_factor(ffm[element],E)
    else:
        struc_keys = list(ff.keys())
        if element not in struc_keys:
            raise NameError(element + " not found in structural form factors")
        F = form_factor(ff[element], E)

    return F

def find_form_factors(element, E, mag):
    """
    Purpose: Retrieve form factor from database
    :param element: String containing element symbol
    :param E: Float or integer of desired energy in units of eV
    :param mag: Boolean
                    True - Magnetic form factor
                    False - Non-magnetic form factor
    :return: Return form factor at energy 'E'
    """
    F = 0  # pre-initialized form factor

    if mag:  # looking for magnetic form factors
        my_dir = os.getcwd() + r'\Magnetic_Scattering_Factor'
    else:  # looking for non-magnetic form factors
        my_dir = os.getcwd() + r'\Scattering_Factor'

    for ifile in os.listdir(my_dir):

        if ifile.endswith(element + '.txt'):
            F = form_factor(np.loadtxt(my_dir +  "\\" + ifile),E)
    return F

def MOC(rho, sfm, E, n):
    """
    Purpose: computes the magneto-optical constant for the energy scan
    :param rho: dictionary containing the element symbol as the key and a numpy array as the value
    :param sfm: dictionary that contains the element symbol as the key and the absorptive and dispersive form factor components
    :param E: a numpy array containing energy values in eV
    :return: The absorptive and dispersive magnetic-optical constants
    """
    # Constants
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    k0 = 2 * pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]

    constant = 2 * pi * re * (avocado) / (k0 ** 2)  # constant for density sum



    elements = list(rho.keys())  # retrieves all the magnetic elements in the layer

    delta_m = np.array([np.zeros(n) for x in range(len(E))])  # pre-initialization
    beta_m = np.array([np.zeros(n) for x in range(len(E))])  # pre-initialization
    # Computes the dispersive and absorptive components of the magnetic-optical constant using list comprehensions
    for element in elements:
        delta_m = delta_m + np.array(
            [constant[x] * sfm[element][x, 0] * rho[element] for x in range(len(sfm[element][:, 0]))])
        beta_m = beta_m + np.array(
            [constant[x] * sfm[element][x, 1] * rho[element] for x in range(len(sfm[element][:, 1]))])


    return delta_m, beta_m

def magnetic_optical_constant(rho, sfm, E):
    """
    Purpose: Calculate the magnetic optical constants
    :param rho: Magnetic density in mol/cm^3
    :param sf: dictionary relates elements to scattering factor sf = {'ele1':'ffm1',...,'eleN':'ffmN'}
    :param E: Desired energy in units of eV
    :return: delta_m - magneto-optic dispersive component
             beta_m  - magneto-optic absorptive component
    """

    # Constants
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    k0 = 2 * pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]

    constant = 2 * pi * re * (avocado) / (k0 ** 2)  # constant for density sum

    f1 = 0  # for dispersive component computation
    f2 = 0  # for absorptive component computation

    elements = list(rho.keys())  # retrieves all the magnetic elements in the layer

    # Computes the dispersive and absorptive components of the index of refraction
    for element in elements:
        f1 = f1 + sfm[element][0] * rho[element]
        f2 = f2 + sfm[element][1] * rho[element]



    delta_m = constant * f1  # dispersive component
    beta_m = constant * f2  # absorptive component


    return delta_m, beta_m

def IoR(rho,sf,E):
    """
    Purpose: compute the refractive index for multiple energies
    :param rho: dictionary containing element symbol as key and numpy array as value
    :param sf: dictionary containing element symbol as key and numpy array of dispersive and absorptive form factors
    :param E: numpy array of energies (eV)
    :return: The absorptive and dispersive components of the refractive index
    """
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    k0 = 2 * pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]

    constant = 2 * pi * re * (avocado) / (k0 ** 2)  # constant for density sum

    elements = list(rho.keys())  # retrieves all the magnetic elements in the layer
    delta = np.array([np.zeros(len(rho[elements[0]])) for x in range(len(E))])  # initialization
    beta = np.array([np.zeros(len(rho[elements[0]])) for x in range(len(E))])  # initialization

    # Computes the dispersive and absorptive components of the index of refraction using list comprehensions
    for element in elements:
        delta = delta + np.array([constant[x]*sf[element][x, 0] * rho[element] for x in range(len(sf[element][:, 0]))])
        beta = beta + np.array([constant[x] * sf[element][x, 1] * rho[element] for x in range(len(sf[element][:, 1]))])
    return delta, beta

def index_of_refraction(rho, sf, E):
    """
    Purpose: Calculates the dispersive and absorptive components of the index of refraction
    :param rho: Dictionary containing density profile of elements {'ele1':rho1,...,'eleN':rhoN}
    :param sf: Dictionary of scattering factors {'ele1':ff1, ... , 'eleN':ffN}
    :param E: Desired energy in units of eV
    :return: delta - dispersive component of the refractive index
             beta - absorptive component of the refractive index
    """
    mag = False  # statement for retrieval of non=magnetic form factors
    # Constants
    h = 4.135667696e-15  # Plank's Constant [eV s]
    #h = 4.1357e-15
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    #c = 2.9979e10
    #re = 2.8179e-13
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    #avocado = 6.0221e23
    k0 = 2 * pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]
    constant = 2 * pi * re * (avocado) / (k0 ** 2)  # constant for density sum

    f1 = 0  # dispersive form factor
    f2 = 0  # absorptive form factor

    elements = list(rho.keys())  # retrieves element symbols within layer
    """
    F = dict()  # dictionary used to contain the form factors
    for element in elements:
        F[ element] = find_form_factor(sf[element],E, mag)
    """
    #  Computes the dispersive and absorptive components
    for element in elements:
        f1 = f1 + sf[element][0] * rho[element]
        f2 = f2 + sf[element][1] * rho[element]

    delta = f1*constant  # dispersive component
    beta = f2*constant  # absorptive component

    return delta, beta




if __name__ == "__main__":

    """
    E = 800 # Xray Energy
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e18  # Speed of light in vacuum [A/s]

    Theta = np.linspace(0.1, 89.9, 899)  # Angles
    wavelength = (h*c)/E  # Wavelength (same unit as roughness) (Angstroms or nm)
    test = np.loadtxt('test_example.txt')  # Data from ReMagX

    sample = slab(2)

    sample.addlayer(0, 'Fe', 50, density=1.56366)
    sample.addlayer(1, 'Fe', 38, density=1.56366)

    thickness, density, mag_density = sample.density_profile()
    eps = dielectric_constant(density, E, mag)

    A = pr.Generate_structure(2)  # initializes slab structure
    A[0].seteps(eps[0])  # creates the substrate layer
    A[1].seteps(eps[0])  # creates film layer
    A[1].setd(38)  # sets thickness


    R1 = pr.Reflectivity(A, Theta, wavelength, MultipleScattering=True)  # Computes the reflectivity

    plt.figure()
    qz = (0.001013546247)*E*sin(Theta*pi/180)
    Sigma, = plt.plot(qz, R1[0], 'k-',label='Python')
    plt.yscale("log")
    plt.xlabel('qz')
    plt.ylabel('Reflectivity')
    plt.title('ReMagX vs. Python Script (800 eV)')
    plt.show()
    """
    #my_dir = os.getcwd() + r'\Magnetic_Scattering_Factor'
    #file = my_dir + "\\" + "Ni.txt"
    #print(file)
    #np.loadtxt(file)

    resetAlteredSF()

