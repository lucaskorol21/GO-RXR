import numpy as np
import matplotlib.pyplot as plt
import Pythonreflectivity as pr




print("One interface")

A=pr.Generate_structure(1)
A[0].seteps(1.01+0.01j) # material epsilon
A[0].setsigma(2) # interface roughness
Theta=np.linspace(0.1, 89.9, 899) # Angles
wavelength=10 # Wavelength (same unit as roughness)

R=pr.Reflectivity(A, Theta, wavelength)

Sigma, = plt.plot(Theta, np.log(R[0]), label='Sigma')
Pi, = plt.plot(Theta, np.log(R[1]), label='Pi')
plt.legend(handles=[Sigma, Pi])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('One interface')
plt.show()

print("Two layers on a substrate")

A=pr.Generate_structure(3)


A[0].seteps(1.01+0.01j) # material epsilon
A[1].seteps(0.99+0.015j)
A[2].seteps(0.998+0.002j)
A[1].setd(40) # Layer thickness
A[2].setd(500)
A[0].setsigma(0) # interface roughness (0 is actually the default value so this is not needed)
A[1].setsigma(0)
A[2].setsigma(0)

R=pr.Reflectivity(A, Theta, wavelength)

Sigma, = plt.plot(Theta, np.log(R[0]), label='Sigma')
Pi, = plt.plot(Theta, np.log(R[1]), label='Pi')
plt.legend(handles=[Sigma, Pi])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('Two layers on a substrate')
plt.show()

print("An anisotropic layer")
A=pr.Generate_structure(2)
A[0].seteps(1.01+0.01j) # material epsilon
A[1].seteps([1.01+0.01j, 1.01+0.01j, 0.99+0.015j])
A[1].setd(40)

R=pr.Reflectivity(A, Theta, wavelength)

Sigma, = plt.plot(Theta, np.log(R[0]), label='Sigma')
Pi, = plt.plot(Theta, np.log(R[1]), label='Pi')
plt.legend(handles=[Sigma, Pi])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('An anisotropic layer')
plt.show()

print("A magnetic layer")

A=pr.Generate_structure(2)
A[0].seteps(1.01+0.01j) 
A[1].setmag("z")
A[1].seteps([0.99+0.015j, 0.99+0.015j, 0.99+0.015j, 0.003+0.003j])
A[1].setd(40)

R=pr.Reflectivity(A, Theta, wavelength)

Sigma, = plt.plot(Theta, np.log(R[0]), label='Sigma')
Pi, = plt.plot(Theta, np.log(R[1]), label='Pi')
Left, = plt.plot(Theta, np.log(R[2]), label='Left circular')
Right, = plt.plot(Theta, np.log(R[3]), label='Right circular')
plt.legend(handles=[Sigma, Pi, Left, Right])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('A magnetic layer')
plt.show()


print("Polarization analysis")

A=pr.Generate_structure(2)
A[0].seteps(1.01+0.01j) 
A[1].setmag("z")
A[1].seteps([0.99+0.015j, 0.99+0.015j, 0.99+0.015j, 0.003+0.003j])
A[1].setd(40)

R=pr.Reflectivity(A, Theta, wavelength, Output="t")

RSS, = plt.plot(Theta, np.log(abs(R[0])**2), label='Sigma to Sigma')
RPP, = plt.plot(Theta, np.log(abs(R[3])**2), label='Pi to Pi')
RSP, = plt.plot(Theta, np.log(abs(R[1])**2), label='Sigma to Pi = Pi to Sigma')

plt.legend(handles=[RSS, RPP, RSP])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('Polarization analysis')
plt.show()


print("Multilayers")

A=pr.Generate_structure(3, "0, 4*(1,2)")
A[0].seteps(1.01+0.01j) 
A[1].seteps([0.995+0.01j, 0.995+0.01j, 0.99+0.015j])
A[1].seteps([0.92+0.03j, 0.92+0.03j, 0.91+0.011j])
A[1].setd(15)
A[2].setd(25)

R=pr.Reflectivity(A, Theta, wavelength)

Sigma, = plt.plot(Theta, np.log(R[0]), label='Sigma')
Pi, = plt.plot(Theta, np.log(R[1]), label='Pi')
plt.legend(handles=[Sigma, Pi])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('Multilayers')
plt.show()

print("With and without multiple scattering")

A=pr.Generate_structure(3)
A[0].seteps(0.82+0.001)

A[1].seteps([0.82+0.04j])
A[2].seteps([0.95+0.02j])
A[1].setd(90)
A[2].setd(35)
R1=pr.Reflectivity(A, Theta, wavelength, MultipleScattering=True)
R2=pr.Reflectivity(A, Theta, wavelength, MultipleScattering=False)

Sigma, = plt.plot(Theta, np.log(R1[0]), label='Sigma exact')
Pi, = plt.plot(Theta, np.log(R1[1]), label='Pi exact')
Left, = plt.plot(Theta, np.log(R2[0]), label='Sigma approximated')
Right, = plt.plot(Theta, np.log(R2[1]), label='Pi approximated')
plt.legend(handles=[Sigma, Pi, Left, Right])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('With and without multiple scattering')
plt.show()





