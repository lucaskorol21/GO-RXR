import Pythonreflectivity as pr
import numpy as np
import matplotlib.pyplot as plt





print( "One interface" )

A=pr.Generate_structure(1)
A[0].setchi(0.01+0.01j) # material chi
A[0].setsigma(2) # interface roughness
Theta=np.linspace(0.1, 89.9, 899) # Angles in degrees
wavelength=10 # Wavelength (same unit as roughness and thickness)

R=pr.Reflectivity(A, Theta, wavelength)

Sigma, = plt.plot(Theta, np.log(R[0]), label='Sigma')
Pi, = plt.plot(Theta, np.log(R[1]), label='Pi')
plt.legend(handles=[Sigma, Pi])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('One interface')
plt.show()

print( "Two layers on a substrate" )

A=pr.Generate_structure(3)


A[0].setchi(0.01+0.01j) # material chi (susceptibility)
A[1].setchi(-0.01+0.015j)
A[2].setchi(-0.002+0.002j)
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

print( "An anisotropic layer" )
A=pr.Generate_structure(2)
A[0].setchi(0.01+0.01j) # material chi
A[1].setchi([0.01+0.01j, 0.01+0.01j, -0.01+0.015j]) # Along x, y and z
A[1].setd(40)

R=pr.Reflectivity(A, Theta, wavelength)

Sigma, = plt.plot(Theta, np.log(R[0]), label='Sigma')
Pi, = plt.plot(Theta, np.log(R[1]), label='Pi')
plt.legend(handles=[Sigma, Pi])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('An anisotropic layer')
plt.show()

print( "A magnetic layer" )

A=pr.Generate_structure(2)
A[0].setchi(0.01+0.01j) 
A[1].setmag("z") # Magnetization in this direction
A[1].setchi([-0.01+0.015j, -0.01+0.015j, -0.01+0.015j, 0.003+0.003j]) # Along x, y and z and the magnetic contribution (Polar magnetooptical Kerr effect in this case)
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


print( "Polarization analysis" )

A=pr.Generate_structure(2)
A[0].setchi(0.01+0.01j) 
A[1].setmag("z") # Magnetization in this direction
A[1].setchi([-0.01+0.015j, -0.01+0.015j, -0.01+0.015j, 0.003+0.003j])# Along x, y and z and the magnetic contribution (Polar magnetooptical Kerr effect in this case)
A[1].setd(40)

R=pr.Reflectivity(A, Theta, wavelength, Output="t") #Output="t" means complex amplitudes instead of intensities

RSS, = plt.plot(Theta, np.log(abs(R[0])**2), label='Sigma to Sigma')
RPP, = plt.plot(Theta, np.log(abs(R[3])**2), label='Pi to Pi')
RSP, = plt.plot(Theta, np.log(abs(R[1])**2), label='Sigma to Pi = Pi to Sigma')

plt.legend(handles=[RSS, RPP, RSP])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('Polarization analysis')
plt.show()


print( "Multilayers" )

A=pr.Generate_structure(3, "0, 1000*(1,2)") # Repeat this combination (Layer 1, Layer2) 1000 times
A[0].setchi(0.01+0.01j) 
A[1].setchi([-0.005+0.01j, -0.005+0.01j, -0.001+0.015j])
A[2].setchi([-0.015+0.003j, -0.015+0.003j, -0.015+0.003j])
A[1].setd(45)
A[2].setd(5)

R=pr.Reflectivity(A, Theta, wavelength)

Sigma, = plt.plot(Theta, np.log(R[0]), label='Sigma')
Pi, = plt.plot(Theta, np.log(R[1]), label='Pi')
plt.legend(handles=[Sigma, Pi])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('Multilayers')
plt.show()

print( "With and without multiple scattering" )

A=pr.Generate_structure(3)
A[0].setchi(-0.002+0.001j)

A[1].setchi([-0.004+0.04j])
A[2].setchi([-0.005+0.02j])
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


print( "Phase analysis" )
A=pr.Generate_structure(3)
A[0].setchi( -0.002+0.001j )
A[1].setchi(-0.004+0.04j)
A[2].setchi(-0.005+0.02j)
A[1].setd(60)
A[2].setd(35)
r=pr.Reflectivity(A, Theta, wavelength, Output='t') #Output="t" means complex amplitudes instead of intensities
phase_sigma=np.arctan2( np.imag(r[0]), np.real(r[0]) )
phase_pi=np.arctan2( np.imag(r[1]), np.real(r[1]) )
phase_sigma-=phase_sigma[0] #choose a reference phase
phase_pi-=phase_pi[0]

Sigma, = plt.plot(Theta, phase_sigma, label='Sigma')
Pi, = plt.plot(Theta, phase_pi, label='Pi')

yaxisticks=[ 0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi ]

plt.yticks( yaxisticks, ( r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$') )
plt.legend(handles=[Sigma, Pi])
plt.xlabel('Theta')
plt.ylabel('Phase change during scattering')
plt.title("Phase analysis")
plt.show()


print( "Full Matrix" )
A=pr.Generate_structure(4, MLstructure="0,4*(1,2,3,2),5*(3,1),2")

chi1=0.003+0.004j
#set chixx, chixy, chixz, chiyx, chiyy, chiyz, chizx, chizy, chizz
chi2=np.array([0.003+0.001j,-0.0002+0.0001j , -0.0001-0.0004j , \
     0.0003+0.0001j, 0.001+0.001j, -0.0004-0.0003j,\
     0.0001+0.0002j ,0.0003+0.0008j, 0.002+0.001j ])

chi3=[0.004+0.001j]
chi4=0.25*chi2

A[0].setchi(chi1)
A[0].setsigma(3)
A[1].setchi(chi2)
A[1].setd(4)
A[1].setsigma(3.5)
A[2].setchi(chi3)
A[2].setd(10)
A[2].setsigma(4)
A[3].setchi(chi4)
A[3].setd(5)
A[3].setsigma(4.5)
R=pr.Reflectivity(A, Theta, wavelength)
Sigma, = plt.plot(Theta, np.log(R[0]), label='Sigma')
Pi, = plt.plot(Theta, np.log(R[1]), label='Pi')
Left, = plt.plot(Theta, np.log(R[2]), label='Left circular')
Right, = plt.plot(Theta, np.log(R[3]), label='Right circular')
plt.legend(handles=[Sigma, Pi, Left, Right])
plt.xlabel('Theta')
plt.ylabel('log(Reflectivity)')
plt.title('Full matrix')
plt.show()

print( "Kramers-Kronig transformations" )
N=1001
x=np.linspace(920,1080,N)
s1=3
s2=4
t1=980
t2=1010 #make some Cauchy-Lorenz peaks
y=10*s1/(s1*s1+(t1-x)**2 )+10*s2/(s2*s2+(t2-x)**2 )

z=pr.KramersKronig(x, y) ##only imaginary to real part is possible

Imaginary, = plt.plot(x, y, label='Imag. part')
Real, = plt.plot(x, z, label='Real part')
plt.xlabel('Energy')
plt.ylabel('Form factors')

plt.legend(handles=[Imaginary, Real])
plt.title("Kramers-Kronig Transformation")
plt.show()
