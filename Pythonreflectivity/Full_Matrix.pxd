#<Pythonreflectivity: A Python Package for simulation of x-ray reflectivities of Heterostructures>
#    Copyright (C) <2017>  <Martin Zwiebler>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.



from Structural cimport CLayer, Heterostructure, MatrixSafer
from Mathematical_Functions_Reflectivity cimport *
from Reflectivity_Sigma cimport CalculateVZsigma
from Reflectivity_Pi cimport CalculateVZpi
from MOKE_transversal cimport CalculateVZpi_m
from Multilayer_Functions_Reflectivity cimport Matrixexp
#from scipy.linalg cimport eig as scipydiagonalization
#from numpy.linalg cimport eig als npeig
cdef extern from "math.h":
    double complex sqrt (double complex) nogil
    double complex pow(double complex, double) nogil
    double complex exp(double complex) nogil
    double sin(double) nogil
    double cos(double) nogil

cdef void Reduce_complexity_of_chi(CLayer *Layer, double Cutoff, int *allx, int *ally, int *allz)

cdef void Fill_Matrixsafer( MatrixSafer *MS, CLayer L )

cdef void NormalizePHI(double complex (*PHI)[4][4] )

cdef void NormalizePSI(double complex (*PSI)[4][4] )

cdef void Calculate_Phi_and_Psi(CLayer L, MatrixSafer *MS, double vy, double vzvz, double vyvy, double complex (*vz)[4], double complex (*PHI)[4][4], double complex (*PSI)[4][4])

cdef void Full_Matrix(Heterostructure* HS, MatrixSafer *AllMS, int* Layer_type_to_Matrixsafe, double th, double wavelength, double complex (*rtot)[2][2])
