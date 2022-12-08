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


cdef extern from "math.h":
    double complex sqrt(double complex)  nogil
    double complex exp(double complex) nogil
    double cos(double) nogil
    double sin(double) nogil

from Mathematical_Functions_Reflectivity cimport *
from Structural cimport CLayer, Heterostructure
from Multilayer_Functions_Reflectivity cimport Calculate_Multilayer


cdef inline double complex CalculateVZsigma(double vyvy, double complex cx):
    return sqrt(1.+cx-vyvy)

cdef inline double complex Calculate_rsigma_precisely(double complex vz1, double complex vz2, double complex cx1, double complex cx2):
    return (cx1-cx2)/cquadr(vz1+vz2)


cdef double complex LinDicParatt_Sigma(Heterostructure* HS, double th, double wavelength)

cdef double complex LinDicParatt_Sigma_MS(Heterostructure* HS, double th, double wavelength)
