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


cdef inline double complex CalculateVZpi(double vyvy, double complex cy, double complex cz):
    return sqrt((1.-vyvy/(1.+cz))*(1+cy))

cdef double complex Calculate_rpi_precisely(double vyvy, double complex vz1, double complex vz2, double complex cy1,double complex cy2, double complex cz1, double complex cz2)

cdef double complex LinDicParatt_Pi(Heterostructure* HS, double th, double wavelength)

cdef double complex LinDicParatt_Pi_MS(Heterostructure* HS, double th, double wavelength)
