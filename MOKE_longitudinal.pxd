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


from Mathematical_Functions_Reflectivity cimport *
cdef extern from "math.h":
    double complex sqrt(double complex)  nogil
    double complex exp(double complex) nogil
    double cos(double) nogil
    double sin(double) nogil

from Structural cimport *
from Mathematical_Functions_Reflectivity cimport *
from Multilayer_Functions_Reflectivity cimport *
from Reflectivity_Sigma cimport *
from Reflectivity_Pi cimport *

cdef void Fill_rMemory(rMemory *Mem, double vy, double vyvy, double vy4, double complex chix, double complex chiy, double complex chiz, double complex chig)

cdef void Calculate_rt_y(rMemory *Mem1, rMemory *Mem2, double vy, double vyvy, double vy4, double complex chix1, double complex chiy1, double complex chiz1, double complex chig1, double complex chix2, double complex chiy2, double complex chiz2, double complex chig2, \
                    int IsMagnetic1, int IsMagnetic2, double complex (*r)[2][2], double complex (*rprime)[2][2], double complex (*t)[2][2], double complex (*tprime)[2][2], double sigma, double k0)

cdef void Paratt_magnetic_y(Heterostructure* HS, double th, double wavelength, double complex (*rtot)[2][2])

cdef void Paratt_magnetic_y_MS(Heterostructure* HS, double th, double wavelength, double complex (*rtot)[2][2])
