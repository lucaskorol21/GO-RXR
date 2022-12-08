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
    double log(double)  nogil

from Mathematical_Functions_Reflectivity cimport dabsvalue
from numpy cimport ndarray


cdef double KKT_Onepoint(  int j,  double *e,  double *im, int Lred )

cdef void Kramers_Kronig_Transformation_internal( ndarray[double, ndim=1, mode="c"] e, \
                                         ndarray[double, ndim=1, mode="c"] im,
                                         ndarray[double, ndim=1, mode="c"] re, int L)
