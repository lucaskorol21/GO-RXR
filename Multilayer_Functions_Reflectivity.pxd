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

cdef void FillC0(double complex (*C0)[2][2], double complex  (*rprime)[2][2], double complex (*rtot)[2][2],double complex  (*p)[2][2])


cdef void Calculate_ANXBN(double complex (*A)[2][2], double complex (*B)[2][2], double complex (*X)[2][2], int N)

cdef void Calculate_Multilayer_equation(double complex  (*A)[2][2], double complex  (*B)[2][2], double complex (*X)[2][2], double complex  (*result)[2][2], int N)

cdef void Calculate_Multilayer(double complex *t_comp1_up, double complex *t_comp2_up, double complex *t_comp1_do, double complex *t_comp2_do, double complex *r_ML_in1, double complex *r_ML_in2, double complex *r_ML_ba1, double complex *r_ML_ba2, int N)

cdef void Calculate_Multilayer_with_Matrices(double complex (*t_comp1_up)[2][2], double complex (*t_comp2_up)[2][2], double complex (*t_comp1_do)[2][2], double complex (*t_comp2_do)[2][2], double complex (*r_ML_in1)[2][2], double complex (*r_ML_in2)[2][2], double complex (*r_ML_ba1)[2][2], double complex (*r_ML_ba2)[2][2], int N)

cdef void Matrixexp( double complex (*A)[4][4], int N)
