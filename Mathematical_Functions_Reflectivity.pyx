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


cdef double Cmaxnorm(double complex x):
    if( x.imag < 0 ):
        x.imag=-x.imag
    if( x.real < 0 ):
        x.real=-x.real
    if( x.imag>x.real):
        return x.imag
    else:
        return x.real

cdef double complex ipow(double complex base, int exp):
    cdef double complex result = 1.
    while (exp):
        if (exp%2): #If exp is Uneven
            result *= base
        exp=exp/2 #exp halbieren
        base=base*base
    return result




cdef void Mult2x2_rightside(double complex (*A)[2][2], double complex (*B)[2][2]):

    cdef double complex R11,R12,R21,R22
    R11=(A[0])[0][0]*(B[0])[0][0]+(A[0])[0][1]*(B[0])[1][0]
    R12=(A[0])[0][0]*(B[0])[0][1]+(A[0])[0][1]*(B[0])[1][1]
    R21=(A[0])[1][0]*(B[0])[0][0]+(A[0])[1][1]*(B[0])[1][0]
    R22=(A[0])[1][0]*(B[0])[0][1]+(A[0])[1][1]*(B[0])[1][1]
    (B[0])[0][0]=R11
    (B[0])[0][1]=R12
    (B[0])[1][0]=R21
    (B[0])[1][1]=R22

cdef void Mult2x2_leftside(double complex (*A)[2][2], double complex (*B)[2][2]):

    cdef double complex R11,R12,R21,R22
    R11=(A[0])[0][0]*(B[0])[0][0]+(A[0])[0][1]*(B[0])[1][0]
    R12=(A[0])[0][0]*(B[0])[0][1]+(A[0])[0][1]*(B[0])[1][1]
    R21=(A[0])[1][0]*(B[0])[0][0]+(A[0])[1][1]*(B[0])[1][0]
    R22=(A[0])[1][0]*(B[0])[0][1]+(A[0])[1][1]*(B[0])[1][1]
    (A[0])[0][0]=R11
    (A[0])[0][1]=R12
    (A[0])[1][0]=R21
    (A[0])[1][1]=R22

cdef void Elimination_4x4(double complex (*A)[4][4], double complex (*B)[4]):
     # Calculates for matrix a and vector b the solution of the system of equations ax=b and stores the result in b.
    cdef double complex x, sum
    cdef int n=4
    cdef int i,j,k
    for i in range(n-1):
        for k in range(i+1,n):
            x=(A[0])[k][i]/(A[0])[i][i]
            for j in range(i+1,n):
                (A[0])[k][j] = (A[0])[k][j] -(A[0])[i][j]*x
            (B[0])[k] = (B[0])[k] - (B[0])[i]*x

    (B[0])[n-1]=(B[0])[n-1]/(A[0])[n-1][n-1]

    k=n-2
    while k>=0:
        sum = (B[0])[k]
        for j in range(k+1,n):
            sum=sum-(A[0])[k][j]*(B[0])[j]

        (B[0])[k]=sum/(A[0])[k][k]
        k=k-1



cdef void Invert2x2(double complex (*M)[2][2]):
    cdef double complex dconj=cconj( (M[0])[0][0]*(M[0])[1][1]-(M[0])[0][1]*(M[0])[1][0] )
    cdef double complex safe=(M[0])[0][0]
    dconj=dconj/cabsquadr(dconj)
    (M[0])[0][0]=(M[0])[1][1]*dconj
    (M[0])[0][1]=-(M[0])[0][1]*dconj
    (M[0])[1][0]=-(M[0])[1][0]*dconj
    (M[0])[1][1]=safe*dconj


cdef double twobytwocondi(double complex (*Matr)[2][2] ):
    ## Condition number of a 2x2 matrix in Frobenius norm
    cdef double sum1, sum2
    cdef double complex det
    sum1=  cabsquadr( (Matr[0])[0][0] ) \
         + cabsquadr( (Matr[0])[0][1] ) \
         + cabsquadr( (Matr[0])[1][0] ) \
         + cabsquadr( (Matr[0])[1][1] )
    sum1=dsqrt(sum1)

    det= (Matr[0])[0][0]*(Matr[0])[1][1]-(Matr[0])[0][1]*(Matr[0])[1][0]

    if( det == 0):
        det=2.0e-16

    sum2=  cabsquadr( (Matr[0])[1][1]/det ) \
         + cabsquadr( (Matr[0])[0][1]/det  ) \
         + cabsquadr( (Matr[0])[1][0]/det  ) \
         + cabsquadr( (Matr[0])[0][0]/det  )
    sum2=dsqrt(sum2)

    return sum1*sum2



cdef void Elimination2x2( double complex (*Matr)[2][2], double complex (*x)[2], double complex (*b)[2]  ):

    det= (Matr[0])[0][0]*(Matr[0])[1][1]-(Matr[0])[0][1]*(Matr[0])[1][0]

    if( det == 0):
        det=2.0e-16
    cdef double c1, c2
    c1=Cmaxnorm( (Matr[0])[0][0] )
    c2=Cmaxnorm( (Matr[0])[1][1] )
    if( c1==0 and c2==0 ):

        if( (Matr[0])[1][0]==0  ):
            (x[0])[0]=0
        else:
            (x[0])[0]=(b[0])[1]/(Matr[0])[1][0]
        if( (Matr[0])[0][1]==0 ):
            (x[0])[1]=0
        else:
            (x[0])[1]=(b[0])[0]/(Matr[0])[0][1]
    elif( c1>c2 ):

        (x[0])[1]=( (Matr[0])[0][0]*(b[0])[1]-(Matr[0])[1][0]*(b[0])[0] )/det
        (x[0])[0]=( (b[0])[1]-(Matr[0])[0][1]*(x[0])[1])/(Matr[0])[0][0]
    else:

        (x[0])[0]=( (Matr[0])[1][1]*(b[0])[0]-(Matr[0])[0][1]*(b[0])[1] )/det
        (x[0])[1]=( (b[0])[1]-(Matr[0])[1][0]*(x[0])[0])/(Matr[0])[1][1]

  #  return x

cdef void Normalize4vector( double complex (*vec)[4] ):
    cdef double sum1

    sum1=0
    for i in range(4 ):
        sum1+=cabsquadr(  (vec[0])[i] )
    sum1=dsqrt(sum1)
    for i in range(4):
        (vec[0])[i]/=sum1



cdef void Matrix4Invert(double complex (*a)[4][4], double complex (*invmat)[4][4] ):
    cdef double complex Dtot


    cdef double complex safe1[6]
    cdef double complex divide[4]

    ###Special case that occurs often:
    if(   (a[0])[2][0]==0 and ( (a[0])[3][0]==0 ) and (a[0])[2][2]==0 and ( (a[0])[3][2]==0 ) \
      and (a[0])[0][1]==0 and ( (a[0])[1][1]==0 ) and (a[0])[0][3]==0 and ( (a[0])[1][3]==0 )   ):


        Dtot= (a[0])[0][2]*(a[0])[1][0]-(a[0])[0][0]*(a[0])[1][2]
        (invmat[0])[0][0] = -(a[0])[1][2]/Dtot
        (invmat[0])[0][1] =  (a[0])[0][2]/Dtot
        (invmat[0])[0][2] = 0
        (invmat[0])[0][3] = 0

        (invmat[0])[2][0] =  (a[0])[1][0]/Dtot
        (invmat[0])[2][1] = -(a[0])[0][0]/Dtot
        (invmat[0])[2][2] = 0
        (invmat[0])[2][3] = 0


        Dtot= (a[0])[2][1]*(a[0])[3][3]-(a[0])[2][3]*(a[0])[3][1]
        (invmat[0])[1][0] = 0
        (invmat[0])[1][1] = 0
        (invmat[0])[1][2] =   (a[0])[3][3]/Dtot
        (invmat[0])[1][3] =  -(a[0])[2][3]/Dtot
        (invmat[0])[3][0] = 0
        (invmat[0])[3][1] = 0
        (invmat[0])[3][2] =  -(a[0])[3][1]/Dtot
        (invmat[0])[3][3] =   (a[0])[2][1]/Dtot
    else:


        safe1[0]=(a[0])[2][2]*(a[0])[3][3]-(a[0])[2][3]*(a[0])[3][2]
        safe1[1]=(a[0])[2][1]*(a[0])[3][3]-(a[0])[2][3]*(a[0])[3][1]
        safe1[2]=(a[0])[2][1]*(a[0])[3][2]-(a[0])[2][2]*(a[0])[3][1]
        safe1[3]=(a[0])[2][0]*(a[0])[3][3]-(a[0])[2][3]*(a[0])[3][0]
        safe1[4]=(a[0])[2][0]*(a[0])[3][2]-(a[0])[2][2]*(a[0])[3][0]
        safe1[5]=(a[0])[2][0]*(a[0])[3][1]-(a[0])[2][1]*(a[0])[3][0]

        divide[0]= (a[0])[1][1]*safe1[0]-(a[0])[1][2]*safe1[1]+(a[0])[1][3]*safe1[2]
        divide[1]=-(a[0])[1][0]*safe1[0]+(a[0])[1][2]*safe1[3]-(a[0])[1][3]*safe1[4]
        divide[2]= (a[0])[1][0]*safe1[1]-(a[0])[1][1]*safe1[3]+(a[0])[1][3]*safe1[5]
        divide[3]=-(a[0])[1][0]*safe1[2]+(a[0])[1][1]*safe1[4]-(a[0])[1][2]*safe1[5]

        Dtot=(a[0])[0][0]*divide[0]+(a[0])[0][1]*divide[1]+(a[0])[0][2]*divide[2]+(a[0])[0][3]*divide[3]

        (invmat[0])[0][0]=divide[0]/Dtot
        (invmat[0])[1][0]=divide[1]/Dtot
        (invmat[0])[2][0]=divide[2]/Dtot
        (invmat[0])[3][0]=divide[3]/Dtot

        (invmat[0])[0][1]=(-(a[0])[0][1]*safe1[0]+(a[0])[0][2]*safe1[1]-(a[0])[0][3]*safe1[2])/Dtot
        (invmat[0])[1][1]=( (a[0])[0][0]*safe1[0]-(a[0])[0][2]*safe1[3]+(a[0])[0][3]*safe1[4])/Dtot
        (invmat[0])[2][1]=(-(a[0])[0][0]*safe1[1]+(a[0])[0][1]*safe1[3]-(a[0])[0][3]*safe1[5])/Dtot
        (invmat[0])[3][1]=( (a[0])[0][0]*safe1[2]-(a[0])[0][1]*safe1[4]+(a[0])[0][2]*safe1[5])/Dtot

        safe1[0]=(a[0])[1][2]*(a[0])[3][3]-(a[0])[1][3]*(a[0])[3][2]
        safe1[1]=(a[0])[1][1]*(a[0])[3][3]-(a[0])[1][3]*(a[0])[3][1]
        safe1[2]=(a[0])[1][1]*(a[0])[3][2]-(a[0])[1][2]*(a[0])[3][1]
        safe1[3]=(a[0])[1][0]*(a[0])[3][3]-(a[0])[1][3]*(a[0])[3][0]
        safe1[4]=(a[0])[1][0]*(a[0])[3][2]-(a[0])[1][2]*(a[0])[3][0]
        safe1[5]=(a[0])[1][0]*(a[0])[3][1]-(a[0])[1][1]*(a[0])[3][0]

        (invmat[0])[0][2]=( (a[0])[0][1]*safe1[0]-(a[0])[0][2]*safe1[1]+(a[0])[0][3]*safe1[2])/Dtot
        (invmat[0])[1][2]=(-(a[0])[0][0]*safe1[0]+(a[0])[0][2]*safe1[3]-(a[0])[0][3]*safe1[4])/Dtot
        (invmat[0])[2][2]=( (a[0])[0][0]*safe1[1]-(a[0])[0][1]*safe1[3]+(a[0])[0][3]*safe1[5])/Dtot
        (invmat[0])[3][2]=(-(a[0])[0][0]*safe1[2]+(a[0])[0][1]*safe1[4]-(a[0])[0][2]*safe1[5])/Dtot

        safe1[0]=(a[0])[1][2]*(a[0])[2][3]-(a[0])[1][3]*(a[0])[2][2]
        safe1[1]=(a[0])[1][1]*(a[0])[2][3]-(a[0])[1][3]*(a[0])[2][1]
        safe1[2]=(a[0])[1][1]*(a[0])[2][2]-(a[0])[1][2]*(a[0])[2][1]
        safe1[3]=(a[0])[1][0]*(a[0])[2][3]-(a[0])[1][3]*(a[0])[2][0]
        safe1[4]=(a[0])[1][0]*(a[0])[2][2]-(a[0])[1][2]*(a[0])[2][0]
        safe1[5]=(a[0])[1][0]*(a[0])[2][1]-(a[0])[1][1]*(a[0])[2][0]

        (invmat[0])[0][3]=(-(a[0])[0][1]*safe1[0]+(a[0])[0][2]*safe1[1]-(a[0])[0][3]*safe1[2])/Dtot
        (invmat[0])[1][3]=(+(a[0])[0][0]*safe1[0]-(a[0])[0][2]*safe1[3]+(a[0])[0][3]*safe1[4])/Dtot
        (invmat[0])[2][3]=(-(a[0])[0][0]*safe1[1]+(a[0])[0][1]*safe1[3]-(a[0])[0][3]*safe1[5])/Dtot
        (invmat[0])[3][3]=(+(a[0])[0][0]*safe1[2]-(a[0])[0][1]*safe1[4]+(a[0])[0][2]*safe1[5])/Dtot


cdef void Mult4x4_leftside(double complex (*A)[4][4], double complex (*B)[4][4]):

    cdef double complex R[4][4]
    cdef int i, j, k
    for i in range(4):
        for j in range(4):
            R[i][j]=0
            for k in range(4):
                R[i][j]+=(A[0])[i][k]*(B[0])[k][j]

    for i in range(4):
        for j in range(4):
            (A[0])[i][j]=R[i][j]

cdef void Mult4x4_leftside_diag(double complex (*A)[4][4], double complex (*B)[4]):

    cdef int i, j
    for i in range(4):
        for j in range(4):
            (A[0])[i][j]*=(B[0])[j]

cdef double dabsvalue(double x) nogil:
    if(x<0):
        return -x
    else:
        return x

cdef double complex Maxcomplex(double complex a, double complex b):
    if(Cmaxnorm(a)>Cmaxnorm(b) ):
        return a
    else:
        return b


