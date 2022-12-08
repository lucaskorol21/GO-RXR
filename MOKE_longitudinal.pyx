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


cdef void Fill_rMemory(rMemory *Mem, double vy, double vyvy, double vy4, double complex chix, double complex chiy, double complex chiz, double complex chig):

    #works!
    cdef double complex B1mag, B1nonmag, xld, xldsquared, magrest, diftester, diftester2

    (Mem[0]).cgcg=chig*chig
    (Mem[0]).epsy=chiy+1.0
    (Mem[0]).epsz=chiz+1.0
    xld=(chix-chiy) + vyvy*(chiy-chiz)/(Mem[0]).epsz
    xldsquared = cquadr( xld )
    (Mem[0]).C1=chix+chiy+chix*chiy
    (Mem[0]).C2=(1.0+1.0/(Mem[0]).epsz)
    (Mem[0]).Delta31=-vy*chig/(Mem[0]).epsz
    B1mag=-(Mem[0]).cgcg/(Mem[0]).epsz
    B1nonmag=-chiy-chix+vyvy*chiy/(Mem[0]).epsz
    (Mem[0]).B1=B1mag+B1nonmag
    (Mem[0]).B2=-2.0+vyvy*(Mem[0]).C2
    magrest = B1mag*(2*B1nonmag + B1mag + 2.0*(  vyvy*(Mem[0]).C2 + 2*chiy   ) )
    (Mem[0]).root=sqrt( magrest + xldsquared )


    (Mem[0]).B=(Mem[0]).B1+(Mem[0]).B2
    (Mem[0]).vz1=sqrt((-(Mem[0]).B-(Mem[0]).root)/2)
    (Mem[0]).vz2=sqrt((-(Mem[0]).B+(Mem[0]).root)/2)
    diftester = xld-(Mem[0]).root
    diftester2= xld+(Mem[0]).root
    ##Choose the equation that has the higher numerical stability
    if( Cmaxnorm(diftester) < Cmaxnorm( diftester2 ) ):

        (Mem[0]).PHI1=(  diftester2+(Mem[0]).cgcg/(Mem[0]).epsz)/( 2*(Mem[0]).Delta31)
    else:

        (Mem[0]).PHI1=( -magrest/diftester +  (Mem[0]).cgcg/(Mem[0]).epsz  )/( 2*(Mem[0]).Delta31)

    (Mem[0]).PHI2=(Mem[0]).PHI1/(Mem[0]).epsy
    (Mem[0]).IsFilled=1


cdef void Calculate_rt_y(rMemory *Mem1, rMemory *Mem2, double vy, double vyvy, double vy4, double complex chix1, double complex chiy1, double complex chiz1, double complex chig1, double complex chix2, double complex chiy2, double complex chiz2, double complex chig2, \
                    int IsMagnetic1, int IsMagnetic2, double complex (*r)[2][2], double complex (*rprime)[2][2], double complex (*t)[2][2], double complex (*tprime)[2][2], double sigma, double k0):
    #seems to work


    cdef double complex vdifterm
    cdef double complex J[2][4]
    cdef double complex roughmat[2][4]
    cdef double complex sumv11v12
    cdef double complex sumv11v22
    cdef double complex sumv12v21
    cdef double complex sumv21v22
    cdef double complex difvz1vz3
    cdef double complex difvz1vz4
    cdef double complex difvz2vz3
    cdef double complex difvz2vz4

    cdef double complex divide1
    cdef double complex divide2
    cdef double complex fac1
    cdef double roughfac=-0.5*quadr(sigma*  k0)
    cdef double complex v1diff, v2diff
    cdef double complex J11030113, J00121002, J02010300, J13101211

    if( IsMagnetic1 & IsMagnetic2 ):
        Fill_rMemory(Mem2, vy, vyvy, vy4, chix2, chiy2, chiz2, chig2)
    #    Ndiv=2*((Mem1[0]).PHI1*(Mem1[0]).PHI2-1.0)
        vdifterm=vyvy*(chiz1-chiz2)/(  (Mem1[0]).epsz*(Mem2[0]).epsz )
        sumv11v12=((Mem1[0]).vz1+(Mem2[0]).vz1)
        sumv11v22=((Mem1[0]).vz1+(Mem2[0]).vz2)
        sumv12v21=((Mem1[0]).vz2+(Mem2[0]).vz1)
        sumv21v22=((Mem1[0]).vz2+(Mem2[0]).vz2)
        difvz1vz3=0.5*( (Mem2[0]).B1+(Mem2[0]).root-(Mem1[0]).B1-(Mem1[0]).root+vdifterm )/sumv11v12
        difvz1vz4=0.5*( (Mem2[0]).B1-(Mem2[0]).root-(Mem1[0]).B1-(Mem1[0]).root+vdifterm )/sumv11v22
        difvz2vz3=0.5*( (Mem2[0]).B1+(Mem2[0]).root-(Mem1[0]).B1+(Mem1[0]).root+vdifterm )/sumv12v21
        difvz2vz4=0.5*( (Mem2[0]).B1-(Mem2[0]).root-(Mem1[0]).B1+(Mem1[0]).root+vdifterm )/sumv21v22
        v1diff = (Mem1[0]).root/( (Mem1[0]).vz1 + (Mem1[0]).vz2 )
        v2diff = (Mem2[0]).root/( (Mem2[0]).vz1 + (Mem2[0]).vz2 )
#        divide1=1.0/(Ndiv*(Mem1[0]).vz1)
#        divide2=1.0/(Ndiv*(Mem1[0]).vz2)
        divide1=1.0/((Mem1[0]).vz1)
        divide2=1.0/((Mem1[0]).vz2)
        fac1=((Mem2[0]).PHI1*(Mem1[0]).PHI2/(Mem2[0]).epsy)

        J[0][0]=divide1*( fac1*(sumv11v12+chiy2*(Mem1[0]).vz1+chiy1*(Mem2[0]).vz1)-sumv11v12 )
        J[0][2]=divide1*( fac1*(difvz1vz3+chiy2*(Mem1[0]).vz1-chiy1*(Mem2[0]).vz1)-difvz1vz3 )

        fac1=(-(Mem2[0]).PHI1/(Mem2[0]).epsy)
        J[1][0]=divide2*( fac1*(sumv12v21+chiy2*(Mem1[0]).vz2+chiy1*(Mem2[0]).vz1)+(Mem1[0]).PHI1*sumv12v21 )
        J[1][2]=divide2*( fac1*(difvz2vz3+chiy2*(Mem1[0]).vz2-chiy1*(Mem2[0]).vz1)+(Mem1[0]).PHI1*difvz2vz3 )

        fac1=( (Mem1[0]).PHI2 )/(Mem2[0]).epsy
        J[0][1]=divide1*( fac1*(sumv11v22+chiy2*(Mem1[0]).vz1+chiy1*(Mem2[0]).vz2)-(Mem2[0]).PHI2*sumv11v22 )
        J[0][3]=divide1*( fac1*(difvz1vz4+chiy2*(Mem1[0]).vz1-chiy1*(Mem2[0]).vz2)-(Mem2[0]).PHI2*difvz1vz4 )

        fac1=(Mem2[0]).PHI2*(Mem1[0]).PHI1
        J[1][1]=divide2*(fac1*(sumv21v22)-(sumv21v22+chiy2*(Mem1[0]).vz2+chiy1*(Mem2[0]).vz2)/(Mem2[0]).epsy )
        J[1][3]=divide2*(fac1*(difvz2vz4)-(difvz2vz4+chiy2*(Mem1[0]).vz2-chiy1*(Mem2[0]).vz2)/(Mem2[0]).epsy )
#        print("facs")
#        print( ((Mem2[0]).PHI1*(Mem1[0]).PHI2/(Mem2[0]).epsy) )
#        print( (-(Mem2[0]).PHI1/(Mem2[0]).epsy) )
#        print( ( (Mem1[0]).PHI2 )/(Mem2[0]).epsy )
#        print( (Mem2[0]).PHI2*(Mem1[0]).PHI1 )


      #  print( "old" ,(J[1][1]*J[0][3]-J[0][1]*J[1][3]) )

#        J11030113 = (J[1][1]*J[0][3]-J[0][1]*J[1][3])
#        J00121002 = (J[0][0]*J[1][2]-J[1][0]*J[0][2])
#        J02010300 = (J[0][2]*J[0][1]-J[0][3]*J[0][0])
#        J13101211 = (J[1][3]*J[1][0]-J[1][2]*J[1][1])


#        J11030113 = divide1*divide2*(\
#    + 2*(Mem1[0]).PHI1*(Mem1[0]).PHI2*(Mem1[0]).vz1*(Mem2[0]).PHI2*(Mem2[0]).vz2*chiy2/(Mem2[0]).epsy \
#    + 2*(Mem1[0]).PHI1*(Mem1[0]).PHI2*(Mem1[0]).vz1*(Mem2[0]).PHI2*(Mem2[0]).vz2/(Mem2[0]).epsy \
#    - 2*(Mem1[0]).PHI1*(Mem1[0]).PHI2*(Mem1[0]).vz2*(Mem2[0]).PHI2*(Mem2[0]).vz2*chiy1/(Mem2[0]).epsy \
#    - 2*(Mem1[0]).PHI1*(Mem1[0]).PHI2*(Mem1[0]).vz2*(Mem2[0]).PHI2*(Mem2[0]).vz2/(Mem2[0]).epsy \
#    - 2*(Mem1[0]).PHI1*(Mem1[0]).vz1*(Mem2[0]).PHI2**2*(Mem2[0]).vz2 \
#    + 2*(Mem1[0]).PHI1*(Mem1[0]).vz2*(Mem2[0]).PHI2**2*(Mem2[0]).vz2 \
#    - 2*(Mem1[0]).PHI2*(Mem1[0]).vz1*(Mem2[0]).vz2*chiy1*chiy2/((Mem2[0]).epsy*(Mem2[0]).epsy)\
#    - 2*(Mem1[0]).PHI2*(Mem1[0]).vz1*(Mem2[0]).vz2*chiy1/((Mem2[0]).epsy*(Mem2[0]).epsy)\
#    - 2*(Mem1[0]).PHI2*(Mem1[0]).vz1*(Mem2[0]).vz2*chiy2/((Mem2[0]).epsy*(Mem2[0]).epsy)\
#    - 2*(Mem1[0]).PHI2*(Mem1[0]).vz1*(Mem2[0]).vz2/((Mem2[0]).epsy*(Mem2[0]).epsy)\
#    + 2*(Mem1[0]).PHI2*(Mem1[0]).vz2*(Mem2[0]).vz2*chiy1*chiy2/((Mem2[0]).epsy*(Mem2[0]).epsy)\
#    + 2*(Mem1[0]).PHI2*(Mem1[0]).vz2*(Mem2[0]).vz2*chiy1/((Mem2[0]).epsy*(Mem2[0]).epsy)\
#    + 2*(Mem1[0]).PHI2*(Mem1[0]).vz2*(Mem2[0]).vz2*chiy2/((Mem2[0]).epsy*(Mem2[0]).epsy)\
#    + 2*(Mem1[0]).PHI2*(Mem1[0]).vz2*(Mem2[0]).vz2/((Mem2[0]).epsy*(Mem2[0]).epsy)\
#    + 2*(Mem1[0]).vz1*(Mem2[0]).PHI2*(Mem2[0]).vz2*chiy1/(Mem2[0]).epsy \
#    + 2*(Mem1[0]).vz1*(Mem2[0]).PHI2*(Mem2[0]).vz2/(Mem2[0]).epsy \
#    - 2*(Mem1[0]).vz2*(Mem2[0]).PHI2*(Mem2[0]).vz2*chiy2/(Mem2[0]).epsy \
#    - 2*(Mem1[0]).vz2*(Mem2[0]).PHI2*(Mem2[0]).vz2/(Mem2[0]).epsy \
#        )
#        J11030113 = divide1*divide2*2*(Mem2[0]).vz2*( \
#       (Mem2[0]).PHI2*(1.0/(Mem2[0]).epsy)*0.5*( (Mem1[0]).vz1 + (Mem1[0]).vz2)*(chiy2-chiy1)*( (Mem1[0]).PHI1*(Mem1[0]).PHI2 - 1 )\
#                 +v1diff*( +( (Mem1[0]).PHI2/((Mem2[0]).epsy) )*( (Mem1[0]).epsy - (Mem1[0]).PHI1*(Mem2[0]).PHI2)\
#                    +  (Mem2[0]).PHI2*( (Mem1[0]).PHI1*(Mem2[0]).PHI2 - 1.0/(Mem2[0]).epsy -0.5*(chiy1+chiy2)*( 1+(Mem1[0]).PHI1*(Mem1[0]).PHI2  )/(Mem2[0]).epsy ) ) )
#        print( "new", J11030113)


    elif( IsMagnetic1 ):
        (Mem2[0]).IsFilled=1
     #   Ndiv=2*((Mem1[0]).PHI1*(Mem1[0]).PHI2-1.0)
        (Mem2[0]).epsy=1.0+chiy2
        (Mem2[0]).epsz=1.0+chiz2
        (Mem2[0]).vz1=sqrt(1.+chix2-vyvy)
        (Mem2[0]).vz2=sqrt((1.-vyvy/(Mem2[0]).epsz)*(Mem2[0]).epsy)
        sumv11v12=((Mem1[0]).vz1+(Mem2[0]).vz1)
        sumv11v22=((Mem1[0]).vz1+(Mem2[0]).vz2)
        sumv12v21=((Mem1[0]).vz2+(Mem2[0]).vz1)
        sumv21v22=((Mem1[0]).vz2+(Mem2[0]).vz2)
       # vz1diff = (Mem1[0]).root/( (Mem1[0]).vz1 + (Mem1[0]).vz2 )

        #vdifterm=vyvy*(chiz1-chiz2)/(  (Mem1[0]).epsz*(Mem2[0]).epsz )
        difvz1vz3=( 0.5*(-(Mem1[0]).B1-(Mem1[0]).root)-chix2+0.5*vyvy*chiz1/(Mem1[0]).epsz )/sumv11v12
        difvz1vz4=( 0.5*(-(Mem1[0]).B1-(Mem1[0]).root)-(vyvy*(chiz2-chiz1)/((Mem2[0]).epsz*(Mem1[0]).epsz) +vyvy*chiz1/(2*(Mem1[0]).epsz) + (1-vyvy/(Mem2[0]).epsz)*chiy2)   )/sumv11v22
        difvz2vz3=( 0.5*(-(Mem1[0]).B1+(Mem1[0]).root)-chix2+0.5*vyvy*chiz1/(Mem1[0]).epsz )/sumv12v21
        difvz2vz4=( 0.5*(-(Mem1[0]).B1+(Mem1[0]).root)-(vyvy*(chiz2-chiz1)/((Mem2[0]).epsz*(Mem1[0]).epsz) +vyvy*chiz1/(2*(Mem1[0]).epsz) + (1-vyvy/(Mem2[0]).epsz)*chiy2)   )/sumv21v22

        divide1=1.0/((Mem1[0]).vz1)
        divide2=1.0/((Mem1[0]).vz2)
#        divide1=1.0/(Ndiv*(Mem1[0]).vz1)
#        divide2=1.0/(Ndiv*(Mem1[0]).vz2)
        J[0][0]=-sumv11v12*divide1
        J[0][2]=-difvz1vz3*divide1

        J[1][0]=divide2*(Mem1[0]).PHI1*sumv12v21
        J[1][2]=divide2*(Mem1[0]).PHI1*difvz2vz3

        fac1=( (Mem1[0]).PHI2 )/(Mem2[0]).epsy
        J[0][1]=divide1*( fac1*(sumv11v22+chiy2*(Mem1[0]).vz1+chiy1*(Mem2[0]).vz2) )
        J[0][3]=divide1*( fac1*(difvz1vz4+chiy2*(Mem1[0]).vz1-chiy1*(Mem2[0]).vz2) )


        J[1][1]=divide2*(-(sumv21v22+chiy2*(Mem1[0]).vz2+chiy1*(Mem2[0]).vz2)/(Mem2[0]).epsy )
        J[1][3]=divide2*(-(difvz2vz4+chiy2*(Mem1[0]).vz2-chiy1*(Mem2[0]).vz2)/(Mem2[0]).epsy )
#        v1diff = (Mem1[0]).root/( (Mem1[0]).vz1 + (Mem1[0]).vz2 )
#        v2diff = ( (chiy2-chix2) + vyvy*(chiz2-chiy2)/(Mem2[0]).epsz )/( (Mem2[0]).vz1 + (Mem2[0]).vz2 )

#        J11030113 = divide1*divide2*fac1*( \
#             +   ( (Mem1[0]).epsy+chiy2*(1.+1.*chiy1)  )*2*(Mem2[0]).vz2*v1diff  \
#             +   chiy2*Mem1[0].vz1*difvz2vz4    )


#        J11030113 = divide1*divide2*fac1*2*(Mem2[0]).vz2*v1diff*(Mem1[0]).epsy
#
#
#
#      #  print( "old" ,(J[0][0]*J[1][2]-J[1][0]*J[0][2]) )
#        J00121002 =  -divide1*divide2*(Mem1[0]).PHI1*2*(Mem2[0]).vz1*v1diff
#      #  print( "old2" ,J00121002 )
#     #   print( "old" ,(J[0][2]*J[0][1]-J[0][3]*J[0][0]) )
#        J02010300 = divide1*divide1*fac1*( \
#        2*(Mem1[0]).vz1*( chiy2*(Mem2[0]).vz1-chiy1*(Mem2[0]).vz2  - v2diff ) )
#      #  print( "old2" ,J02010300 )
#     #   print( "old" ,(J[1][3]*J[1][0]-J[1][2]*J[1][1]) )
#        J13101211 = divide2*divide2*( (Mem1[0]).PHI1/(Mem2[0]).epsy )*( \
#              +2*(Mem1[0]).vz2*(  chiy1*(Mem2[0]).vz2 - chiy2*(Mem2[0]).vz1  +v2diff) )




       # print("old2",J13101211 )


       # print(J13101211)

    elif( IsMagnetic2 ):
        (Mem2[0]).IsFilled=1
        Fill_rMemory(Mem2, vy, vyvy, vy4, chix2, chiy2, chiz2, chig2)
        sumv11v12=((Mem1[0]).vz1+(Mem2[0]).vz1)
        sumv11v22=((Mem1[0]).vz1+(Mem2[0]).vz2)
        sumv12v21=((Mem1[0]).vz2+(Mem2[0]).vz1)
        sumv21v22=((Mem1[0]).vz2+(Mem2[0]).vz2)
        difvz1vz3=( 0.5*((Mem2[0]).B1+(Mem2[0]).root)+chix1-0.5*vyvy*chiz2/(Mem2[0]).epsz )/sumv11v12
        difvz1vz4=( 0.5*((Mem2[0]).B1-(Mem2[0]).root)+chix1-0.5*vyvy*chiz2/(Mem2[0]).epsz )/sumv11v22
        difvz2vz3=( 0.5*((Mem2[0]).B1+(Mem2[0]).root)+vyvy*(chiz1-chiz2)/((Mem1[0]).epsz*(Mem2[0]).epsz) +vyvy*chiz2/(2*(Mem2[0]).epsz) + (1-vyvy/(Mem1[0]).epsz)*chiy1   )/sumv12v21
        difvz2vz4=( 0.5*((Mem2[0]).B1-(Mem2[0]).root)+vyvy*(chiz1-chiz2)/((Mem1[0]).epsz*(Mem2[0]).epsz) +vyvy*chiz2/(2*(Mem2[0]).epsz) + (1-vyvy/(Mem1[0]).epsz)*chiy1   )/sumv21v22
        divide1=0.5/(Mem1[0]).vz1
        divide2=0.5/( (Mem1[0]).vz2*(Mem2[0]).epsy )

        J[0][0]=divide1*sumv11v12
        J[0][2]=divide1*difvz1vz3

        J[1][0]=divide2*(Mem2[0]).PHI1*(sumv12v21+chiy2*(Mem1[0]).vz2+chiy1*(Mem2[0]).vz1)
        J[1][2]=divide2*(Mem2[0]).PHI1*(difvz2vz3+chiy2*(Mem1[0]).vz2-chiy1*(Mem2[0]).vz1)

        J[0][1]=divide1*(Mem2[0]).PHI2*sumv11v22
        J[0][3]=divide1*(Mem2[0]).PHI2*difvz1vz4

        J[1][1]=divide2*(sumv21v22+chiy2*(Mem1[0]).vz2+chiy1*(Mem2[0]).vz2)
        J[1][3]=divide2*(difvz2vz4+chiy2*(Mem1[0]).vz2-chiy1*(Mem2[0]).vz2)

#        v1diff = ( (chiy1-chix1) + vyvy*(chiz1-chiy1)/(Mem1[0]).epsz )/( (Mem1[0]).vz1 + (Mem1[0]).vz2 )
#        v2diff = (Mem2[0]).root/( (Mem2[0]).vz1 + (Mem2[0]).vz2 )
#
#      #  print( J[1][1]*J[0][3]-J[0][1]*J[1][3] )
#        J11030113 = divide1*divide2*(Mem2[0]).PHI2*( \
#                         +2*(Mem2[0]).vz2*(-v1diff + ( (Mem1[0]).vz1*chiy1-(Mem1[0]).vz2*chiy2 ) ) )
#       # print( J11030113)
#       # print( J[0][0]*J[1][2]-J[1][0]*J[0][2] )
#        J00121002 = divide1*divide2*(Mem2[0]).PHI1*( \
#                        2*(Mem2[0]).vz1*(v1diff + ( (Mem1[0]).vz2*chiy2 - (Mem1[0]).vz1*chiy1  ) )  )
#      #  print( J00121002)
#       # print( J[0][2]*J[0][1]-J[0][3]*J[0][0] )
#
#
#        J02010300 = divide1*divide1*(Mem2[0]).PHI2*( \
#           2*(Mem1[0]).vz1*v2diff   )
#      #  print( J02010300 )
#       # print( J[1][3]*J[1][0]-J[1][2]*J[1][1] )
#        J13101211 = -divide2*divide2*(Mem2[0]).PHI1*2*(Mem1[0]).vz2*v2diff*(Mem2[0]).epsy*(Mem1[0]).epsy

       # print( J13101211 )
    else:
        (Mem2[0]).IsFilled=1
        (Mem2[0]).epsy=1.0+chiy2
        (Mem2[0]).epsz=1.0+chiz2
        (Mem2[0]).vz1=sqrt(1.+chix2-vyvy)
        (Mem2[0]).vz2=sqrt((1.-vyvy/(Mem2[0]).epsz)*(Mem2[0]).epsy)
        sumv11v12=((Mem1[0]).vz1+(Mem2[0]).vz1)
        J[0][0]=sumv11v12/(Mem1[0]).vz1
        difvz1vz3=(chix1-chix2)/sumv11v12
        J[0][2]=difvz1vz3/(Mem1[0]).vz1
        J[1][1]=1+(Mem1[0]).epsy*(Mem2[0]).vz2/( (Mem2[0]).epsy*(Mem1[0]).vz2 )
        J[1][3]=( (Mem1[0]).epsz*(Mem2[0]).epsz*(chiy2-chiy1)+vyvy*(chiz1+chiy1+chiz1*chiy1-chiz2-chiy2-chiy2*chiz2) ) *(Mem1[0]).epsy/( (Mem1[0]).vz2*(Mem1[0]).epsz*(Mem2[0]).epsz*((Mem1[0]).vz2*(Mem2[0]).epsy+(Mem2[0]).vz2*(Mem1[0]).epsy)  )

        difvz2vz4=( -vyvy*(chiz2-chiz1)/( (Mem2[0]).epsz*(Mem1[0]).epsz ) + (1-vyvy/(Mem1[0]).epsz)*chiy1 - (1-vyvy/(Mem2[0]).epsz)*chiy2 )/( (Mem1[0]).vz2+(Mem2[0]).vz2 )

        J[1][3]=( 1/( (Mem2[0]).epsy*(Mem1[0]).vz2 ) )*( difvz2vz4 + chiy2*(Mem1[0]).vz2 -chiy1*(Mem2[0]).vz2 )

        J[0][0]*=exp(roughfac*cquadr(difvz1vz3) )
        J[0][2]*=exp(roughfac*cquadr((Mem1[0]).vz1+(Mem2[0]).vz1))
        J[1][1]*=exp(roughfac*cquadr(difvz2vz4) )
        J[1][3]*=exp(roughfac*cquadr((Mem1[0]).vz2+(Mem2[0]).vz2))
        (r[0])[0][0]=-J[0][2]/J[0][0]
        (r[0])[0][1]=0
        (r[0])[1][0]=0
        (r[0])[1][1]=-J[1][3]/J[1][1]
        (t[0])[0][0]=J[0][2]*(r[0])[0][0]+J[0][0] # Incoming 1 transmitted 1
        (t[0])[0][1]=0 # Incoming 2 transmitted 1
        (t[0])[1][0]=0 # Incoming 1 transmitted 2
        (t[0])[1][1]=J[1][3]*(r[0])[1][1]+J[1][1] # Incoming 2 transmitted 2
        (tprime[0])[0][0]=1.0/J[0][0]
        (tprime[0])[0][1]=0
        (tprime[0])[1][0]=0
        (tprime[0])[1][1]=1.0/J[1][1]
        (rprime[0])[0][0]=-(r[0])[0][0]
        (rprime[0])[0][1]=0
        (rprime[0])[1][0]=0
        (rprime[0])[1][1]=-(r[0])[1][1]
        return


    roughmat[0][0]=exp(roughfac*cquadr(difvz1vz3) )
    roughmat[0][2]=exp(roughfac*cquadr((Mem1[0]).vz1+(Mem2[0]).vz1))
    roughmat[1][1]=exp(roughfac*cquadr(difvz2vz4) )
    roughmat[1][3]=exp(roughfac*cquadr((Mem1[0]).vz2+(Mem2[0]).vz2))

    roughmat[0][1]=exp(roughfac*cquadr(difvz1vz4))
    roughmat[0][3]=exp(roughfac*cquadr((Mem1[0]).vz1+(Mem2[0]).vz2))
    roughmat[1][0]=exp(roughfac*cquadr(difvz2vz3) )
    roughmat[1][2]=exp(roughfac*cquadr((Mem1[0]).vz2+(Mem2[0]).vz1))
    J[0][0]*=roughmat[0][0]
    J[0][2]*=roughmat[0][2]
    J[1][1]*=roughmat[1][1]
    J[1][3]*=roughmat[1][3]
    J[0][1]*=roughmat[0][1]
    J[0][3]*=roughmat[0][3]
    J[1][0]*=roughmat[1][0]
    J[1][2]*=roughmat[1][2]

#    if(IsMagnetic1):
#        print( "before1", J[0][2]*J[0][1]-J[0][3]*J[0][0])
#        print( "before2", J02010300 )
#        print( "before3", J[1][3]*J[1][0]-J[1][2]*J[1][1])
#        print( "before4", J13101211 )

#    J11030113 = 0.5*( \
#     (J[1][1]*J[0][3]+J[0][1]*J[1][3])*( roughmat[1][1]*roughmat[0][3] - roughmat[0][1]*roughmat[1][3] ) \
#    +J11030113*( roughmat[1][1]*roughmat[0][3] + roughmat[0][1]*roughmat[1][3] ) )
#
#    J00121002 = 0.5*( \
#    (J[0][0]*J[1][2]+J[1][0]*J[0][2])*( \
#      ( roughmat[0][0]*roughmat[1][2] - roughmat[1][0]*roughmat[0][2] ) ) \
#    + J00121002*( \
#      ( roughmat[0][0]*roughmat[1][2] + roughmat[1][0]*roughmat[0][2] ) ) )
#
#    J02010300=0.5*( \
#     (J[0][2]*J[0][1]+J[0][3]*J[0][0])*(roughmat[0][2]*roughmat[0][1]-roughmat[0][3]*roughmat[0][0])  \
#                            +J02010300*(roughmat[0][2]*roughmat[0][1]+roughmat[0][3]*roughmat[0][0]) )
#
#    J13101211 = 0.5*( \
#     (J[1][3]*J[1][0]+J[1][2]*J[1][1])*( roughmat[1][3]*roughmat[1][0]-roughmat[1][2]*roughmat[1][1] ) \
#                            +J13101211*( roughmat[1][3]*roughmat[1][0]+roughmat[1][2]*roughmat[1][1] )
#    )





#    if(IsMagnetic1):
#        print( "after1", J[0][2]*J[0][1]-J[0][3]*J[0][0])
#        print( "after2", J02010300 )
#        print( "after3", J[1][3]*J[1][0]-J[1][2]*J[1][1])
#        print( "after4", J13101211 )
#    if(IsMagnetic1):
#        print("x")
#        print( "rough3", (J[0][0]*J[1][2]-J[1][0]*J[0][2]))
#        print( "rough4", J00121002)

    divide1=1.0/(J[0][1]*J[1][0]-J[1][1]*J[0][0])

    (r[0])[0][0]=(J[1][1]*J[0][2]-J[0][1]*J[1][2])*divide1 # Incoming 1 reflected 1
    (r[0])[0][1]=(J[1][1]*J[0][3]-J[0][1]*J[1][3])*divide1 # Incoming 2 reflected 1
    (r[0])[1][0]=(J[0][0]*J[1][2]-J[1][0]*J[0][2])*divide1 # Incoming 1 reflected 2
    (r[0])[1][1]=(J[0][0]*J[1][3]-J[1][0]*J[0][3])*divide1 # Incoming 2 reflected 2

#    (r[0])[0][0]=(J[1][1]*J[0][2]-J[0][1]*J[1][2])*divide1 # Incoming 1 reflected 1
#    (r[0])[0][1]=J11030113*divide1 # Incoming 2 reflected 1
#    (r[0])[1][0]=J00121002*divide1 # Incoming 1 reflected 2
#    (r[0])[1][1]=(J[0][0]*J[1][3]-J[1][0]*J[0][3])*divide1 # Incoming 2 reflected 2




#    if(IsMagnetic1):
#        print("1", (r[0])[0][1], (r[0])[1][0])
#    elif( IsMagnetic2):
#        print("2", (r[0])[0][1], (r[0])[1][0])

    (t[0])[0][0]=J[0][2]*(r[0])[0][0]+J[0][3]*(r[0])[1][0]+J[0][0] # Incoming 1 transmitted 1
    (t[0])[0][1]=J[0][2]*(r[0])[0][1]+J[0][3]*(r[0])[1][1]+J[0][1] # Incoming 2 transmitted 1
    (t[0])[1][0]=J[1][2]*(r[0])[0][0]+J[1][3]*(r[0])[1][0]+J[1][0] # Incoming 1 transmitted 2
    (t[0])[1][1]=J[1][2]*(r[0])[0][1]+J[1][3]*(r[0])[1][1]+J[1][1] # Incoming 2 transmitted 2

    (tprime[0])[0][0]=-J[1][1]*divide1
    (tprime[0])[0][1]=J[0][1]*divide1
    (tprime[0])[1][0]=J[1][0]*divide1
    (tprime[0])[1][1]=-J[0][0]*divide1



    (rprime[0])[0][0]=J[0][2]*((tprime[0])[0][0])+J[0][3]*((tprime[0])[1][0])
    (rprime[0])[0][1]=J[0][2]*((tprime[0])[0][1])+J[0][3]*((tprime[0])[1][1])
    (rprime[0])[1][0]=J[1][2]*((tprime[0])[0][0])+J[1][3]*((tprime[0])[1][0])
    (rprime[0])[1][1]=J[1][2]*((tprime[0])[0][1])+J[1][3]*((tprime[0])[1][1])





cdef void Paratt_magnetic_y(Heterostructure* HS, double th, double wavelength, double complex (*rtot)[2][2]):

    cdef double k0=6.283185307179586/wavelength
    cdef double sintheta=sin(two_pi_div_360()*th)

    cdef double vy=cos(two_pi_div_360()*th)
    cdef double vyvy=quadr(vy)
    cdef double vy4=quadr(vyvy)
    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
    cdef CLayer UpperLayer, LowerLayer
    cdef rMemory Memory1, Memory2
    cdef rMemory* Mempointer1
    cdef rMemory* Mempointer2
    cdef int ML_is_diagonal
    cdef double complex r[2][2]
    cdef double complex rprime[2][2]
    cdef double complex t[2][2]
    cdef double complex tprime[2][2]
    cdef double complex p[2][2]

    cdef double complex t_ML_in_1[2][2] #Transmission through one Compound layer, down
    cdef double complex t_ML_back_1[2][2] #Transmission through one Compound layer, up
    cdef double complex t_ML_in_2[2][2] #Transmission through all Compound layer, down
    cdef double complex t_ML_back_2[2][2] #Transmission through all Compound layer, up
    cdef double complex r_ML_in_1[2][2] # Upper reflectivity from one Compound layer
    cdef double complex r_ML_in_2[2][2] # Upper reflectivity from the Multilayer
    cdef double complex ev1, ev2, ev3, ev4
    cdef double product, product2, absevmax1, absevmax2
    cdef double normalizator
    p[0][1]=0
    p[1][0]=0
   # print "0"
  #  cdef double complex test=LR[MLCOMP[0][0]].ey

    cdef int Cap=NLAYERS-1

   # print "1"
    LowerLayer=LR[MLCOMP[0][0]]

    if(LowerLayer.magdir):

        Fill_rMemory(&Memory1, vy,vyvy,vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)
    else:
        Memory1.epsy=1.0+LowerLayer.cy
        Memory1.epsz=1.0+LowerLayer.cz
        Memory1.vz1=sqrt(1.+LowerLayer.cx-vyvy)
        Memory1.vz2=sqrt((1.-vyvy/Memory1.epsz)*Memory1.epsy)

    if(NLAYERS==1):
        Calculate_rt_y(&Memory1, &Memory2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0, 0, 0, 0, \
                    LowerLayer.magdir, 0, rtot, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
    else:

        UpperLayer=LR[MLCOMP[1][0]]
        Calculate_rt_y(&Memory1, &Memory2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, rtot, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
        Memory1.IsFilled=0
        p[0][0]=exp(1j*k0*UpperLayer.Thickness*Memory2.vz1)
        p[1][1]=exp(1j*k0*UpperLayer.Thickness*Memory2.vz2)
    i=1
    while i<NLAYERS:
       # print "loop start", i
        if(MLLENGTH[i]==1):

            #Check which Memory to use:
            if(Memory1.IsFilled):
                Mempointer1=&Memory1
                Mempointer2=&Memory2
            else:
                Mempointer2=&Memory1
                Mempointer1=&Memory2
            LowerLayer=LR[MLCOMP[i][0]]
            if(i!=Cap):
                #Upper=MLCOMP[i+1][0]
                UpperLayer=LR[MLCOMP[i+1][0]]
               # print "i, Upper", i, Upper
                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0
            else:
                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0,0,0,0, \
                    LowerLayer.magdir, 0, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)


            Mult2x2_rightside(&p, rtot)

            Mult2x2_leftside(rtot, &p)

            Mult2x2_rightside(&tprime, rtot)

            Mult2x2_leftside(rtot, &t)

            (rtot[0])[0][0]+=r[0][0]
            (rtot[0])[1][0]+=r[1][0]
            (rtot[0])[0][1]+=r[0][1]
            (rtot[0])[1][1]+=r[1][1]

            if(i!=Cap):
                p[0][0]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz1)
                p[1][1]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz2)




        else: #A Multilayer

            LowerLayer=LR[MLCOMP[i][0]]
            UpperLayer=LR[MLCOMP[i][1]]
            if(Memory1.IsFilled):
                Mempointer1=&Memory1
                Mempointer2=&Memory2
            else:
                Mempointer2=&Memory1
                Mempointer1=&Memory2

            Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, LowerLayer.Roughness, k0)
            (Mempointer1[0]).IsFilled=0
            if(LowerLayer.magdir):
                ML_is_diagonal=0

            Mult2x2_leftside(&t_ML_back_1, &p) # t'(AB)*p(A)
            Mult2x2_rightside(&p, &t_ML_in_1) # p(A) * t(AB)

            p[0][0]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz1)
            p[1][1]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz2)

            j=1
            while j<MLLENGTH[i]:
               # Upper=MLCOMP[i][(j+1)%MLLENGTH[i]]

                UpperLayer=LR[ MLCOMP[i][(j+1)%MLLENGTH[i]] ]
                LowerLayer=LR[ MLCOMP[i][j] ]
                if(LowerLayer.magdir):
                    ML_is_diagonal=0
                if(Memory1.IsFilled):
                    Mempointer1=&Memory1
                    Mempointer2=&Memory2
                else:
                    Mempointer2=&Memory1
                    Mempointer1=&Memory2

                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0

                Mult2x2_rightside(&p, &r_ML_in_1) #p_B r p_B and so on

                Mult2x2_leftside(&r_ML_in_1, &p)

                Mult2x2_rightside(&tprime, &r_ML_in_1) #t' p_B r p_b t and so on

                Mult2x2_leftside(&r_ML_in_1, &t)


                r_ML_in_1[0][0]+=r[0][0]
                r_ML_in_1[1][0]+=r[1][0]
                r_ML_in_1[0][1]+=r[0][1]
                r_ML_in_1[1][1]+=r[1][1]

                Mult2x2_leftside(&t_ML_in_1, &p)
                Mult2x2_rightside(&p, &t_ML_back_1)

                Mult2x2_leftside(&t_ML_in_1, &t) # p(A) * t(AB) p_B t(BC) and so on
                Mult2x2_rightside(&tprime, &t_ML_back_1)# t'(BC) p_B t'(AB)*p(A) and so on


                p[0][0]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz1)
                p[1][1]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz2)
           #     if(j==1):
            #        print "p C components"
            #        print Upper, LR[Upper].Thickness, vz3, vz4

                j=j+1






            if(ML_is_diagonal):
                r_ML_in_2[0][0]=r_ML_in_1[0][0]*(1-ipow(t_ML_in_1[0][0]*t_ML_back_1[0][0], MLREP[i]-1))/(1-t_ML_in_1[0][0]*t_ML_back_1[0][0])
                r_ML_in_2[1][1]=r_ML_in_1[1][1]*(1-ipow(t_ML_in_1[1][1]*t_ML_back_1[1][1], MLREP[i]-1))/(1-t_ML_in_1[1][1]*t_ML_back_1[1][1])

                product=cabsvalue(( t_ML_in_1[0][0]*t_ML_back_1[0][0] ) )
                normalizator=0.5*(1.+product)
                product=normalizator/cabsvalue(( t_ML_in_1[0][0] ) )

                t_ML_in_1[0][0]*=product
                t_ML_back_1[0][0]/=product

                product=cabsvalue(( t_ML_in_1[1][1]*t_ML_back_1[1][1] ) )
                normalizator=0.5*(1.+product)
                product=normalizator/cabsvalue(( t_ML_in_1[1][1] ) )

                t_ML_in_1[1][1]*=product
                t_ML_back_1[1][1]/=product


                t_ML_in_1[0][0]=ipow(t_ML_in_1[0][0],  MLREP[i]-1)
                t_ML_in_1[1][1]=ipow(t_ML_in_1[1][1],  MLREP[i]-1)
                t_ML_back_1[0][0]=ipow(t_ML_back_1[0][0],  MLREP[i]-1)
                t_ML_back_1[1][1]=ipow(t_ML_back_1[1][1],  MLREP[i]-1)

            else:
                ev1=0.5*( -sqrt( t_ML_in_1[1][1]*t_ML_in_1[1][1]-2*t_ML_in_1[0][0]*t_ML_in_1[1][1]+4*t_ML_in_1[1][0]*t_ML_in_1[0][1]+t_ML_in_1[0][0]*t_ML_in_1[0][0] )+t_ML_in_1[0][0]+t_ML_in_1[1][1]  )
                ev2=0.5*( sqrt( t_ML_in_1[1][1]*t_ML_in_1[1][1]-2*t_ML_in_1[0][0]*t_ML_in_1[1][1]+4*t_ML_in_1[1][0]*t_ML_in_1[0][1]+t_ML_in_1[0][0]*t_ML_in_1[0][0] )+t_ML_in_1[0][0]+t_ML_in_1[1][1]  )


                absevmax1=cabsvalue(( ev1 ) )
                absevmax2=cabsvalue(( ev2 ) )
                if(absevmax2>absevmax1):
                    absevmax1=absevmax2

                ev3=0.5*( -sqrt( t_ML_back_1[1][1]*t_ML_back_1[1][1]-2*t_ML_back_1[0][0]*t_ML_back_1[1][1]+4*t_ML_back_1[1][0]*t_ML_back_1[0][1]+t_ML_back_1[0][0]*t_ML_back_1[0][0] )+t_ML_back_1[0][0]+t_ML_back_1[1][1]  )
                ev4=0.5*( sqrt( t_ML_back_1[1][1]*t_ML_back_1[1][1]-2*t_ML_back_1[0][0]*t_ML_back_1[1][1]+4*t_ML_back_1[1][0]*t_ML_back_1[0][1]+t_ML_back_1[0][0]*t_ML_back_1[0][0] )+t_ML_back_1[0][0]+t_ML_back_1[1][1]  )


              #  print(ev3, abs(ev3), ev4, abs(ev4) )

                product= cabsvalue(( ev1*ev3 ) )
              #  print(product)
                product2= cabsvalue(( ev1*ev4 ) )
              #  print(product2)
                if( product2>product  ):
                    product=product2
                product2= cabsvalue(( ev2*ev3 ) )
               # print(product2)
                if( product2>product  ):
                    product=product2
                product2= cabsvalue(( ev2*ev4 ) )
              #  print(product2)
                if( product2>product  ):
                    product=product2
             #   print(product)
                normalizator=0.5*(1.+product)

                t_ML_in_1[0][0]=(normalizator/absevmax1)*( t_ML_in_1[0][0] )
                t_ML_in_1[0][1]=(normalizator/absevmax1)*( t_ML_in_1[0][1] )
                t_ML_in_1[1][0]=(normalizator/absevmax1)*( t_ML_in_1[1][0] )
                t_ML_in_1[1][1]=(normalizator/absevmax1)*( t_ML_in_1[1][1] )
                t_ML_back_1[0][0]=(absevmax1/normalizator)*( t_ML_back_1[0][0] )
                t_ML_back_1[0][1]=(absevmax1/normalizator)*( t_ML_back_1[0][1] )
                t_ML_back_1[1][0]=(absevmax1/normalizator)*( t_ML_back_1[1][0] )
                t_ML_back_1[1][1]=(absevmax1/normalizator)*( t_ML_back_1[1][1] )

                Calculate_Multilayer_equation(&t_ML_back_1, &t_ML_in_1, &r_ML_in_1, &r_ML_in_2, MLREP[i]-1)


            ML_is_diagonal=1


            Mult2x2_rightside(&t_ML_back_1, rtot)
            Mult2x2_leftside(rtot, &t_ML_in_1) # (t'(CA) p_C t'(BC) p_B t'(AB)*p(A))^N rtot (p(A) * t(AB) p_B t(BC) p_C t(CA))^N



            (rtot[0])[0][0]+=r_ML_in_2[0][0]
            (rtot[0])[1][0]+=r_ML_in_2[1][0]
            (rtot[0])[0][1]+=r_ML_in_2[0][1]
            (rtot[0])[1][1]+=r_ML_in_2[1][1]


            j=1
         #   print "9"
            while j<MLLENGTH[i]:
                LowerLayer=LR[MLCOMP[i][j-1]]
                UpperLayer=LR[MLCOMP[i][j]]
                if(Memory1.IsFilled):
                    Mempointer1=&Memory1
                    Mempointer2=&Memory2
                else:
                    Mempointer2=&Memory1
                    Mempointer1=&Memory2
                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0
                Mult2x2_rightside(&p, rtot)
                Mult2x2_leftside(rtot, &p)
                Mult2x2_rightside(&tprime, rtot)
                Mult2x2_leftside(rtot, &t)
                (rtot[0])[0][0]+=r[0][0]
                (rtot[0])[1][0]+=r[1][0]
                (rtot[0])[0][1]+=r[0][1]
                (rtot[0])[1][1]+=r[1][1]

                p[0][0]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz1)
                p[1][1]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz2)
                j=j+1
            LowerLayer=LR[MLCOMP[i][j-1]]
            if(i==Cap):
                if(Memory1.IsFilled):
                    Mempointer1=&Memory1
                    Mempointer2=&Memory2
                else:
                    Mempointer2=&Memory1
                    Mempointer1=&Memory2
                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0,0,0,0, \
                    LowerLayer.magdir, 0, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)

            else:
                UpperLayer=LR[MLCOMP[i+1][0]]
                if(Memory1.IsFilled):
                    Mempointer1=&Memory1
                    Mempointer2=&Memory2
                else:
                    Mempointer2=&Memory1
                    Mempointer1=&Memory2
                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0

             #   print "16"
         #   print "17"
            Mult2x2_rightside(&p, rtot)
            Mult2x2_leftside(rtot, &p)
            Mult2x2_rightside(&tprime, rtot)
            Mult2x2_leftside(rtot, &t)
            (rtot[0])[0][0]+=r[0][0]
            (rtot[0])[1][0]+=r[1][0]
            (rtot[0])[0][1]+=r[0][1]
            (rtot[0])[1][1]+=r[1][1]

            p[0][0]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz1)
            p[1][1]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz2)

        i=i+1



cdef void Paratt_magnetic_y_MS(Heterostructure* HS, double th, double wavelength, double complex (*rtot)[2][2]):

    cdef double k0=6.283185307179586/wavelength
    cdef double sintheta=sin(two_pi_div_360()*th)

    cdef double vy=cos(two_pi_div_360()*th)
    cdef double vyvy=quadr(vy)
    cdef double vy4=quadr(vyvy)
    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
    cdef CLayer UpperLayer, LowerLayer
    cdef rMemory Memory1, Memory2
    cdef rMemory* Mempointer1
    cdef rMemory* Mempointer2

    cdef int ML_is_diagonal=1

    cdef double complex r[2][2]
    cdef double complex rprime[2][2]
    cdef double complex t[2][2]
    cdef double complex tprime[2][2]
    cdef double complex p[2][2]

    cdef double complex t_ML_in_1[2][2] #Transmission through one Compound layer, down
    cdef double complex t_ML_back_1[2][2] #Transmission through one Compound layer, up
    cdef double complex t_ML_in_2[2][2] #Transmission through all Compound layer, down
    cdef double complex t_ML_back_2[2][2] #Transmission through all Compound layer, up
    cdef double complex r_ML_in_1[2][2] # Upper reflectivity from one Compound layer
    cdef double complex r_ML_in_2[2][2] # Upper reflectivity from the Multilayer
    cdef double complex r_ML_back_1[2][2] # Upper reflectivity from one Compound layer
    cdef double complex r_ML_back_2[2][2] # Upper reflectivity from the Multilayer
    cdef double complex C0[2][2]
    cdef double complex C1[2][2]
    cdef double complex ev1, ev2, ev3, ev4
    cdef double product, product2, absevmax1, absevmax2
    cdef double normalizator

    p[0][1]=0
    p[1][0]=0
   # print "0"
   # cdef double complex test=LR[MLCOMP[0][0]].ey

    cdef int Cap=NLAYERS-1
    LowerLayer=LR[MLCOMP[0][0]]

    if(LowerLayer.magdir):

        Fill_rMemory(&Memory1, vy,vyvy,vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)
    else:
        Memory1.epsy=1.0+LowerLayer.cy
        Memory1.epsz=1.0+LowerLayer.cz
        Memory1.vz1=sqrt(1.+LowerLayer.cx-vyvy)
        Memory1.vz2=sqrt((1.-vyvy/Memory1.epsz)*Memory1.epsy)

    if(NLAYERS==1):
        Calculate_rt_y(&Memory1, &Memory2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0, 0, 0, 0, \
                    LowerLayer.magdir, 0, rtot, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
    else:

        UpperLayer=LR[MLCOMP[1][0]]
        Calculate_rt_y(&Memory1, &Memory2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, rtot, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
        Memory1.IsFilled=0
        p[0][0]=exp(1j*k0*UpperLayer.Thickness*Memory2.vz1)
        p[1][1]=exp(1j*k0*UpperLayer.Thickness*Memory2.vz2)
    i=1
    while i<NLAYERS:
       # print "loop start", i
        if(MLLENGTH[i]==1):

            if(Memory1.IsFilled):
                Mempointer1=&Memory1
                Mempointer2=&Memory2
            else:
                Mempointer2=&Memory1
                Mempointer1=&Memory2
            LowerLayer=LR[MLCOMP[i][0]]
            if(i!=Cap):
                #Upper=MLCOMP[i+1][0]
                UpperLayer=LR[MLCOMP[i+1][0]]
               # print "i, Upper", i, Upper
                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0
            else:
                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0,0,0,0, \
                    LowerLayer.magdir, 0, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)

            FillC0(&C0, &rprime, rtot, &p)

            Mult2x2_rightside(&p, rtot)

            Mult2x2_leftside(rtot, &p)

            Mult2x2_rightside(&tprime, rtot)

            Mult2x2_leftside(rtot, &C0)


            Mult2x2_leftside(rtot, &t)

            (rtot[0])[0][0]+=r[0][0]
            (rtot[0])[1][0]+=r[1][0]
            (rtot[0])[0][1]+=r[0][1]
            (rtot[0])[1][1]+=r[1][1]


            if(i!=Cap):
                p[0][0]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz1)
                p[1][1]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz2)
        else: #Multilayer
            LowerLayer=LR[MLCOMP[i][0]]
            UpperLayer=LR[MLCOMP[i][1]]
            if(Memory1.IsFilled):
                Mempointer1=&Memory1
                Mempointer2=&Memory2
            else:
                Mempointer2=&Memory1
                Mempointer1=&Memory2

            Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, LowerLayer.Roughness, k0)



            (Mempointer1[0]).IsFilled=0
            if(LowerLayer.magdir):
                ML_is_diagonal=0


            Mult2x2_leftside(&t_ML_back_1, &p) # t'(AB)*p(A)


            p[0][0]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz1)
            p[1][1]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz2)

            j=1
            while j<MLLENGTH[i]:
               # Upper=MLCOMP[i][(j+1)%MLLENGTH[i]]



                UpperLayer=LR[ MLCOMP[i][(j+1)%MLLENGTH[i]] ]
                LowerLayer=LR[ MLCOMP[i][j] ]
                if(LowerLayer.magdir):
                    ML_is_diagonal=0
                if(Memory1.IsFilled):
                    Mempointer1=&Memory1
                    Mempointer2=&Memory2
                else:
                    Mempointer2=&Memory1
                    Mempointer1=&Memory2

                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0
             #   print "2", vz1, vz2, vz3, vz4
                FillC0(&C0, &rprime, &r_ML_in_1, &p)
                FillC0(&C1, &r_ML_in_1, &rprime, &p)
                Mult2x2_rightside(&p, &r_ML_in_1) #p_B r p_B and so on

                Mult2x2_leftside(&r_ML_in_1, &p)

                Mult2x2_rightside(&tprime, &r_ML_in_1) #t' p_B r p_b t and so on
                Mult2x2_leftside(&r_ML_in_1, &C0)
                Mult2x2_leftside(&r_ML_in_1, &t)

                r_ML_in_1[0][0]+=r[0][0]
                r_ML_in_1[1][0]+=r[1][0]
                r_ML_in_1[0][1]+=r[0][1]
                r_ML_in_1[1][1]+=r[1][1]



                Mult2x2_rightside(&C1, &t_ML_back_1)
                Mult2x2_rightside(&p, &t_ML_back_1)
               # Mult2x2_leftside(&t_ML_in_1, &t) # p(A) * t(AB) p_B t(BC) and so on #this comes later now
                Mult2x2_rightside(&tprime, &t_ML_back_1)# t'(BC) p_B t'(AB)*p(A) and so on
                p[0][0]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz1)
                p[1][1]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz2)
                j=j+1

            p[0][0]=exp(1j*k0*LowerLayer.Thickness*(Mempointer1[0]).vz1)
            p[1][1]=exp(1j*k0*LowerLayer.Thickness*(Mempointer1[0]).vz2)
         #   print "p C components"
         #   print MLLENGTH[i]-1, LR[MLCOMP[i][MLLENGTH[i]-1]].Thickness, vz1, vz2
            r_ML_back_1[0][0]=rprime[0][0]
            r_ML_back_1[1][0]=rprime[1][0]
            r_ML_back_1[0][1]=rprime[0][1]
            r_ML_back_1[1][1]=rprime[1][1]
            t_ML_in_1[0][0]=t[0][0]
            t_ML_in_1[0][1]=t[0][1]
            t_ML_in_1[1][0]=t[1][0]
            t_ML_in_1[1][1]=t[1][1]
            Mult2x2_rightside(&p, &t_ML_in_1)
         #   Mult2x2_leftside(&t_ML_in_1, &p)
            j=MLLENGTH[i]-1


            while j>0:
                UpperLayer=LR[ MLCOMP[i][j] ]
                LowerLayer=LR[ MLCOMP[i][j-1] ]
                Mempointer1=&Memory1
                Mempointer2=&Memory2

                if(LowerLayer.magdir):

                    Fill_rMemory(Mempointer1, vy,vyvy,vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)
                else:
                    Memory1.epsy=1.0+LowerLayer.cy
                    Memory1.epsz=1.0+LowerLayer.cz
                    Memory1.vz1=sqrt(1.+LowerLayer.cx-vyvy)
                    Memory1.vz2=sqrt((1.-vyvy/Memory1.epsz)*Memory1.epsy)


                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0


                FillC0(&C0, &r, &r_ML_back_1, &p)
               # FillC0(&C1, &r_ML_back_1, &r, &p)
                C1[0][0]=p[0][0]
                C1[1][0]=0
                C1[0][1]=0
                C1[1][1]=p[1][1]
                Mult2x2_leftside(&C1, &r_ML_back_1)
                Mult2x2_leftside(&C1, &p)
                Mult2x2_leftside(&C1, &r)
                C1[0][0]=1.-C1[0][0]
                C1[1][0]=-C1[1][0]
                C1[0][1]=-C1[0][1]
                C1[1][1]=1.-C1[1][1]
                Invert2x2(&C1)

                Mult2x2_rightside(&p, &r_ML_back_1) #p_B r p_B and so on
                Mult2x2_leftside(&r_ML_back_1, &p)
                Mult2x2_rightside(&t, &r_ML_back_1)
                Mult2x2_leftside(&r_ML_back_1, &C0)
                Mult2x2_leftside(&r_ML_back_1, &tprime)

                r_ML_back_1[0][0]+=rprime[0][0]
                r_ML_back_1[1][0]+=rprime[1][0]
                r_ML_back_1[0][1]+=rprime[0][1]
                r_ML_back_1[1][1]+=rprime[1][1]


                Mult2x2_rightside(&C1, &t_ML_in_1)
             #   Mult2x2_rightside(&p, &t_ML_in_1)
                Mult2x2_rightside(&t, &t_ML_in_1)
                p[0][0]=exp(1j*k0*LowerLayer.Thickness*(Mempointer1[0]).vz1)
                p[1][1]=exp(1j*k0*LowerLayer.Thickness*(Mempointer1[0]).vz2)

                Mult2x2_rightside(&p, &t_ML_in_1)

                j=j-1


            Mult2x2_rightside(&p, &r_ML_back_1)
            Mult2x2_leftside(&r_ML_back_1, &p)



            ev1=0.5*( -sqrt( t_ML_in_1[1][1]*t_ML_in_1[1][1]-2*t_ML_in_1[0][0]*t_ML_in_1[1][1]+4*t_ML_in_1[1][0]*t_ML_in_1[0][1]+t_ML_in_1[0][0]*t_ML_in_1[0][0] )+t_ML_in_1[0][0]+t_ML_in_1[1][1]  )
            ev2=0.5*( sqrt( t_ML_in_1[1][1]*t_ML_in_1[1][1]-2*t_ML_in_1[0][0]*t_ML_in_1[1][1]+4*t_ML_in_1[1][0]*t_ML_in_1[0][1]+t_ML_in_1[0][0]*t_ML_in_1[0][0] )+t_ML_in_1[0][0]+t_ML_in_1[1][1]  )


            absevmax1=cabsvalue(( ev1 ) )
            absevmax2=cabsvalue(( ev2 ) )
            if(absevmax2>absevmax1):
                absevmax1=absevmax2

            ev3=0.5*( -sqrt( t_ML_back_1[1][1]*t_ML_back_1[1][1]-2*t_ML_back_1[0][0]*t_ML_back_1[1][1]+4*t_ML_back_1[1][0]*t_ML_back_1[0][1]+t_ML_back_1[0][0]*t_ML_back_1[0][0] )+t_ML_back_1[0][0]+t_ML_back_1[1][1]  )
            ev4=0.5*( sqrt( t_ML_back_1[1][1]*t_ML_back_1[1][1]-2*t_ML_back_1[0][0]*t_ML_back_1[1][1]+4*t_ML_back_1[1][0]*t_ML_back_1[0][1]+t_ML_back_1[0][0]*t_ML_back_1[0][0] )+t_ML_back_1[0][0]+t_ML_back_1[1][1]  )


          #  print(ev3, abs(ev3), ev4, abs(ev4) )

            product= cabsvalue(( ev1*ev3 ) )
          #  print(product)
            product2= cabsvalue(( ev1*ev4 ) )
          #  print(product2)
            if( product2>product  ):
                product=product2
            product2= cabsvalue(( ev2*ev3 ) )
           # print(product2)
            if( product2>product  ):
                product=product2
            product2= cabsvalue(( ev2*ev4 ) )
          #  print(product2)
            if( product2>product  ):
                product=product2
         #   print(product)
            normalizator=0.5*(1.+product)

            t_ML_in_1[0][0]=(normalizator/absevmax1)*( t_ML_in_1[0][0] )
            t_ML_in_1[0][1]=(normalizator/absevmax1)*( t_ML_in_1[0][1] )
            t_ML_in_1[1][0]=(normalizator/absevmax1)*( t_ML_in_1[1][0] )
            t_ML_in_1[1][1]=(normalizator/absevmax1)*( t_ML_in_1[1][1] )
            t_ML_back_1[0][0]=(absevmax1/normalizator)*( t_ML_back_1[0][0] )
            t_ML_back_1[0][1]=(absevmax1/normalizator)*( t_ML_back_1[0][1] )
            t_ML_back_1[1][0]=(absevmax1/normalizator)*( t_ML_back_1[1][0] )
            t_ML_back_1[1][1]=(absevmax1/normalizator)*( t_ML_back_1[1][1] )


            Calculate_Multilayer_with_Matrices(&t_ML_back_1, &t_ML_back_2,&t_ML_in_1, &t_ML_in_2, &r_ML_in_1, &r_ML_in_2, &r_ML_back_1, &r_ML_back_2, MLREP[i]-1)
          #  print t_ML_back_2[1][1], t_ML_in_2[1][1], r_ML_in_2[1][1], r_ML_back_2[1][1]



            ML_is_diagonal=1


            C0[0][0]=(rtot[0])[0][0]
            C0[0][1]=(rtot[0])[0][1]
            C0[1][0]=(rtot[0])[1][0]
            C0[1][1]=(rtot[0])[1][1]
            Mult2x2_rightside(&r_ML_back_2, &C0)
            C0[0][0]=1-C0[0][0]
            C0[0][1]=-C0[0][1]
            C0[1][0]=-C0[1][0]
            C0[1][1]=1-C0[1][1]
            Invert2x2(&C0)


            Mult2x2_rightside(&t_ML_back_2, rtot)
            Mult2x2_leftside(rtot, &C0)
            Mult2x2_leftside(rtot, &t_ML_in_2) # (t'(CA) p_C t'(BC) p_B t'(AB)*p(A))^N rtot (p(A) * t(AB) p_B t(BC) p_C t(CA))^N
            (rtot[0])[0][0]+=r_ML_in_2[0][0]
            (rtot[0])[1][0]+=r_ML_in_2[1][0]
            (rtot[0])[0][1]+=r_ML_in_2[0][1]
            (rtot[0])[1][1]+=r_ML_in_2[1][1]


            LowerLayer=LR[MLCOMP[i][0]]
            if(LowerLayer.magdir):
                Fill_rMemory(&Memory1, vy,vyvy,vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)
            else:
                Memory1.epsy=1.0+LowerLayer.cy
                Memory1.epsz=1.0+LowerLayer.cz
                Memory1.vz1=sqrt(1.+LowerLayer.cx-vyvy)
                Memory1.vz2=sqrt((1.-vyvy/Memory1.epsz)*Memory1.epsy)
            Mempointer1=&Memory1
            Mempointer2=&Memory2
            Memory1.IsFilled=1
            j=1

            while j<MLLENGTH[i]:
                LowerLayer=LR[MLCOMP[i][j-1]]
                UpperLayer=LR[MLCOMP[i][j]]
                if(Memory1.IsFilled):
                    Mempointer1=&Memory1
                    Mempointer2=&Memory2
                else:
                    Mempointer2=&Memory1
                    Mempointer1=&Memory2
                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0
                FillC0(&C0, &rprime, rtot, &p)
                Mult2x2_rightside(&p, rtot)
                Mult2x2_leftside(rtot, &p)
                Mult2x2_rightside(&tprime, rtot)
                Mult2x2_leftside(rtot, &C0)
                Mult2x2_leftside(rtot, &t)
                (rtot[0])[0][0]+=r[0][0]
                (rtot[0])[1][0]+=r[1][0]
                (rtot[0])[0][1]+=r[0][1]
                (rtot[0])[1][1]+=r[1][1]

                p[0][0]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz1)
                p[1][1]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz2)
                j=j+1
              #  print "12"
            LowerLayer=LR[MLCOMP[i][j-1]]
            if(i==Cap):
                if(Memory1.IsFilled):
                    Mempointer1=&Memory1
                    Mempointer2=&Memory2
                else:
                    Mempointer2=&Memory1
                    Mempointer1=&Memory2
                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0,0,0,0, \
                    LowerLayer.magdir, 0, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)

            else:
                UpperLayer=Layer=LR[MLCOMP[i+1][0]]
                if(Memory1.IsFilled):
                    Mempointer1=&Memory1
                    Mempointer2=&Memory2
                else:
                    Mempointer2=&Memory1
                    Mempointer1=&Memory2
                Calculate_rt_y(Mempointer1, Mempointer2, vy, vyvy, vy4, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0
             #   print "5", LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness
          #  print "5", vz1, vz2, vz3, vz4
         #   print "17"
            FillC0(&C0, &rprime, rtot, &p)
            Mult2x2_rightside(&p, rtot)
            Mult2x2_leftside(rtot, &p)
            Mult2x2_rightside(&tprime, rtot)
            Mult2x2_leftside(rtot, &C0)
            Mult2x2_leftside(rtot, &t)
            (rtot[0])[0][0]+=r[0][0]
            (rtot[0])[1][0]+=r[1][0]
            (rtot[0])[0][1]+=r[0][1]
            (rtot[0])[1][1]+=r[1][1]


            p[0][0]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz1)
            p[1][1]=exp(1j*k0*UpperLayer.Thickness*(Mempointer2[0]).vz2)

        i=i+1
