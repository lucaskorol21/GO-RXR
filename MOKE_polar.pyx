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


cdef void Fill_rMemory_z(rMemory *Mem, double vy, double vyvy, double omvyvy, double complex chix, double complex chiy, double complex chiz, double complex chig):
    cdef double complex xld, magrest, xldsquared, diftester, diftester2

    (Mem[0]).epsy=chiy+1.0
    (Mem[0]).epsz=chiz+1.0
    xld=(chix-chiy) + vyvy*(chiy-chiz)/(Mem[0]).epsz
    xldsquared = cquadr(xld)

    (Mem[0]).C1=(1-vyvy/(Mem[0]).epsz)
    (Mem[0]).C2=(Mem[0]).epsy/(Mem[0]).epsz
    (Mem[0]).C3=(Mem[0]).C1 #(1-vyvy/(Mem[0]).epsz)
    (Mem[0]).B1=(((Mem[0]).epsy/(Mem[0]).epsz+1)*omvyvy)/2
    (Mem[0]).B2=(chiz*(Mem[0]).C2+chix)/2
    (Mem[0]).root=sqrt(0.25*xldsquared-chig*chig*(Mem[0]).C1)



    (Mem[0]).B=(Mem[0]).B1+(Mem[0]).B2
    (Mem[0]).vz1=sqrt(((Mem[0]).B+(Mem[0]).root))
    (Mem[0]).vz2=sqrt(((Mem[0]).B-(Mem[0]).root))

    diftester= xld/2 -(Mem[0]).root
    diftester2=xld/2 + (Mem[0]).root

    if( Cmaxnorm(diftester2) < Cmaxnorm(diftester)  ):
        (Mem[0]).PHI1=( diftester )/chig
    else:
        (Mem[0]).PHI1= chig*(Mem[0]).C1/( diftester2 )


    (Mem[0]).PHI2=(Mem[0]).PHI1/(Mem[0]).vz2
    #(Mem[0]).PHI1=-((Mem[0]).B2+(Mem[0]).root+   ( ((chiy+chiz)*omvyvy)/2 -chix-(1+chix-vyvy)*chiz )/(Mem[0]).epsz )/chig


    #(Mem[0]).PHI2=-((Mem[0]).root-(Mem[0]).B2+( chiz*(Mem[0]).epsy+ 0.5*(chiy-chiz)*omvyvy  )/(Mem[0]).epsz   )/(chig*(Mem[0]).vz2)



    (Mem[0]).PHI3=(Mem[0]).C3/(Mem[0]).vz2
    (Mem[0]).IsFilled=1


cdef void Calculate_rt_z(rMemory *Mem1, rMemory *Mem2, double vy, double vyvy, double omvyvy, double complex chix1, double complex chiy1, double complex chiz1, double complex chig1, double complex chix2, double complex chiy2, double complex chiz2, double complex chig2, \
                    int IsMagnetic1, int IsMagnetic2, double complex (*r)[2][2], double complex (*rprime)[2][2], double complex (*t)[2][2], double complex (*tprime)[2][2], double sigma, double k0):


    cdef double complex J[2][4]
    cdef double complex difB
    cdef double complex difvz1vz3
    cdef double complex difvz1vz4
    cdef double complex difvz2vz3
    cdef double complex difvz2vz4


    cdef double complex divide2
    cdef double complex fac1
    cdef double roughfac=-0.5*quadr(sigma*  k0)
    cdef double complex rough
    cdef double complex save1
    cdef double complex save3
    cdef double complex save4
    cdef double complex ez1ez2
    if( IsMagnetic1 & IsMagnetic2 ):
        Fill_rMemory_z(Mem2, vy, vyvy, omvyvy, chix2, chiy2, chiz2, chig2)

        difB=omvyvy*(chiy1+chiz2*(1+chiy1)-chiy2-chiz1*(1+chiy2))/(2*(Mem1[0]).epsz*(Mem2[0]).epsz)+(Mem1[0]).B2-(Mem2[0]).B2
        difvz1vz3=(difB+(Mem1[0]).root-(Mem2[0]).root)/((Mem1[0]).vz1+(Mem2[0]).vz1)
        difvz1vz4=(difB+(Mem1[0]).root+(Mem2[0]).root)/((Mem1[0]).vz1+(Mem2[0]).vz2)
        difvz2vz3=(difB-(Mem1[0]).root-(Mem2[0]).root)/((Mem1[0]).vz2+(Mem2[0]).vz1)
        difvz2vz4=(difB-(Mem1[0]).root+(Mem2[0]).root)/((Mem1[0]).vz2+(Mem2[0]).vz2)

    #    divide1=1.0/(2*(  (Mem1[0]).PHI3    - (Mem1[0]).PHI2 *(Mem1[0]).PHI1))
        ez1ez2=(Mem1[0]).epsz*(Mem2[0]).epsz
        J[0][0]=( (Mem1[0]).PHI3*(1+(Mem2[0]).vz1/(Mem1[0]).vz1) - (Mem2[0]).PHI1*(Mem1[0]).PHI2*(1+((Mem2[0]).vz1/(Mem1[0]).vz1)*( (Mem1[0]).C1/(Mem2[0]).C1 ) ) )
        J[0][2]=( (Mem1[0]).PHI3*difvz1vz3/(Mem1[0]).vz1 - (Mem2[0]).PHI1*(Mem1[0]).PHI2/((Mem2[0]).C1*(Mem1[0]).vz1)*(difvz1vz3-(difvz1vz3+chiz1*(Mem1[0]).vz1-(Mem2[0]).vz1*chiz2)*vyvy/(ez1ez2)  ) )

        J[1][1]=( -(Mem1[0]).PHI1*(Mem2[0]).PHI2*((Mem2[0]).vz2/(Mem1[0]).vz2+1)+((Mem1[0]).C1/(Mem1[0]).vz2+(Mem2[0]).C1/(Mem2[0]).vz2) )
        J[1][3]=( (Mem1[0]).PHI1*(Mem2[0]).PHI2*difvz2vz4/(Mem1[0]).vz2+(-difvz2vz4-( -difvz2vz4+(Mem2[0]).vz2*chiz2-(Mem1[0]).vz2*chiz1 )*vyvy/(ez1ez2) )/((Mem1[0]).vz2*(Mem2[0]).vz2) )

        J[0][1]=( ((Mem2[0]).vz2/(Mem1[0]).vz1+1)*(Mem2[0]).PHI2*(Mem1[0]).C1/(Mem1[0]).vz2 - (Mem1[0]).PHI2*((Mem1[0]).C1/(Mem1[0]).vz1+(Mem2[0]).C1/(Mem2[0]).vz2))
        J[0][3]=( -difvz1vz4*(Mem2[0]).PHI2*(Mem1[0]).C1/((Mem1[0]).vz2*(Mem1[0]).vz1) - (-difvz1vz4-vyvy*( -difvz1vz4+chiz2*(Mem2[0]).vz2-chiz1*(Mem1[0]).vz1 )/(ez1ez2)  )*(Mem1[0]).PHI2/((Mem1[0]).vz1*(Mem2[0]).vz2)   )

        J[1][0]=( (Mem2[0]).PHI1*(1+(Mem1[0]).C1*(Mem2[0]).vz1/((Mem2[0]).C1*(Mem1[0]).vz2) )-(Mem1[0]).PHI1*(1+(Mem2[0]).vz1/(Mem1[0]).vz2) )
        J[1][2]=( (Mem2[0]).PHI1*(difvz2vz3- ( difvz2vz3+chiz1*(Mem1[0]).vz2-chiz2*(Mem2[0]).vz1 )*vyvy/(ez1ez2) )/((Mem2[0]).C1*(Mem1[0]).vz2)-(Mem1[0]).PHI1*difvz2vz3/(Mem1[0]).vz2 )
#        J[0][0]=( (Mem1[0]).PHI3*(1+(Mem2[0]).vz1/(Mem1[0]).vz1) - (Mem2[0]).PHI1*(Mem1[0]).PHI2*(1+((Mem2[0]).vz1/(Mem1[0]).vz1)*( (Mem1[0]).C1/(Mem2[0]).C1 ) ) )*divide1
#        J[0][2]=( (Mem1[0]).PHI3*difvz1vz3/(Mem1[0]).vz1 - (Mem2[0]).PHI1*(Mem1[0]).PHI2/((Mem2[0]).C1*(Mem1[0]).vz1)*(difvz1vz3-(difvz1vz3+chiz1*(Mem1[0]).vz1-(Mem2[0]).vz1*chiz2)*vyvy/(ez1ez2)  ) )*divide1
#
#        J[1][1]=( -(Mem1[0]).PHI1*(Mem2[0]).PHI2*((Mem2[0]).vz2/(Mem1[0]).vz2+1)+((Mem1[0]).C1/(Mem1[0]).vz2+(Mem2[0]).C1/(Mem2[0]).vz2) )*divide1
#        J[1][3]=( (Mem1[0]).PHI1*(Mem2[0]).PHI2*difvz2vz4/(Mem1[0]).vz2+(-difvz2vz4-( -difvz2vz4+(Mem2[0]).vz2*chiz2-(Mem1[0]).vz2*chiz1 )*vyvy/(ez1ez2) )/((Mem1[0]).vz2*(Mem2[0]).vz2) )*divide1
#
#        J[0][1]=( ((Mem2[0]).vz2/(Mem1[0]).vz1+1)*(Mem2[0]).PHI2*(Mem1[0]).C1/(Mem1[0]).vz2 - (Mem1[0]).PHI2*((Mem1[0]).C1/(Mem1[0]).vz1+(Mem2[0]).C1/(Mem2[0]).vz2))*divide1
#        J[0][3]=( -difvz1vz4*(Mem2[0]).PHI2*(Mem1[0]).C1/((Mem1[0]).vz2*(Mem1[0]).vz1) - (-difvz1vz4-vyvy*( -difvz1vz4+chiz2*(Mem2[0]).vz2-chiz1*(Mem1[0]).vz1 )/(ez1ez2)  )*(Mem1[0]).PHI2/((Mem1[0]).vz1*(Mem2[0]).vz2)   )*divide1
#
#        J[1][0]=( (Mem2[0]).PHI1*(1+(Mem1[0]).C1*(Mem2[0]).vz1/((Mem2[0]).C1*(Mem1[0]).vz2) )-(Mem1[0]).PHI1*(1+(Mem2[0]).vz1/(Mem1[0]).vz2) )*divide1
#        J[1][2]=( (Mem2[0]).PHI1*(difvz2vz3- ( difvz2vz3+chiz1*(Mem1[0]).vz2-chiz2*(Mem2[0]).vz1 )*vyvy/(ez1ez2) )/((Mem2[0]).C1*(Mem1[0]).vz2)-(Mem1[0]).PHI1*difvz2vz3/(Mem1[0]).vz2 )*divide1


    elif(IsMagnetic1):

        (Mem2[0]).epsy=chiy2+1.0
        (Mem2[0]).epsz=chiz2+1.0
        (Mem2[0]).vz1=sqrt(omvyvy+chix2)
        (Mem2[0]).vz2=sqrt( (1-vyvy/(Mem2[0]).epsz)*(Mem2[0]).epsy )
        (Mem2[0]).IsFilled=1
        ez1ez2=(Mem1[0]).epsz*(Mem2[0]).epsz
        difB=(((chiz1+chiy1)*omvyvy)*0.5-chiz1*omvyvy)/(Mem1[0]).epsz-chix2
        difvz1vz3=( (Mem1[0]).B2+(Mem1[0]).root+difB )/((Mem1[0]).vz1+(Mem2[0]).vz1)
        difvz2vz3=( (Mem1[0]).B2-(Mem1[0]).root+difB )/((Mem1[0]).vz2+(Mem2[0]).vz1)
        difB=  ( chiz2*(omvyvy-(Mem2[0]).epsy)+0.5*(chiy1+chiz1)*omvyvy*(Mem2[0]).epsz- omvyvy*chiy2)/ez1ez2 -chiz1*(1-vyvy/(Mem2[0]).epsz)*(Mem2[0]).epsy/(Mem1[0]).epsz
        difvz1vz4=( (Mem1[0]).B2+(Mem1[0]).root+difB )/((Mem1[0]).vz1+(Mem2[0]).vz2)
        difvz2vz4=( (Mem1[0]).B2-(Mem1[0]).root+difB )/((Mem1[0]).vz2+(Mem2[0]).vz2)

        save3=(1-vyvy/(Mem1[0]).epsz)
        save4=(1-vyvy/(Mem2[0]).epsz)

#        divide1=1.0/(2*(  (Mem1[0]).PHI3    - (Mem1[0]).PHI2 *(Mem1[0]).PHI1))
#        save1=save3*divide1/((Mem1[0]).vz1*(Mem1[0]).vz2)
#        J[0][0]=save1*( (Mem1[0]).vz1+(Mem2[0]).vz1 )
#        J[0][2]=save1*difvz1vz3
#
#        save1=divide1/((Mem2[0]).vz2*(Mem1[0]).vz2)
#
#        J[1][1]=save1*( save3*(Mem2[0]).vz2+save4*(Mem1[0]).vz2  )
#        J[1][3]=save1*( -difvz2vz4-(vyvy/ez1ez2)*( -difvz2vz4+ (Mem2[0]).vz2*chiz2-(Mem1[0]).vz2*chiz1 ) )
#
#        save1=-(Mem1[0]).PHI2*divide1/( (Mem2[0]).vz2*(Mem1[0]).vz1 )
#        J[0][1]=save1*( save3*(Mem2[0]).vz2 + save4*(Mem1[0]).vz1 )
#        J[0][3]=save1*( -difvz1vz4-(vyvy/ez1ez2)*( -difvz1vz4+ (Mem2[0]).vz2*chiz2-(Mem1[0]).vz1*chiz1 ) )
#
#        save1=-(Mem1[0]).PHI1*divide1/(Mem1[0]).vz2
#        J[1][0]=save1*( (Mem1[0]).vz2+(Mem2[0]).vz1 )
#        J[1][2]=save1*difvz2vz3
        save1=save3/((Mem1[0]).vz1*(Mem1[0]).vz2)
        J[0][0]=save1*( (Mem1[0]).vz1+(Mem2[0]).vz1 )
        J[0][2]=save1*difvz1vz3

        save1=1.0/((Mem2[0]).vz2*(Mem1[0]).vz2)

        J[1][1]=save1*( save3*(Mem2[0]).vz2+save4*(Mem1[0]).vz2  )
        J[1][3]=save1*( -difvz2vz4-(vyvy/ez1ez2)*( -difvz2vz4+ (Mem2[0]).vz2*chiz2-(Mem1[0]).vz2*chiz1 ) )

        save1=-(Mem1[0]).PHI2/( (Mem2[0]).vz2*(Mem1[0]).vz1 )
        J[0][1]=save1*( save3*(Mem2[0]).vz2 + save4*(Mem1[0]).vz1 )
        J[0][3]=save1*( -difvz1vz4-(vyvy/ez1ez2)*( -difvz1vz4+ (Mem2[0]).vz2*chiz2-(Mem1[0]).vz1*chiz1 ) )

        save1=-(Mem1[0]).PHI1/(Mem1[0]).vz2
        J[1][0]=save1*( (Mem1[0]).vz2+(Mem2[0]).vz1 )
        J[1][2]=save1*difvz2vz3

    elif(IsMagnetic2):
        Fill_rMemory_z(Mem2, vy, vyvy, omvyvy, chix2, chiy2, chiz2, chig2)

        difB=chix1+( omvyvy*chiz2-(chiy2+chiz2)*omvyvy*0.5 )/(Mem2[0]).epsz
        ez1ez2=(Mem1[0]).epsz*(Mem2[0]).epsz
        difvz1vz3=( -(Mem2[0]).B2-(Mem2[0]).root+difB )/((Mem1[0]).vz1+(Mem2[0]).vz1)
        difvz1vz4=( -(Mem2[0]).B2+(Mem2[0]).root+difB )/((Mem1[0]).vz1+(Mem2[0]).vz2)
        difB=-( ( chiz1*(omvyvy-(Mem1[0]).epsy)+0.5*(chiy2+chiz2)*omvyvy*(Mem1[0]).epsz- omvyvy*chiy1)/ez1ez2 -chiz2*(1-vyvy/(Mem1[0]).epsz)*(Mem1[0]).epsy/(Mem2[0]).epsz)
        difvz2vz3=( -(Mem2[0]).B2-(Mem2[0]).root+difB )/((Mem1[0]).vz2+(Mem2[0]).vz1)
        difvz2vz4=( -(Mem2[0]).B2+(Mem2[0]).root+difB )/((Mem1[0]).vz2+(Mem2[0]).vz2)
        save3=(1-vyvy/(Mem1[0]).epsz)
        save4=(1-vyvy/(Mem2[0]).epsz)

#        save1=0.5/(Mem1[0]).vz1
#        J[0][0]=save1*((Mem1[0]).vz1+(Mem2[0]).vz1)
#        J[0][2]=save1*difvz1vz3
#
#
#        save1=0.5/( save3*(Mem2[0]).vz2 )
#        J[1][1]=save1*( save3*(Mem2[0]).vz2 + save4*(Mem1[0]).vz2 )
#        J[1][3]=save1*( -difvz2vz4-(vyvy/ez1ez2)*( -difvz2vz4+ (Mem2[0]).vz2*chiz2 - (Mem1[0]).vz2*chiz1 ) )
#
#
#
#        save1=0.5*(Mem2[0]).PHI2/(Mem1[0]).vz1
#
#        J[0][1]=save1*((Mem2[0]).vz2+(Mem1[0]).vz1)
#        J[0][3]=-save1*difvz1vz4
#
#
#        save1=(Mem2[0]).PHI1/(2*save3*save4)
#        J[1][0]=save1*( (Mem1[0]).vz2*save4+(Mem2[0]).vz1*save3 )
#        J[1][2]=save1*( difvz2vz3+  -vyvy/ez1ez2*( difvz2vz3+ (Mem1[0]).vz2*chiz1-(Mem2[0]).vz1*chiz2   ) )

        save1=1.0/(Mem1[0]).vz1
        J[0][0]=save1*((Mem1[0]).vz1+(Mem2[0]).vz1)
        J[0][2]=save1*difvz1vz3


        save1=1.0/( save3*(Mem2[0]).vz2 )
        J[1][1]=save1*( save3*(Mem2[0]).vz2 + save4*(Mem1[0]).vz2 )
        J[1][3]=save1*( -difvz2vz4-(vyvy/ez1ez2)*( -difvz2vz4+ (Mem2[0]).vz2*chiz2 - (Mem1[0]).vz2*chiz1 ) )



        save1=1.0*(Mem2[0]).PHI2/(Mem1[0]).vz1

        J[0][1]=save1*((Mem2[0]).vz2+(Mem1[0]).vz1)
        J[0][3]=-save1*difvz1vz4


        save1=(Mem2[0]).PHI1/(save3*save4)
        J[1][0]=save1*( (Mem1[0]).vz2*save4+(Mem2[0]).vz1*save3 )
        J[1][2]=save1*( difvz2vz3+  -vyvy/ez1ez2*( difvz2vz3+ (Mem1[0]).vz2*chiz1-(Mem2[0]).vz1*chiz2   ) )

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



    J[0][0]*=exp(roughfac*cquadr(difvz1vz3) )
    J[0][2]*=exp(roughfac*cquadr((Mem1[0]).vz1+(Mem2[0]).vz1))
    J[1][1]*=exp(roughfac*cquadr(difvz2vz4) )
    J[1][3]*=exp(roughfac*cquadr((Mem1[0]).vz2+(Mem2[0]).vz2))
    J[0][1]*=exp(roughfac*cquadr(difvz1vz4) )
    J[0][3]*=exp(roughfac*cquadr((Mem1[0]).vz1+(Mem2[0]).vz2))
    J[1][0]*=exp(roughfac*cquadr(difvz2vz3) )
    J[1][2]*=exp(roughfac*cquadr((Mem1[0]).vz2+(Mem2[0]).vz1))

    divide2=1.0/(J[0][1]*J[1][0]-J[1][1]*J[0][0])

    (r[0])[0][0]=(J[1][1]*J[0][2]-J[0][1]*J[1][2])*divide2 # Incoming 1 reflected 1
    (r[0])[0][1]=-(J[1][1]*J[0][3]-J[0][1]*J[1][3])*divide2 # Incoming 2 reflected 1
    (r[0])[1][0]=-(J[0][0]*J[1][2]-J[1][0]*J[0][2])*divide2 # Incoming 1 reflected 2
    (r[0])[1][1]=(J[0][0]*J[1][3]-J[1][0]*J[0][3])*divide2 # Incoming 2 reflected 2



    (t[0])[0][0]=J[0][2]*(r[0])[0][0]+J[0][3]*(r[0])[1][0]+J[0][0] # Incoming 1 transmitted 1
    (t[0])[0][1]=J[0][2]*(r[0])[0][1]+J[0][3]*(r[0])[1][1]+J[0][1] # Incoming 2 transmitted 1
    (t[0])[1][0]=J[1][2]*(r[0])[0][0]+J[1][3]*(r[0])[1][0]+J[1][0] # Incoming 1 transmitted 2
    (t[0])[1][1]=J[1][2]*(r[0])[0][1]+J[1][3]*(r[0])[1][1]+J[1][1] # Incoming 2 transmitted 2



    (tprime[0])[0][0]=-J[1][1]*divide2
    (tprime[0])[0][1]=-J[0][1]*divide2
    (tprime[0])[1][0]=-J[1][0]*divide2
    (tprime[0])[1][1]=-J[0][0]*divide2

    (rprime[0])[0][0]=J[0][2]*(tprime[0])[0][0]+J[0][3]*((tprime[0])[1][0])
    (rprime[0])[0][1]=J[0][2]*((tprime[0])[0][1])+J[0][3]*((tprime[0])[1][1])
    (rprime[0])[1][0]=J[1][2]*((tprime[0])[0][0])+J[1][3]*((tprime[0])[1][0])
    (rprime[0])[1][1]=J[1][2]*((tprime[0])[0][1])+J[1][3]*((tprime[0])[1][1])

cdef void Paratt_magnetic_z(Heterostructure* HS, double th, double wavelength, double complex (*rtot)[2][2]):

    cdef double k0=6.283185307179586/wavelength
    cdef double sintheta=sin(two_pi_div_360()*th)

    cdef double vy=cos(two_pi_div_360()*th)
    cdef double vyvy=quadr(vy)
    cdef double omvyvy=1-vyvy
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

        Fill_rMemory_z(&Memory1, vy,vyvy,omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)
    else:
        Memory1.epsy=1.0+LowerLayer.cy
        Memory1.epsz=1.0+LowerLayer.cz
        Memory1.vz1=sqrt(1.+LowerLayer.cx-vyvy)
        Memory1.vz2=sqrt((1.-vyvy/Memory1.epsz)*Memory1.epsy)

    if(NLAYERS==1):
        Calculate_rt_z(&Memory1, &Memory2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0, 0, 0, 0, \
                    LowerLayer.magdir, 0, rtot, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
    else:

        UpperLayer=LR[MLCOMP[1][0]]
        Calculate_rt_z(&Memory1, &Memory2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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
                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0
            else:
                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0,0,0,0, \
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

            Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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

                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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
                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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
                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0,0,0,0, \
                    LowerLayer.magdir, 0, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)

            else:
                UpperLayer=LR[MLCOMP[i+1][0]]
                if(Memory1.IsFilled):
                    Mempointer1=&Memory1
                    Mempointer2=&Memory2
                else:
                    Mempointer2=&Memory1
                    Mempointer1=&Memory2
                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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


cdef void Paratt_magnetic_z_MS(Heterostructure* HS, double th, double wavelength, double complex (*rtot)[2][2]):

    cdef double k0=6.283185307179586/wavelength
    cdef double sintheta=sin(two_pi_div_360()*th)

    cdef double vy=cos(two_pi_div_360()*th)
    cdef double vyvy=quadr(vy)
    cdef double omvyvy=1-vyvy
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

        Fill_rMemory_z(&Memory1, vy,vyvy,omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)
    else:
        Memory1.epsy=1.0+LowerLayer.cy
        Memory1.epsz=1.0+LowerLayer.cz
        Memory1.vz1=sqrt(1.+LowerLayer.cx-vyvy)
        Memory1.vz2=sqrt((1.-vyvy/Memory1.epsz)*Memory1.epsy)

    if(NLAYERS==1):
        Calculate_rt_z(&Memory1, &Memory2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0, 0, 0, 0, \
                    LowerLayer.magdir, 0, rtot, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
    else:

        UpperLayer=LR[MLCOMP[1][0]]
        Calculate_rt_z(&Memory1, &Memory2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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
                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
                    LowerLayer.magdir, UpperLayer.magdir, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)
                (Mempointer1[0]).IsFilled=0
            else:
                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0,0,0,0, \
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

            Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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

                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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
#            if(wavelength<23):
#                if(th>4):
#                    if( th<5):
#                        print("tml1\n",t_ML_in_1[1][1] )

            while j>0:
                UpperLayer=LR[ MLCOMP[i][j] ]
                LowerLayer=LR[ MLCOMP[i][j-1] ]

                Mempointer1=&Memory1
                Mempointer2=&Memory2

                if(LowerLayer.magdir):

                    Fill_rMemory_z(Mempointer1, vy,vyvy,omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)
                else:
                    Memory1.epsy=1.0+LowerLayer.cy
                    Memory1.epsz=1.0+LowerLayer.cz
                    Memory1.vz1=sqrt(1.+LowerLayer.cx-vyvy)
                    Memory1.vz2=sqrt((1.-vyvy/Memory1.epsz)*Memory1.epsy)


                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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
#                if(wavelength<23):
#                    if(th>4):
#                        if( th<5):
#                            print("tml1\n",t_ML_in_1[1][1] )

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

            Calculate_Multilayer_with_Matrices(&t_ML_back_1, &t_ML_back_2,&t_ML_in_1, &t_ML_in_2, &r_ML_in_1, &r_ML_in_2, &r_ML_back_1, &r_ML_back_2, MLREP[i]-1)
          #  print t_ML_back_2[1][1], t_ML_in_2[1][1], r_ML_in_2[1][1], r_ML_back_2[1][1]

#            if(wavelength<23):
#                if(th>4):
#                    if( th<5):
#                        print("rml\n",r_ML_in_2[1][1] )

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
                Fill_rMemory_z(&Memory1, vy,vyvy,omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)
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
                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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
                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, 0,0,0,0, \
                    LowerLayer.magdir, 0, &r, &rprime, &t, &tprime, LowerLayer.Roughness, k0)

            else:
                UpperLayer=Layer=LR[MLCOMP[i+1][0]]
                if(Memory1.IsFilled):
                    Mempointer1=&Memory1
                    Mempointer2=&Memory2
                else:
                    Mempointer2=&Memory1
                    Mempointer1=&Memory2
                Calculate_rt_z(Mempointer1, Mempointer2, vy, vyvy, omvyvy, LowerLayer.cx, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg, UpperLayer.cx, UpperLayer.cy, UpperLayer.cz, UpperLayer.cg, \
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
