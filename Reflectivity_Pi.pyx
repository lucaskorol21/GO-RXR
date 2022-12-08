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

cdef double complex Calculate_rpi_precisely(double vyvy, double complex vz1, double complex vz2, double complex cy1,double complex cy2, double complex cz1, double complex cz2):
    cdef double complex epsz1=1.+cz1
    cdef double complex epsz2=1.+cz2
    cdef double complex epsy1=1.+cy1
    cdef double complex epsy2=1.+cy2
    return ( epsz1*epsz2*(cy2-cy1)+vyvy*(cy1-cy2+cz1*epsy1-epsy2*cz2) ) *epsy1*epsy2/( epsz1*epsz2*( cquadr(vz1*epsy2+vz2*epsy1) ) )



cdef double complex LinDicParatt_Pi(Heterostructure* HS, double th, double wavelength):

    cdef double k0=6.283185307179586/wavelength
    cdef double sintheta=sin(two_pi_div_360()*th)
    cdef double vyvy=quadr(cos(two_pi_div_360()*th))

    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j

    cdef int Cap=NLAYERS-1
    cdef double complex rtot, rprime, p,t,pquad, vzupper, vzlower,rough
    #cdef double complex vzlower

    cdef double complex r_ML_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_2 # Upper reflectivity from the Multilayer
    cdef double complex ptot
    #cdef int Lower, Upper
    cdef CLayer LowerLayer, UpperLayer

    cdef double roughfak=-2.*quadr(k0)

    #rough=1.
    LowerLayer=LR[MLCOMP[0][0]]
    vzlower=CalculateVZpi(vyvy, LowerLayer.cy, LowerLayer.cz)

    if(NLAYERS==1):
        vzupper=sintheta

        rough=exp(vzlower*vzupper*quadr(LowerLayer.Roughness)*roughfak)
      #  rtot=(vzupper*(1+LR[MLCOMP[0][0]].cy)-vzlower)/(vzupper*(1+LR[MLCOMP[0][0]].cy)+vzlower)*rough

        rtot=Calculate_rpi_precisely(vyvy, vzupper, vzlower, 0,LowerLayer.cy,0,LowerLayer.cz)*rough

    else:
        #Lower=MLCOMP[0][0]
        #Upper=MLCOMP[1][0]

        UpperLayer=LR[MLCOMP[1][0]]
        vzupper=CalculateVZpi(vyvy, UpperLayer.cy, UpperLayer.cz)
        rough=exp(vzlower*vzupper*quadr(LowerLayer.Roughness )*roughfak)
       # rtot=rough*(vzupper*(1+LR[MLCOMP[0][0]].cy)-vzlower*(1+LR[MLCOMP[1][0]].cy))/(vzupper*(1+LR[MLCOMP[0][0]].cy)+vzlower*(1+LR[MLCOMP[1][0]].cy))

        rtot=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough

    i=1

    while i<NLAYERS:
        if(MLLENGTH[i]==1):
            LowerLayer=LR[MLCOMP[i][0]]
            vzlower=vzupper
            if(i!=Cap):
                #Upper=MLCOMP[i+1][0]
                UpperLayer=LR[MLCOMP[i+1][0]]
                vzupper=CalculateVZpi(vyvy, UpperLayer.cy, UpperLayer.cz)
                #Lower=MLCOMP[i][0]
                LowerLayer=LR[MLCOMP[i][0]]
                rough=exp(roughfak*vzlower*vzupper*quadr(LowerLayer.Roughness))
            #    rprime=(vzupper*LR[Lower].cy-vzlower*LR[Upper].cy)/(vzupper*LR[Lower].cy+vzlower*LR[Upper].cy)*rough
                rprime=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough
            else:
                #Lower=MLCOMP[i][0]
                LowerLayer=LR[MLCOMP[i][0]]
                rough=exp(roughfak*vzlower*sintheta*quadr(LowerLayer.Roughness))
             #   rprime=(sintheta*(LR[Lower].cy+1)-vzlower)/(sintheta*(LR[Lower].cy+1)+vzlower)*rough
                rprime=Calculate_rpi_precisely(vyvy, sintheta, vzlower, 0,LowerLayer.cy,0,LowerLayer.cz)*rough


            pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)

           # No Multiple scattering:
            rtot=(rprime+rtot*pquad)

        else:
            LowerLayer=LR[ MLCOMP[i][0] ]
            rough=exp(roughfak*vzlower*vzupper*quadr(   LowerLayer.Roughness ))

            vzlower=vzupper
            UpperLayer=LR[ MLCOMP[i][1] ]

            vzupper=CalculateVZpi(vyvy, UpperLayer.cy, UpperLayer.cz)
           # rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough
           #  r_ML_1=(vzupper*(LR[MLCOMP[i][0]].cy+1)-vzlower*(LR[MLCOMP[i][1]].cy+1))/(vzupper*(LR[MLCOMP[i][0]].cy+1)+vzlower*(LR[MLCOMP[i][1]].cy+1))*rough
          #  print r_ML_1
            r_ML_1=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough
           # print "2"

            j=1
            p=1.


            while j<MLLENGTH[i]:
              #  print "3 ", j
                vzlower=vzupper
                LowerLayer=LR[MLCOMP[i][j]]
                pquad=exp(2j*k0* LowerLayer.Thickness*vzlower)

                p*=pquad
                if(j==(MLLENGTH[i]-1)):
                    #Upper=MLCOMP[i][0]
                    UpperLayer=LR[MLCOMP[i][0]]
                else:
                    #Upper=MLCOMP[i][j+1]
                    UpperLayer=LR[MLCOMP[i][j+1]]
                vzupper=CalculateVZpi(vyvy, UpperLayer.cy, UpperLayer.cz)
               # Lower=MLCOMP[i][j]
                rough=exp(roughfak*vzlower*vzupper*quadr(  LowerLayer.Roughness))


              #  rprime=(vzupper*(LR[Lower].cy+1)-vzlower*(LR[Upper].cy+1))/(vzupper*(LR[Lower].cy+1)+vzlower*(LR[Upper].cy+1))*rough
               # print rprime
                rprime=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough
                r_ML_1=rprime+r_ML_1*pquad
              #  print "5 ", j
             #   print "new rml", r_ML_1
                j+=1
            pquad=exp(2j*k0* UpperLayer.Thickness*vzupper)
            p*=pquad

            if(MLREP[i]<=1):
                ptot=1
                r_ML_2=0
            else:
                ptot=ipow(p, MLREP[i]-1)
                r_ML_2=r_ML_1*(1.-ptot)/(1.-p)
           # print "7 "
            rtot=r_ML_2+ptot*rtot
            j=0
            while(j<MLLENGTH[i]-1):
                vzlower=vzupper
                #Lower=MLCOMP[i][j]
                #Upper=MLCOMP[i][j+1]
                LowerLayer=LR[MLCOMP[i][j]]
                UpperLayer=LR[MLCOMP[i][j+1]]
                vzupper=CalculateVZpi(vyvy, UpperLayer.cy, UpperLayer.cz)
                rough=exp(roughfak*vzlower*vzupper*quadr(  LowerLayer.Roughness )  )
                pquad=exp(2j*k0*( LowerLayer ).Thickness*vzlower)
              #  rprime=(vzupper*(LR[Lower].cy+1)-vzlower*(LR[Upper].cy+1))/(vzupper*(LR[Lower].cy+1)+vzlower*(LR[Upper].cy+1))*rough
            #    print rprime
                rprime=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough
                rtot=(rprime+rtot*pquad)
                j+=1
            #Now the next layer begins:
            vzlower=vzupper
            #Lower=MLCOMP[i][j]
            LowerLayer=LR[MLCOMP[i][j]]
            if(i!=Cap):
               # print "hier1"
               # Upper=MLCOMP[i+1][0]
                UpperLayer=LR[MLCOMP[i+1][0]]

                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
                vzupper=CalculateVZpi(vyvy, UpperLayer.cy, UpperLayer.cz)
                rough=exp(roughfak*vzlower*vzupper*quadr(  LowerLayer.Roughness )  )
              #  rprime=(vzupper*(LR[Lower].cy+1)-vzlower*(LR[Upper].cy+1))/(vzupper*(LR[Lower].cy+1)+vzlower*(LR[Upper].cy+1))*rough
              #  print rprime
                rprime=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough
                rtot=rprime+rtot*pquad
            else:
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
                rough=exp(roughfak*sintheta*vzlower*quadr(  LowerLayer.Roughness )  )
              #  rprime=(sintheta*(LR[Lower].cy+1)-vzlower)/(sintheta*(LR[Lower].cy+1)+vzlower)*rough
              #  print rprime
                rprime=Calculate_rpi_precisely(vyvy, sintheta, vzlower, 0,LowerLayer.cy,0,LowerLayer.cz)*rough
                rtot=rprime+rtot*pquad

        i=i+1
    return rtot


cdef double complex LinDicParatt_Pi_MS(Heterostructure* HS, double th, double wavelength):

    cdef double k0=6.283185307179586/wavelength
    cdef double sintheta=sin(two_pi_div_360()*th)
    cdef double vyvy=quadr(cos(two_pi_div_360()*th))


    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
    cdef double roughfak=-2.*quadr(k0)

    cdef int Cap=NLAYERS-1
    cdef double complex rtot, rprime, p,t,MSfac,pquad, vzupper, vzlower, rough
    #cdef double complex vzlower
    cdef double complex t_comp1_down #Transmission through one Compound layer, down
    cdef double complex t_comp1_up #Transmission through one Compound layer, up
    cdef double complex t_comp2_down #Transmission through one Compound layer, down
    cdef double complex t_comp2_up #Transmission through one Compound layer, up
    cdef double complex r_ML_in_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_in_2 # Upper reflectivity from the Multilayer
    cdef double complex r_ML_back_1 # Back reflectivity from one Compound layer
    cdef double complex r_ML_back_2 # Back reflectivity from the Multilayer
   # cdef int Lower, Upper
    cdef CLayer LowerLayer, UpperLayer
    LowerLayer=LR[MLCOMP[0][0]]
    vzlower=CalculateVZpi(vyvy, LowerLayer.cy, LowerLayer.cz)

    if(NLAYERS==1):
        vzupper=sintheta
        #Lower=MLCOMP[0][0]
        rough=exp(vzlower*vzupper*quadr(LowerLayer.Roughness)*roughfak)
      #  rtot=(vzupper*(1+LR[MLCOMP[0][0]].cy)-vzlower)/(vzupper*(1+LR[MLCOMP[0][0]].cy)+vzlower)*rough

        rtot=Calculate_rpi_precisely(vyvy, vzupper, vzlower, 0,LowerLayer.cy,0,LowerLayer.cz)*rough

    else:

        #Lower=MLCOMP[0][0]
        #Upper=MLCOMP[1][0]
        UpperLayer=LR[MLCOMP[1][0]]
        vzupper=CalculateVZpi(vyvy, UpperLayer.cy, UpperLayer.cz)
        rough=exp(vzlower*vzupper*quadr(LowerLayer.Roughness )*roughfak)
       # rtot=rough*(vzupper*(1+LR[MLCOMP[0][0]].cy)-vzlower*(1+LR[MLCOMP[1][0]].cy))/(vzupper*(1+LR[MLCOMP[0][0]].cy)+vzlower*(1+LR[MLCOMP[1][0]].cy))

        rtot=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough


    i=1
    while i<NLAYERS:
        if(MLLENGTH[i]==1):
            vzlower=vzupper
            if(i!=Cap):
               # Upper=MLCOMP[i+1][0]
              #  Lower=MLCOMP[i][0]
                UpperLayer=LR[MLCOMP[i+1][0]]
                LowerLayer=LR[MLCOMP[i][0]]
                vzupper=CalculateVZpi(vyvy, UpperLayer.cy, UpperLayer.cz)
                rough=exp(roughfak*vzlower*vzupper*quadr(LowerLayer.Roughness))
                #rprime=(vzupper*LR[Lower].cy-vzlower*LR[Upper].cy)/(vzupper*LR[Lower].cy+vzlower*LR[Upper].cy)*rough
                rprime=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough
            else:
                vzupper=sintheta
               # Lower=MLCOMP[i][0]
                LowerLayer=LR[MLCOMP[i][0]]
                rough=exp(roughfak*vzlower*vzupper*quadr(LowerLayer.Roughness))
                #rprime=(vzupper*LR[Lower].cy-vzlower)/(vzupper*LR[Lower].cy+vzlower)*rough
                rprime=Calculate_rpi_precisely(vyvy, vzupper, vzlower, 0,LowerLayer.cy,0,LowerLayer.cz)*rough




            pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
           # if(MultipleScattering):
            rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
           # else:
            #    rtot=rprime+rtot*pquad

        else:
            #Upper=MLCOMP[i][1]
            UpperLayer=LR[MLCOMP[i][1]]
            #Lower=MLCOMP[i][0]
            LowerLayer=LR[MLCOMP[i][0]]
            vzlower=vzupper
            vzupper=CalculateVZpi(vyvy, UpperLayer.cy, UpperLayer.cz )


            #r_ML_in_1=(vzupper*LR[MLCOMP[i][0]].cy-LR[MLCOMP[i][1]].cy*vzlower)/(vzupper*LR[MLCOMP[i][0]].cy+LR[MLCOMP[i][1]].cy*vzlower)
            r_ML_in_1=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)
            rough=exp(roughfak*vzlower*vzupper*quadr(LowerLayer.Roughness))
            r_ML_in_1*=rough

            t_comp1_up=1-r_ML_in_1
            j=1
            while j<MLLENGTH[i]:
                if(j+1<MLLENGTH[i]):
                   # Upper=MLCOMP[i][j+1]
                    UpperLayer=LR[MLCOMP[i][j+1]]
                else:
                  #  Upper=MLCOMP[i][0]
                    UpperLayer=LR[MLCOMP[i][0]]
                #Lower=MLCOMP[i][j]
                LowerLayer=LR[MLCOMP[i][j]]
                vzlower=vzupper
                vzupper=CalculateVZpi(vyvy, UpperLayer.cy , UpperLayer.cz)
                rough=exp(roughfak*vzupper*vzlower*quadr(LowerLayer.Roughness))

                rprime=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough

                p=exp(1j*k0*LowerLayer.Thickness*vzlower)
                pquad=cquadr(p)
                MSfac=1.0/(1+rprime*r_ML_in_1*pquad)
                t=1-rprime
                t_comp1_up*=p*t*MSfac
                r_ML_in_1=(rprime+r_ML_in_1*pquad)*MSfac
                j+=1
            p=exp(1j*k0*UpperLayer.Thickness*vzupper)
            t_comp1_up*=p

            r_ML_back_1=-rprime
            t_comp1_down=1-r_ML_back_1
            j=MLLENGTH[i]-2
            while j>=0:
                vzupper=vzlower
               # Upper=MLCOMP[i][j+1]
               # Lower=MLCOMP[i][j]
                UpperLayer=LR[MLCOMP[i][j+1]]
                LowerLayer=LR[MLCOMP[i][j]]
                vzlower=CalculateVZpi(vyvy, LowerLayer.cy , LowerLayer.cz)
                rough=exp(roughfak*vzupper*vzlower*quadr(LowerLayer.Roughness))


                rprime=Calculate_rpi_precisely(vyvy,  vzlower,vzupper, LowerLayer.cy,UpperLayer.cy,LowerLayer.cz,UpperLayer.cz)*rough
                t=1-rprime
                p=exp(1j*k0*UpperLayer.Thickness*vzupper)
                pquad=cquadr(p)

                MSfac=1.0/(1+rprime*r_ML_back_1*pquad)
                t_comp1_down*=p*t*MSfac

                r_ML_back_1=(rprime+r_ML_back_1*pquad)*MSfac
                j-=1
            p=exp(1j*k0*LowerLayer.Thickness*vzlower)
            t_comp1_down*=p
            r_ML_back_1*=cquadr(p)


            Calculate_Multilayer(&t_comp1_up, &t_comp2_up,&t_comp1_down, &t_comp2_down, &r_ML_in_1, &r_ML_in_2, &r_ML_back_1, &r_ML_back_2, MLREP[i]-1)



            rtot=r_ML_in_2+t_comp2_up*t_comp2_down*rtot/(1.-r_ML_back_2*rtot)
            #Now the next layer begins:
            vzupper=vzlower #CalculateVZsigma(vyvy, LR[MLCOMP[i][0]].cx)
            j=0
            while(j<MLLENGTH[i]-1):
                vzlower=vzupper
               # Lower=MLCOMP[i][j]
              #  Upper=MLCOMP[i][j+1]
                LowerLayer=LR[MLCOMP[i][j]]
                UpperLayer=LR[MLCOMP[i][j+1]]
                vzupper=CalculateVZpi(vyvy, UpperLayer.cy , UpperLayer.cz)
                rough=exp(roughfak*vzlower*vzupper*quadr(  LowerLayer.Roughness )  )
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
               # rprime=(vzupper*LR[Lower].cy-vzlower*LR[Upper].cy)/(vzupper*LR[Lower].cy+vzlower*LR[Upper].cy)*rough
                rprime=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
                j+=1
            #Now the next layer begins:
            vzlower=vzupper

          #  Lower=MLCOMP[i][j]
            LowerLayer=LR[MLCOMP[i][j]]
            if(i!=Cap):
               # print "hier1"

            #    Upper=MLCOMP[i+1][0]
                UpperLayer=LR[MLCOMP[i+1][0]]
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)

                vzupper=CalculateVZpi(vyvy, UpperLayer.cy , UpperLayer.cz)
                rough=exp(roughfak*vzlower*vzupper*quadr(  LowerLayer.Roughness )  )
                #rprime=(vzupper*LR[Lower].cy-vzlower*LR[Upper].cy)/(vzupper*LR[Lower].cy+vzlower*LR[Upper].cy)*rough
                rprime=Calculate_rpi_precisely(vyvy, vzupper, vzlower, UpperLayer.cy,LowerLayer.cy,UpperLayer.cz,LowerLayer.cz)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
            else:

                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
                rough=exp(roughfak*sintheta*vzupper*quadr(  LowerLayer.Roughness )  )
             #   rprime=(sintheta*LR[Lower].cy-vzlower)/(sintheta*LR[Lower].cy+vzlower)*rough
                rprime=Calculate_rpi_precisely(vyvy, sintheta, vzlower, 0,LowerLayer.cy,0,LowerLayer.cz)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)



        i=i+1


    return rtot
