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

cdef double complex LinDicParatt_Sigma(Heterostructure* HS, double th, double wavelength):



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
    cdef double complex rtot, rprime, vz, p,t,pquad, vzupper, vzlower
    #cdef double complex vzlower

    cdef double complex r_ML_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_2 # Upper reflectivity from the Multilayer
    cdef double complex ptot
   # cdef int Lower, Upper
    cdef CLayer LowerLayer, UpperLayer
    cdef double complex rough
    cdef double roughfak=-2.*quadr(k0)


    LowerLayer=LR[ MLCOMP[0][0] ]
    vzlower=CalculateVZsigma(vyvy, LowerLayer.cx )



  #  print NLAYERS

    if(NLAYERS==1):

        rough=exp(vzlower*sintheta*quadr(   LowerLayer.Roughness )*roughfak)
        rtot=Calculate_rsigma_precisely(sintheta, vzlower, 0 , LowerLayer.cx)*rough
    else:
        UpperLayer=LR[ MLCOMP[1][0] ]
        vzupper=CalculateVZsigma(vyvy,   UpperLayer.cx  )
        rough=exp(vzlower*vzupper*quadr(   LowerLayer.Roughness )*roughfak)
        rtot=Calculate_rsigma_precisely(vzupper, vzlower, UpperLayer.cx , LowerLayer.cx)*rough


 #   rtot=(vzupper-vzlower)/(vzupper+vzlower)*rough



  #  print rtot
    i=1
    while i<NLAYERS:

        if(MLLENGTH[i]==1):
            vzlower=vzupper
           # Lower=MLCOMP[i][0]
            LowerLayer=LR[MLCOMP[i][0]]

            if(i!=Cap):
                #Upper=MLCOMP[i+1][0]
                UpperLayer=LR[MLCOMP[i+1][0]]
                vzupper=CalculateVZsigma(vyvy,  UpperLayer.cx  )
                rough=exp(roughfak*vzlower*vzupper*quadr(   LowerLayer.Roughness  ))
                rprime=Calculate_rsigma_precisely(vzupper, vzlower, UpperLayer.cx , LowerLayer.cx)*rough
            else:
                rough=exp(roughfak*vzlower*sintheta*quadr(   LowerLayer.Roughness  ))
                rprime=Calculate_rsigma_precisely(sintheta, vzlower, 0 , LowerLayer.cx)*rough
              #  print rprime


          #  rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough


            pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
           # print pquad
          #  print pquad, rprime
           # No Multiple scattering:
            rtot=(rprime+rtot*pquad)

        else:
            UpperLayer=LR[ MLCOMP[i][1] ]
            LowerLayer=LR[ MLCOMP[i][0] ]
            rough=exp(roughfak*vzlower*vzupper*quadr(  LowerLayer.Roughness ))

            vzlower=vzupper
            vzupper=CalculateVZsigma(vyvy, UpperLayer.cx)
          #  rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough
            rprime=Calculate_rsigma_precisely(vzupper, vzlower,UpperLayer.cx , LowerLayer.cx)*rough

            r_ML_1=rprime
           # print "2"

            j=1
            p=1.

           # print( i, MLLENGTH[i] )
            while j<MLLENGTH[i]:
              #  print "3 ", j
                LowerLayer=LR[MLCOMP[i][j]]
                vzlower=vzupper
                pquad=exp(2j*k0* LowerLayer.Thickness*vzlower)
                p*=pquad
             #   print pquad, p
                if(j==(MLLENGTH[i]-1)):
                   # print("Upper", MLCOMP[i][0])
                    UpperLayer=LR[MLCOMP[i][0]]
                   # print("Thickness", UpperLayer.Thickness)
                   # print("Chi",UpperLayer.cx,UpperLayer.cy,UpperLayer.cz )
                else:
                    #Upper=MLCOMP[i][j+1]
                    UpperLayer=LR[MLCOMP[i][j+1]]
                vzupper=CalculateVZsigma(vyvy, UpperLayer.cx )

                rough=exp(roughfak*vzlower*vzupper*quadr(  LowerLayer.Roughness))

            #    rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough
                rprime=Calculate_rsigma_precisely(vzupper, vzlower, UpperLayer.cx , LowerLayer.cx)*rough

                r_ML_1=rprime+r_ML_1*pquad
              #  print "5 ", j
             #   print "new rml", r_ML_1
                j+=1
            pquad=exp(2j*k0*UpperLayer.Thickness*vzupper)
            p*=pquad
          #  print( pquad, p, UpperLayer.Thickness, vzupper, k0 )
          #  print "6 "

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
                Lower=MLCOMP[i][j]
                LowerLayer=LR[MLCOMP[i][j]]
                #Upper=MLCOMP[i][j+1]
                UpperLayer=LR[MLCOMP[i][j+1]]
                vzupper=CalculateVZsigma(vyvy,  UpperLayer.cx )
                rough=exp(roughfak*vzlower*vzupper*quadr(  LowerLayer.Roughness )  )
                pquad=exp(2j*k0*( LowerLayer ).Thickness*vzlower)
               # rprime=(vzupper-vzlower)/(vzlower+vzupper)*rough
                rprime=Calculate_rsigma_precisely(vzupper, vzlower, UpperLayer.cx , LowerLayer.cx)*rough
                rtot=(rprime+rtot*pquad)
                j+=1
            #Now the next layer begins:
            vzlower=vzupper
           # Lower=MLCOMP[i][j]
            LowerLayer=LR[MLCOMP[i][j]]
            if(i!=Cap):
               # print "hier1"

                #Upper=MLCOMP[i+1][0]
                UpperLayer=LR[MLCOMP[i+1][0]]
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)

                vzupper=CalculateVZsigma(vyvy,UpperLayer.cx  )
                rough=exp(roughfak*vzlower*vzupper*quadr(  LowerLayer.Roughness )  )
              #  rprime=(vzupper-vzlower)/(vzlower+vzupper)*rough
                rprime=Calculate_rsigma_precisely(vzupper, vzlower, UpperLayer.cx , LowerLayer.cx)*rough
                rtot=rprime+rtot*pquad
            else:

                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
                rough=exp(roughfak*sintheta*vzupper*quadr(  LowerLayer.Roughness )  )
                #rprime=(sintheta-vzlower)/(vzlower+sintheta)*rough
                rprime=Calculate_rsigma_precisely(sintheta, vzlower, 0 , LowerLayer.cx)*rough
                rtot=rprime+rtot*pquad
        i=i+1
    return rtot


cdef double complex LinDicParatt_Sigma_MS(Heterostructure* HS, double th, double wavelength):

    if(th==0):
        return 1.0
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
    cdef double complex rtot, rprime, p,t,MSfac,pquad, vzupper, vzlower,rough
    #cdef double complex vzlower
    cdef double complex t_comp1_down #Transmission through one Compound layer, down
    cdef double complex t_comp1_up #Transmission through one Compound layer, up
    cdef double complex t_comp2_down #Transmission through one Compound layer, down
    cdef double complex t_comp2_up #Transmission through one Compound layer, up
    cdef double complex r_ML_in_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_in_2 # Upper reflectivity from the Multilayer
    cdef double complex r_ML_back_1 # Back reflectivity from one Compound layer
    cdef double complex r_ML_back_2 # Back reflectivity from the Multilayer
    #cdef int Lower, Upper
    cdef CLayer LowerLayer, UpperLayer
    LowerLayer=LR[MLCOMP[0][0]]
    vzlower=CalculateVZsigma(vyvy, LowerLayer.cx)

    if(NLAYERS==1):

        rough=exp(vzlower*sintheta*quadr(   LowerLayer.Roughness )*roughfak)
        rtot=Calculate_rsigma_precisely(sintheta, vzlower, 0 , LowerLayer.cx)*rough
    else:
        UpperLayer=LR[ MLCOMP[1][0] ]
        vzupper=CalculateVZsigma(vyvy,   UpperLayer.cx  )
        rough=exp(vzlower*vzupper*quadr(   LowerLayer.Roughness )*roughfak)
        rtot=Calculate_rsigma_precisely(vzupper, vzlower,UpperLayer.cx , LowerLayer.cx)*rough





    i=1
    while i<NLAYERS:
        if(MLLENGTH[i]==1):
            vzlower=vzupper
          #  Lower=MLCOMP[i][0]
            LowerLayer=LR[MLCOMP[i][0]]
            if(i!=Cap):
                UpperLayer=LR[MLCOMP[i+1][0]]
               # Upper=MLCOMP[i+1][0]
                vzupper=CalculateVZsigma(vyvy, UpperLayer.cx)
                rough=exp(roughfak*vzlower*vzupper*quadr(LowerLayer.Roughness))
                rprime=Calculate_rsigma_precisely(vzupper, vzlower, UpperLayer.cx , LowerLayer.cx)*rough
            else:
               # vzupper=sintheta
                rough=exp(roughfak*vzlower*sintheta*quadr(LowerLayer.Roughness))
                rprime=Calculate_rsigma_precisely(sintheta, vzlower, 0 , LowerLayer.cx)*rough


            #rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough

            pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
           # if(MultipleScattering):
            rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
           # else:
            #    rtot=rprime+rtot*pquad

        else:
            vzlower=vzupper
            UpperLayer=LR[MLCOMP[i][1]]
            LowerLayer=LR[MLCOMP[i][0]]
            vzupper=CalculateVZsigma(vyvy, UpperLayer.cx)


           # r_ML_in_1=(vzupper-vzlower)/(vzupper+vzlower)
            r_ML_in_1=Calculate_rsigma_precisely(vzupper, vzlower, UpperLayer.cx , LowerLayer.cx)
            rough=exp(roughfak*vzlower*vzupper*quadr(LowerLayer.Roughness))
            r_ML_in_1*=rough

            t_comp1_up=1-r_ML_in_1

            j=1
            while j<MLLENGTH[i]:
                if(j+1<MLLENGTH[i]):
                   # Upper=MLCOMP[i][j+1]
                    UpperLayer=LR[MLCOMP[i][j+1]]
                else:
                   # Upper=MLCOMP[i][0]
                    UpperLayer=LR[MLCOMP[i][0]]
               # Lower=MLCOMP[i][j]
                LowerLayer=LR[MLCOMP[i][j]]
                vzlower=vzupper
                vzupper=CalculateVZsigma(vyvy, UpperLayer.cx)
                rough=exp(roughfak*vzupper*vzlower*quadr(LowerLayer.Roughness))
                #rprime=rough*(vzupper-vzlower)/(vzupper+vzlower)
                rprime=Calculate_rsigma_precisely(vzupper, vzlower, UpperLayer.cx , LowerLayer.cx)*rough

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
                #Upper=MLCOMP[i][j+1]
                #Lower=MLCOMP[i][j]
                UpperLayer=LR[MLCOMP[i][j+1]]
                LowerLayer=LR[MLCOMP[i][j]]
                vzlower=CalculateVZsigma(vyvy, LowerLayer.cx)
                rough=exp(roughfak*vzupper*vzlower*quadr(LowerLayer.Roughness))

                #rprime=rough*(vzlower-vzupper)/(vzupper+vzlower)
                rprime=Calculate_rsigma_precisely(vzlower, vzupper,  LowerLayer.cx , UpperLayer.cx)*rough
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
               # Upper=MLCOMP[i][j+1]
                LowerLayer=LR[MLCOMP[i][j]]
                UpperLayer=LR[MLCOMP[i][j+1]]
                vzupper=CalculateVZsigma(vyvy,  UpperLayer.cx )
                rough=exp(roughfak*vzlower*vzupper*quadr(  LowerLayer.Roughness )  )
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
                #rprime=(vzupper-vzlower)/(vzlower+vzupper)*rough
                rprime=Calculate_rsigma_precisely(vzupper, vzlower, UpperLayer.cx , LowerLayer.cx)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
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

                vzupper=CalculateVZsigma(vyvy,UpperLayer.cx  )
                rough=exp(roughfak*vzlower*vzupper*quadr( LowerLayer.Roughness )  )
                #rprime=(vzupper-vzlower)/(vzlower+vzupper)*rough
                rprime=Calculate_rsigma_precisely(vzupper, vzlower, UpperLayer.cx , LowerLayer.cx)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
            else:

                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
                rough=exp(roughfak*sintheta*vzlower*quadr(  LowerLayer.Roughness )  )
               # rprime=(sintheta-vzlower)/(vzlower+sintheta)*rough
                rprime=Calculate_rsigma_precisely(sintheta, vzlower, 0 , LowerLayer.cx)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)


        i=i+1


    return rtot


