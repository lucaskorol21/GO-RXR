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


cdef void Relevant_Stuff_for_xmag_precisely(double complex cy1, double complex cy2, double complex cz1, double complex cz2, \
                                  double complex cg1, double complex cg2, double complex *vz1, double complex *vz2, \
                                  double vy, double vyvy, double k0, double sigma, int IsMagnetic1, int IsMagnetic2, \
                                  double complex *r, double complex *rp, double complex *t, double complex *tp):
    cdef double complex epsz1=1.+cz1
    cdef double complex epsz2=1.+cz2
    cdef double complex epsy1=1.+cy1
    cdef double complex epsy2=1.+cy2

    cdef double complex difvz
    cdef double complex cgcgz1
    cdef double complex cgcgz2
    cdef double complex term_a
    cdef double complex term_b
    cdef double complex term_c
    cdef double complex term_d
    cdef double complex term_e

    cdef double complex phidif1
    cdef double complex phidif2
    cdef double complex J11,J12,J21,J22
    cdef double roughfac=-0.5*quadr(sigma)*quadr(k0)
    cdef double complex roughplus
    cdef double complex roughminus

    if( IsMagnetic1&IsMagnetic2 ):
        cgcgz1=cg1*cg1/epsz1
        cgcgz2=cg2*cg2/epsz2
        vz2[0]=sqrt((1.-vyvy/epsz2)*epsy2+cgcgz2)

        difvz=( ( (epsz1-vyvy)*(cy1+cz2*epsy1) - (epsz2-vyvy)*(cy2+cz1*epsy2)+(cz1-cz2) ) -cg2*cg2*epsz1+cg1*cg1*epsz2  )/( (vz1[0]+vz2[0])*epsz1*epsz2  )
        term_a=cy1+cgcgz1
        term_b=cy2+cgcgz2
        term_c=difvz+term_b*vz1[0] - term_a*vz2[0]
        term_d=-(1+term_b)*vy*cg1/epsz1 + (1+term_a)*vy*cg2/epsz2
        term_e=( (1+term_a)*(1+term_b) )
        phidif1=( term_c +term_d )/term_e
        phidif2=( term_c -term_d )/term_e

        prefac=(1+term_a)/(2*vz1[0])
        roughplus=exp(roughfac*cquadr(vz1[0]+vz2[0]))
        roughminus=exp(roughfac*cquadr(vz1[0]-vz2[0]))
        term_d=-vy*cg1/epsz1
        term_e=-vy*cg2/epsz2

        J11=roughminus*( (vz1[0]-term_d)/(1+term_a) + (term_e+vz2[0])/(1+term_b) )
        J22=roughminus*( (vz1[0]+term_d)/(1+term_a) - (term_e-vz2[0])/(1+term_b) )
        J12=roughplus*phidif2
        J21=roughplus*phidif1
        r[0]=(-J12/J11)
        t[0]=J21*r[0]+J22
        rp[0]=J21/J11
        tp[0]=1./J11

    elif(IsMagnetic1):
        cgcgz1=cg1*cg1/epsz1
        vz2[0]=sqrt((1.-vyvy/epsz2)*epsy2)

        difvz=( ( (epsz1-vyvy)*(cy1+cz2*epsy1) - (epsz2-vyvy)*(cy2+cz1*epsy2)+(cz1-cz2) ) +cg1*cg1*epsz2  )/( (vz1[0]+vz2[0])*epsz1*epsz2  )
        term_a=cy1+cgcgz1

        term_c=difvz+cy2*vz1[0] - term_a*vz2[0]
        term_d=-(1+cy2)*vy*cg1/epsz1
        term_e=( (1+term_a)*(1+cy2) )
        phidif1=( term_c +term_d )/term_e
        phidif2=( term_c -term_d )/term_e


        roughplus=exp(roughfac*cquadr(vz1[0]+vz2[0]))
        roughminus=exp(roughfac*cquadr(vz1[0]-vz2[0]))
        term_d=-vy*cg1/epsz1

        J11=roughminus*( (vz1[0]-term_d)/(1+term_a) + vz2[0]/(1+cy2) )
        J22=roughminus*( (vz1[0]+term_d)/(1+term_a) + vz2[0]/(1+cy2) )
        J12=roughplus*phidif2
        J21=roughplus*phidif1
        r[0]=(-J12/J11)
        t[0]=J21*r[0]+J22
        rp[0]=J21/J11
        tp[0]=1./J11
    elif(IsMagnetic2):

        cgcgz2=cg2*cg2/epsz2
        vz2[0]=sqrt((1.-vyvy/epsz2)*epsy2+cgcgz2)

        difvz=( ( (epsz1-vyvy)*(cy1+cz2*epsy1) - (epsz2-vyvy)*(cy2+cz1*epsy2)+(cz1-cz2) ) -cg2*cg2*epsz1 )/( (vz1[0]+vz2[0])*epsz1*epsz2  )

        term_b=cy2+cgcgz2
        term_c=difvz+term_b*vz1[0] - cy1*vz2[0]
        term_d=(  1+cy1)*vy*cg2/epsz2
        term_e=( (1+cy1)*(1+term_b) )
        phidif1=( term_c +term_d )/term_e
        phidif2=( term_c -term_d )/term_e

        roughplus=exp(roughfac*cquadr(vz1[0]+vz2[0]))
        roughminus=exp(roughfac*cquadr(vz1[0]-vz2[0]))

        term_e=-vy*cg2/epsz2

        J11=roughminus*( vz1[0]/(1+cy1) + (term_e+vz2[0])/(1+term_b) )
        J22=roughminus*( vz1[0]/(1+cy1) - (term_e-vz2[0])/(1+term_b) )
        J12=roughplus*phidif2
        J21=roughplus*phidif1
        r[0]=(-J12/J11)
        t[0]=J21*r[0]+J22
        rp[0]=J21/J11
        tp[0]=1./J11
    else:

        vz2[0]=sqrt((1.-vyvy/epsz2)*epsy2)

        difvz=( ( (epsz1-vyvy)*(cy1+cz2*epsy1) - (epsz2-vyvy)*(cy2+cz1*epsy2)+(cz1-cz2) ) )/( (vz1[0]+vz2[0])*epsz1*epsz2  )

        term_c=difvz+cy2*vz1[0] - cy1*vz2[0]

        term_e=( (1+cy1)*(1+cy2) )
        phidif1=term_c/term_e

        roughplus=exp(roughfac*cquadr(vz1[0]+vz2[0]))
        roughminus=exp(roughfac*cquadr(vz1[0]-vz2[0]))

        J11=roughminus*( vz1[0]/(1+cy1) +vz2[0]/(1+cy2) )

        J12=roughplus*phidif1

        r[0]=(-J12/J11)
        t[0]=J12*r[0]+J11
        rp[0]=J12/J11
        tp[0]=1./J11

cdef double complex LinDicParatt_Pi_xmag(Heterostructure* HS, double th, double wavelength):

    cdef double k0=6.283185307179586/wavelength
    cdef double sintheta=sin(two_pi_div_360()*th)
    cdef double vy=cos(two_pi_div_360()*th)
    cdef double vyvy=quadr(vy)



    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j

    cdef int Cap=NLAYERS-1
    cdef double complex rtot, r, rprime, p,t,tp,ttot,pquad, vzupper, vzlower
    #cdef double complex vzlower

    cdef double complex r_ML_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_2 # Upper reflectivity from the Multilayer
    cdef double complex ptot
   # cdef int Lower, Upper
    cdef CLayer LowerLayer, UpperLayer
   # cdef double complex rough, groundrough
  #  cdef double roughfak=-2.*quadr(k0)

    #rough=1.
    LowerLayer=LR[MLCOMP[0][0]]
    vzlower=CalculateVZpi_m(vyvy, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)
    if(NLAYERS==1):
        #vzupper=sintheta
       # rough=exp(vzlower*vzupper*quadr(LR[MLCOMP[0][0]].Roughness)*roughfak)
#        Relevant_Stuff_for_xmag(LR[MLCOMP[0][0]].ey, 1., LR[MLCOMP[0][0]].ez, 1, \
#                                LR[MLCOMP[0][0]].eg, 0., vzlower, sintheta, \
#                                vy, k0, LR[MLCOMP[0][0]].Roughness, \
#                                &rtot, &rprime, &t, &tp)

        Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, 0, LowerLayer.cz, 0, \
                                  LowerLayer.cg, 0, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, 0, \
                                  &rtot, &rprime, &t, &tp)

    #    print rtot



    else:
       # vzupper=CalculateVZpi(vyvy, LR[MLCOMP[1][0]].ey, LR[MLCOMP[1][0]].ez)
     #   rough=exp(vzlower*vzupper*quadr(LR[MLCOMP[0][0]].Roughness)*roughfak)
#        Relevant_Stuff_for_xmag(LR[MLCOMP[0][0]].ey, LR[MLCOMP[1][0]].ey, LR[MLCOMP[0][0]].ez, LR[MLCOMP[1][0]].ez, \
#                                LR[MLCOMP[0][0]].eg, LR[MLCOMP[1][0]].eg, vzlower, vzupper, \
#                                vy, k0, LR[MLCOMP[0][0]].Roughness, \
#                                &rtot, &rprime, &t, &tp)

        UpperLayer=LR[MLCOMP[1][0]]
        Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &rtot, &rprime, &t, &tp)

    #rtot=rtot*rough


    i=1
    while i<NLAYERS:

        if(MLLENGTH[i]==1):
            vzlower=vzupper
            if(i!=Cap):
                LowerLayer=LR[MLCOMP[i][0]]
                UpperLayer=LR[MLCOMP[i+1][0]]
                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r, &rprime, &t, &tp)
              #  print "endif"
            else:
                LowerLayer=LR[MLCOMP[i][0]]
                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, 0, LowerLayer.cz, 0, \
                                  LowerLayer.cg, 0, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, 0, \
                                  &r, &rprime, &t, &tp)

            pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)

            rtot=r+t*tp*pquad*rtot

        else:
           # print "going into ml"

            vzlower=vzupper
#            Upper=MLCOMP[i][1]
#            Lower=MLCOMP[i][0]
#            vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
#
#            Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
#                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper,  \
#                                vy, k0, LR[Lower].Roughness, \
#                                &r_ML_1, &rprime, &t, &tp)
            LowerLayer=LR[MLCOMP[i][0]]
            UpperLayer=LR[MLCOMP[i][1]]
            Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r_ML_1, &rprime, &t, &tp)
            j=1
            p=1.
            ttot=t*tp

            while j<MLLENGTH[i]:
              #  print "3 ", j
                LowerLayer=LR[MLCOMP[i][j]]
                vzlower=vzupper
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)

                p*=pquad
                if(j==(MLLENGTH[i]-1)):
#                    Lower=MLCOMP[i][j]
#                    Upper=MLCOMP[i][0]
#                    vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
#
#                    Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
#                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
#                                vy, k0, LR[Lower].Roughness, \
#                                &r, &rprime, &t, &tp)

                    UpperLayer=LR[MLCOMP[i][0]]
                    Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r, &rprime, &t, &tp)
                else:
#                    Lower=MLCOMP[i][j]
#                    Upper=MLCOMP[i][j+1]
#                    vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
#
#                    Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
#                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
#                                vy, k0, LR[Lower].Roughness, \
#                                &r, &rprime, &t, &tp)

                    UpperLayer=LR[MLCOMP[i][j+1]]
                    Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r, &rprime, &t, &tp)
              #  r=r*rough
                ttot*=t*tp
                r_ML_1=r+t*tp*r_ML_1*pquad

                j=j+1
            pquad=exp(2j*k0*UpperLayer.Thickness*vzupper)
            p*=pquad
          #  print "6 "
            if(MLREP[i]<=1):
                ptot=1
                r_ML_2=0
            else:
                ptot=ipow(p*ttot, MLREP[i]-1)

                r_ML_2=r_ML_1*(1.-ptot)/(1.-p*ttot)
            #ttot=ipow(, MLREP[i])
          #  r_ML_2=r_ML_1*(1.-ptot)/(1.-p*ttot)
           # print "7 "
            #print "rml2", r_ML_2

            rtot=r_ML_2+ptot*rtot




            j=0
            while(j<MLLENGTH[i]-1):
                vzlower=vzupper
#                Lower=MLCOMP[i][j]
#                Upper=MLCOMP[i][j+1]
#                vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
#                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
#                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
#                                vy, k0, LR[Lower].Roughness, \
#                                &r, &rprime, &t, &tp)
                LowerLayer=LR[MLCOMP[i][j]]
                UpperLayer=LR[MLCOMP[i][j+1]]
                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r, &rprime, &t, &tp)
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
                rtot=(r+rtot*pquad)
                j+=1
            #Now the next layer begins:
            vzlower=vzupper

            if(i!=Cap):
               # print "hier1"
                #Lower=MLCOMP[i][j]
            #    Upper=MLCOMP[i+1][0]
                LowerLayer=LR[MLCOMP[i][j]]
                UpperLayer=LR[MLCOMP[i+1][0]]
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)

              #  vzupper=CalculateVZpi_m(vyvy, UpperLayer.cy, LR[Upper].ez, LR[Upper].eg)
#                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
#                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
#                                vy, k0, LR[Lower].Roughness, \
#                                &r, &rprime, &t, &tp)

                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r, &rprime, &t, &tp)
                rtot=r+rtot*pquad
            else:
             #   Lower=MLCOMP[i][j]
                LowerLayer=LR[MLCOMP[i][j]]
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
#                Relevant_Stuff_for_xmag(LR[Lower].ey, 1, LR[Lower].ez, 1, \
#                                LR[Lower].eg, 0, vzlower, sintheta, \
#                                vy, k0, LR[Lower].Roughness, \
#                                &r, &rprime, &t, &tp)
                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, 0, LowerLayer.cz, 0, \
                                  LowerLayer.cg, 0, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, 0, \
                                  &r, &rprime, &t, &tp)
                rtot=r+rtot*pquad

        i=i+1
   # print "8 "

    return rtot





cdef double complex LinDicParatt_Pi_xmag_MS(Heterostructure* HS, double th, double wavelength):


    cdef double k0=6.283185307179586/wavelength
    cdef double sintheta=sin(two_pi_div_360()*th)
    cdef double vy=cos(two_pi_div_360()*th)
    cdef double vyvy=quadr(vy)


    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
   # cdef double roughfak=-2.*quadr(k0)
    cdef double complex rough, rough2
    cdef int Cap=NLAYERS-1
    cdef double complex rtot, rprime, vz,r, p,t,tp,MSfac,pquad, vzupper, vzlower
    #cdef double complex vzlower
    cdef double complex t_ML_in_1 #Transmission through one Compound layer, down
    cdef double complex t_ML_back_1 #Transmission through one Compound layer, up
    cdef double complex t_ML_in_2 #Transmission through all Compound layer, down
    cdef double complex t_ML_back_2 #Transmission through all Compound layer, up
    cdef double complex r_ML_in_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_in_2 # Upper reflectivity from the Multilayer
    cdef double complex r_ML_back_1 # Back reflectivity from one Compound layer
    cdef double complex r_ML_back_2 # Back reflectivity from the Multilayer
    cdef CLayer LowerLayer, UpperLayer
   # cdef double complex rough, groundrough
  #  cdef double roughfak=-2.*quadr(k0)

    #rough=1.
    LowerLayer=LR[MLCOMP[0][0]]
    vzlower=CalculateVZpi_m(vyvy, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)

    if(NLAYERS==1):
        Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, 0, LowerLayer.cz, 0, \
                                  LowerLayer.cg, 0, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, 0, \
                                  &rtot, &rprime, &t, &tp)
    else:
        UpperLayer=LR[MLCOMP[1][0]]
        Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &rtot, &rprime, &t, &tp)
    i=1
    while i<NLAYERS:

        if(MLLENGTH[i]==1):
            vzlower=vzupper
            if(i!=Cap):
                LowerLayer=LR[MLCOMP[i][0]]
                UpperLayer=LR[MLCOMP[i+1][0]]
                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r, &rprime, &t, &tp)
              #  print "endif"
            else:
                LowerLayer=LR[MLCOMP[i][0]]
                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, 0, LowerLayer.cz, 0, \
                                  LowerLayer.cg, 0, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, 0, \
                                  &r, &rprime, &t, &tp)

            pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)
            rtot=r+t*tp*pquad*rtot/(1.-rtot*rprime*pquad)
          #  print rtot, vzlower, vzupper
        else:

            vzlower=vzupper

            LowerLayer=LR[ MLCOMP[i][0] ]
            UpperLayer=LR[ MLCOMP[i][1] ]

#            Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
#                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
#                                vy, k0, LR[Lower].Roughness, \
#                                &r_ML_in_1, &rprime, &t, &t_ML_back_1)
            Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r_ML_in_1, &rprime, &t, &t_ML_back_1)

            j=1
            while j<MLLENGTH[i]:
                if(j+1<MLLENGTH[i]):
                    UpperLayer=LR[ MLCOMP[i][j+1] ]
                else:
                    UpperLayer=LR[ MLCOMP[i][0] ]
                LowerLayer=LR[ MLCOMP[i][j] ]

                vzlower=vzupper
             #   vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
#                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
#                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
#                                vy, k0, LR[Lower].Roughness, \
#                                &r, &rprime, &t, &tp)
                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r, &rprime, &t, &tp)
                p=exp(1j*k0*LowerLayer.Thickness*vzlower)
                pquad=cquadr(p)


                MSfac=1.0/(1.-rprime*r_ML_in_1*pquad)

                t_ML_back_1=t_ML_back_1*p*tp*MSfac
                r_ML_in_1=r+t*tp*pquad*r_ML_in_1*MSfac

                j+=1
            p=exp(1j*k0*UpperLayer.Thickness*vzupper)
            t_ML_back_1*=p


            r_ML_back_1=rprime
            t_ML_in_1=t




            j=MLLENGTH[i]-2
            while j>=0:
              #  vzupper=vzlower
                UpperLayer=LR[MLCOMP[i][j+1]]
                LowerLayer=LR[MLCOMP[i][j]]
                vzlower=CalculateVZpi_m(vyvy, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)
#                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
#                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
#                                vy, k0, LR[Lower].Roughness, \
#                                &r, &rprime, &t, &tp)
                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r, &rprime, &t, &tp)
                p=exp(1j*k0*UpperLayer.Thickness*vzupper)
                pquad=cquadr(p)

                MSfac=1.0/(1.-r*r_ML_back_1*pquad)

                t_ML_in_1*=t*p*MSfac
                r_ML_back_1=rprime+t*tp*pquad*r_ML_back_1*MSfac
                j-=1

            p=exp(1j*k0*LowerLayer.Thickness*vzlower)
            t_ML_in_1*=p
            r_ML_back_1*=cquadr(p)



            Calculate_Multilayer(&t_ML_back_1, &t_ML_back_2,&t_ML_in_1, &t_ML_in_2, &r_ML_in_1, &r_ML_in_2, &r_ML_back_1, &r_ML_back_2, MLREP[i]-1)




            rtot=r_ML_in_2+t_ML_back_2*t_ML_in_2*rtot/(1.-r_ML_back_2*rtot)




        #    Bis hier scheints zu stimmen...
            #Now the next layer begins:

            LowerLayer=LR[ MLCOMP[i][0] ]
            vzlower=CalculateVZpi_m(vyvy, LowerLayer.cy, LowerLayer.cz, LowerLayer.cg)

            vzupper=vzlower #CalculateVZsigma(vyvy, LR[MLCOMP[i][0]].ex)
            j=0
            while(j<MLLENGTH[i]-1):
                vzlower=vzupper
                UpperLayer=LR[MLCOMP[i][j+1]]
                LowerLayer=LR[MLCOMP[i][j]]

               # vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
#                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
#                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
#                                vy, k0, LR[Lower].Roughness, \
#                                &r, &rprime, &t, &tp)
                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r, &rprime, &t, &tp)
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)

                rtot=r+t*tp*pquad*rtot/(1.-rtot*rprime*pquad)
#                if(wavelength<23):
#                    if(th>40):
#                        if( th<41):
#                            print("rtot",j+1, rtot, r, t*tp, pquad)


             #   print rtot, vzlower, vzupper
                j+=1
            #Now the next layer begins:
            vzlower=vzupper
            LowerLayer=LR[MLCOMP[i][j]]
            if(i!=Cap):
               # print "hier1"

                UpperLayer=LR[MLCOMP[i+1][0]]
                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)

            #    vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, UpperLayer.cy, LowerLayer.cz, UpperLayer.cz, \
                                  LowerLayer.cg, UpperLayer.cg, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, UpperLayer.magdir, \
                                  &r, &rprime, &t, &tp)
                rtot=r+t*tp*pquad*rtot/(1.-rtot*rprime*pquad)
            else:

                pquad=exp(2j*k0*LowerLayer.Thickness*vzlower)

                Relevant_Stuff_for_xmag_precisely(LowerLayer.cy, 0, LowerLayer.cz, 0, \
                                  LowerLayer.cg, 0, &vzlower, &vzupper, \
                                  vy, vyvy, k0, LowerLayer.Roughness, LowerLayer.magdir, 0, \
                                  &r, &rprime, &t, &tp)
                rtot=r+t*tp*pquad*rtot/(1.-rtot*rprime*pquad)
#            if(wavelength<23):
#                if(th>40):
#                    if( th<41):
#                        print("rtot end", rtot)


     #print rtot, vzlower, vzupper
        i=i+1

    return rtot


