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




cdef void Reduce_complexity_of_chi(CLayer *Layer, double Cutoff, int *allx, int *ally, int *allz):
    cdef int nox, noy, noz
    cdef int gyro


    gyro=1
    if( Cmaxnorm((Layer[0]).cxy)< Cutoff ):
        (Layer[0]).cxy=0
    if( Cmaxnorm((Layer[0]).cyx)< Cutoff ):
        (Layer[0]).cyx=0
    if( Cmaxnorm((Layer[0]).cxz)< Cutoff ):
        (Layer[0]).cxz=0
    if( Cmaxnorm((Layer[0]).czx)< Cutoff ):
        (Layer[0]).czx=0
    if( Cmaxnorm((Layer[0]).cyz)< Cutoff ):
        (Layer[0]).cyz=0
    if( Cmaxnorm((Layer[0]).czy)< Cutoff ):
        (Layer[0]).czy=0

    if( ((Layer[0]).cxy == 0 ) and ((Layer[0]).cyx==0) ):
        noz=1
    else:
       # print("hallo, ", (Layer[0]).cxy, (Layer[0]).cyx )
        noz=0
        if( ((Layer[0]).cxy == -(Layer[0]).cyx ) ):
            gyro=1
        else:
            gyro=0
    if( ((Layer[0]).cxz == 0 ) and ((Layer[0]).czx==0) ):
        noy=1
    else:
        noy=0
        if( ((Layer[0]).cxz == -(Layer[0]).czx ) ):
            gyro=1
        else:
            gyro=0
    if( ((Layer[0]).cyz == 0 ) and ((Layer[0]).czy==0) ):
        nox=1
    else:
        nox=0
        if( ((Layer[0]).cyz == -(Layer[0]).czy ) ):
            gyro=1
        else:
            gyro=0


  #  print( "checker:", nox, noy, noz, gyro)
  #  print( "checker2:", (Layer[0]).cxy, (Layer[0]).cyx, (Layer[0]).cxz, (Layer[0]).czx, (Layer[0]).cyz, (Layer[0]).czy)
    if( (nox+noy+noz)<2 or gyro==0 ):
     #   print("Full chi", nox, noy, noz, gyro)
        return
    else:
        if( nox and noy and noz ):
            (Layer[0]).type=2
            (Layer[0]).magdir=0
        #    print( "diagonal chi")
            return
        if( nox and noy ):
            (Layer[0]).cg=(Layer[0]).cxy
            (Layer[0]).type=3
            (Layer[0]).magdir=3
            allz[0]=1
            #print( "Polar Kerr")
            return
        if( nox and noz ):
            (Layer[0]).cg=(Layer[0]).cxz
            (Layer[0]).type=3
            (Layer[0]).magdir=2
            ally[0]=1
           # print( "Longitudinal Ker")
            return
        if( noy and noz ):
            (Layer[0]).cg=(Layer[0]).czy
            (Layer[0]).type=3
            (Layer[0]).magdir=1
            allx[0]=1
           # print( "Transversal Kerr")
            return


cdef void Fill_Matrixsafer( MatrixSafer *MS, CLayer L ):

    cdef double complex cx, cy, cz, cxy, cyx, cxz, czx, cyz, czy, cg
    cx=L.cx
    cy=L.cy
    cz=L.cz
    cxy=L.cxy
    cyx=L.cyx
    cxz=L.cxz
    czx=L.czx
    cyz=L.cyz
    czy=L.czy

   # print("hi")
    (MS[0]).IsFilled=1
    (MS[0]).exx=cx+1.
    (MS[0]).eyy=cy+1.
    (MS[0]).ezz=cz+1.

    (MS[0]).Mx=czy+cyz
    (MS[0]).exyyx=cxy*cyx
    (MS[0]).exzzx=cxz*czx
    (MS[0]).eyzzy=cyz*czy
    (MS[0]).crossmag=cxy*cyz*czx + cxz*cyx*czy



    (MS[0]).summag=(MS[0]).exzzx + (MS[0]).eyzzy
    (MS[0]).mixmag=cxy*czx + cxz*cyx

    if( (cxy==cyx) and (cxz==-czx) ):
        (MS[0]).mixmag=0
    elif( (cxy==-cyx) and (cxz==czx) ):
        (MS[0]).mixmag=0

    (MS[0]).inverseezz = 1.0/(MS[0]).ezz

    (MS[0]).D21ic=(MS[0]).exx-(MS[0]).exzzx*(MS[0]).inverseezz
    (MS[0]).D23=cxy-cxz*czy*(MS[0]).inverseezz
    (MS[0]).D24ic=-cxz*(MS[0]).inverseezz

    (MS[0]).D31ic=-czx*(MS[0]).inverseezz
    (MS[0]).D33ic=-czy*(MS[0]).inverseezz


    (MS[0]).D41=cyx-cyz*czx*(MS[0]).inverseezz
    (MS[0]).D43=(MS[0]).eyy-(MS[0]).eyzzy*(MS[0]).inverseezz
    (MS[0]).D44ic=-cyz*(MS[0]).inverseezz


cdef void NormalizePHI(double complex (*PHI)[4][4] ):

    cdef int i, j
    cdef double sum1

    for i in range(4):
        sum1=0
        for j in range(4):
            sum1+=cabsquadr( (PHI[0])[j][i] )
        sum1=dsqrt(sum1)
        for j in range(4):
            (PHI[0])[j][i]/=sum1


cdef void NormalizePSI(double complex (*PSI)[4][4] ):

    cdef int i, j
    cdef double sum1

    for i in range(4):
        sum1=0
        for j in range(4):
            sum1+=cabsquadr( (PSI[0])[i][j] )
        sum1=dsqrt(sum1)
        for j in range(4):
            (PSI[0])[i][j]/=sum1


cdef void Calculate_Phi_and_Psi(CLayer L, MatrixSafer *MS, double vy, double vzvz, double vyvy, double complex (*vz)[4], double complex (*PHI)[4][4], double complex (*PSI)[4][4]):

    cdef double complex cx, cy, cz, cxy, cyx, cxz, czx, cyz, czy, cg
    cx=L.cx
    cy=L.cy
    cz=L.cz
    cxy=L.cxy
    cyx=L.cyx
    cxz=L.cxz
    czx=L.czx
    cyz=L.cyz
    czy=L.czy
    cg=L.cg
    cdef double complex a,b,c,d # For the forth order polynomial
    cdef double complex bu, bm, du, dm
    cdef double complex vsigmasquared, vpisquared, vsigma, vpi, vpisquared_special1, vpisquared_special2
    cdef double complex one_minus_vyvy_div_ezz, eyzzy_div_ezz
    cdef double complex xld, xld_special1, xld_special2
    cdef double complex aquadr, bquadr, dquadr
    cdef double complex bd_splitter
    cdef double complex min_a_div_4
    cdef double complex ri_bd, ri_ac, ri, root
    cdef double complex Delta0, Delta0_ac, Delta0_bd, Delta1, Delta1_ac, Delta1_bd
    cdef double complex Qinternal, Q, p, qu, qm, q, Qsquared,Qrest, q_div_S
    cdef double complex pm, S, F1, F2, M
    cdef double complex vrem[4]

    cdef double complex D21_min_SS
    cdef double complex D21_min_vzsq[4]
    cdef double complex D34D43_min_vD33vD44[4]
    cdef double complex phipre[4]


    cdef double complex rel34u[4]
    cdef double complex rel34m[4]
    cdef double complex rel34[4]

    cdef double complex rel13

    cdef double complex bs1, bs2, bs3, bs4, bs5, bs6, bs7, bs8
    cdef double complex denom
    cdef double comp1, comp2, comp3, comp4, comp5, comp6
    cdef double maxdiag
    cdef double complex D21, D31, D33, D44, D24
    cdef int use_vz, use_eq, maxcon, mincon
    cdef double complex testzero, testzerosq, small_c_correction, small_c_error, small_a_correction, small_a_error
    cdef int use_perturbative_aczero_case
    cdef int switch
    cdef double complex xld_min=1.0e-25
    ###DEBUG
    cdef double evch[4][4]
    cdef double complex mag_vs_xld

    if( L.type==4 ):
        eyzzy_div_ezz = (MS[0]).eyzzy*(MS[0]).inverseezz

        xld = ( (cy-cx) + vyvy*(cz-cy)*(MS[0]).inverseezz ) - eyzzy_div_ezz



        one_minus_vyvy_div_ezz=(1.-vyvy*(MS[0]).inverseezz)


        a=vy*(MS[0]).Mx*(MS[0]).inverseezz


        vsigmasquared = vzvz +cx
     #   print("start")

      #  vsigmasquared = (MS[0]).exx -vyvy

        vpisquared = (cz + vzvz)*(MS[0]).eyy*(MS[0]).inverseezz-eyzzy_div_ezz
       # print( vpisquared )
       # vpisquared = one_minus_vyvy_div_ezz*(MS[0]).eyy-eyzzy_div_ezz
       # print( vpisquared )
        bu = -vsigmasquared-vpisquared
        bm= +(MS[0]).exzzx*(MS[0]).inverseezz


        b= bu+bm

        c= -a*vsigmasquared + vy*(MS[0]).mixmag*(MS[0]).inverseezz

        du=vpisquared*vsigmasquared
        dm=  - (MS[0]).eyy*bm - (MS[0]).exyyx*one_minus_vyvy_div_ezz +(MS[0]).crossmag*(MS[0]).inverseezz
        d=du+dm


        bquadr=b*b
        aquadr=a*a
        dquadr=d*d
        min_a_div_4 = -a/4
        D21= (MS[0]).D21ic - vyvy
        D24=vy*(MS[0]).D24ic
        D31=vy*(MS[0]).D31ic
        D33=vy*(MS[0]).D33ic
        D44=vy*(MS[0]).D44ic

#        print("start of pp")
#        print(cyz, czy)
#        print(a)
#        print(b)
#        print(c)
#        print(d)



#
#        for i in range(4):
#            print(i,  (vz[0])[i] )
        maxdiag = Cmaxnorm(cx)
        comp6=Cmaxnorm(  vpisquared-(1-vyvy)  )
        if( comp6>maxdiag ):
            maxdiag=comp6

        if(xld==0):
            xld=maxdiag*xld_min

        mag_vs_xld=sqrt( 0.25*xld*xld + 0.25*bm*(2*bu + bm) - dm  )

        if( a==0 and c==0):
            use_perturbative_aczero_case=1
            small_c_correction=0
            small_a_correction=0
        else:
            bs1=mag_vs_xld
            if(bs1==0):
                use_perturbative_aczero_case=0
            else:
                testzerosq = -0.5*b + bs1
                testzero=sqrt(testzerosq)
                small_c_correction = -1./(4*bs1)
                small_a_correction = testzerosq*small_c_correction
                small_c_error = c*( small_c_correction*(6*testzerosq*small_c_correction + b*small_c_correction + 1) )
                small_a_error = a*( small_a_correction*(6*testzerosq*small_a_correction + b*small_a_correction + 3*testzerosq ) )
              #  av**2*deva *( (6*x01*x01+b)*deva + 3*x01*x01 )

                if( Cmaxnorm( (c*small_c_error + a*small_a_error)/testzero)  < 2.0e-16*maxdiag ):
                    use_perturbative_aczero_case=1
                else:
                    use_perturbative_aczero_case=0

#        if(use_perturbative_aczero_case):
#            print("use perturbative aczero")
#        else:
#            print("no perturbation")
#            print(a,c)
#            print( cyz, czy , bs1)



        if( ( (D31==0) and ( (MS[0]).D41==0) ) and ( ( (MS[0]).D23==0 )  and (D24==0) ) ):
        #    print("1")
            (vz[0])[0]=sqrt(vsigmasquared)
            (vz[0])[1]=( sqrt( -min_a_div_4*a+vpisquared ) -a/2 )
            (vz[0])[2]=-(vz[0])[0]
            (vz[0])[3]=-( sqrt( -min_a_div_4*a+vpisquared ) +a/2 )

            (PHI[0])[0][0] = 1
            (PHI[0])[1][0] = (vz[0])[0]
            (PHI[0])[2][0] = 0
            (PHI[0])[3][0] = 0

            (PHI[0])[2][1] = 1
            (PHI[0])[3][1] = ( (vz[0])[1] - D33 )/one_minus_vyvy_div_ezz
            (PHI[0])[0][1] = 0
            (PHI[0])[1][1] = 0


            (PHI[0])[0][2] = 1
            (PHI[0])[1][2] = (vz[0])[2]
            (PHI[0])[2][2] = 0
            (PHI[0])[3][2] = 0

            (PHI[0])[2][3] = 1
            (PHI[0])[3][3] = ( (vz[0])[3] -D33 )/one_minus_vyvy_div_ezz
            (PHI[0])[0][3] = 0
            (PHI[0])[1][3] = 0
            NormalizePHI(PHI)

            Matrix4Invert(PHI, PSI)


        else:
            if(use_perturbative_aczero_case):
#                print("a and c", a, c)
#                print("hallo")



               # bs1=sqrt( 0.25*xld*xld + 0.25*bm*(2*bu + bm) - dm  )
                bs1=mag_vs_xld
                bs2=( + 0.5*( xld +bm ) - bs1   )


                bs4= 0.5*( xld +bm ) + bs1

                if(bs2==0):
                    bs2= ( 0.5*bm*(xld-bu)-dm )/bs4
                if( Cmaxnorm(dm)<2e-16*maxdiag and Cmaxnorm(bm)<2e-16*maxdiag ):
                  #  print("case a")
#                    print("bs1 is", bs1)
#                    print(xld, bm, dm)
                    vsigma=sqrt(vsigmasquared)
                    vpi=sqrt(vpisquared)
                    vrem[0]=- (c+a*vsigmasquared)*0.25/bs1
                    vrem[1]=+ (c+a*vpisquared)*0.25/bs1

                    (vz[0])[0]= vsigma + vrem[0]
                    (vz[0])[1]= vpi  +vrem[1]
                    (vz[0])[2]=-vsigma+vrem[0]
                    (vz[0])[3]=-vpi   +vrem[1]

#                    print("vzs are")
#                    for i in range(4):
#                        print( (vz[0])[i] )

                  #  i=0
                   # print( D21 - (vz[0])[i]*(vz[0])[i]  )
                    D21_min_vzsq[0]=-bm -vrem[0]*(  +2*vsigma + vrem[0] )
                  #  print(D21_min_vzsq[i])

                    #i=1
                   # print( D21 - (vz[0])[i]*(vz[0])[i] )
                    D21_min_vzsq[1] = -bm -xld -vrem[1]*(+2*vpi + vrem[1] )
                   # print(D21_min_vzsq[i])

                   # i=2
                 #   print( D21 - (vz[0])[i]*(vz[0])[i]  )
                    D21_min_vzsq[2]=-bm +vrem[0]*(  +2*vsigma - vrem[0] )
                  #  print(D21_min_vzsq[i])

                   # i=3
                  #  print( D21 - (vz[0])[i]*(vz[0])[i] )
                    D21_min_vzsq[3] = -bm -xld +vrem[1]*(+2*vpi - vrem[1] )
                  #  print(D21_min_vzsq[i])
                    for i in range(4):
                        D34D43_min_vD33vD44[i]= D21_min_vzsq[i] -(vz[0])[i]*a +xld +bm


                else:
                 #   print("case b")
                #    print(cx, cy, cz)
                 #   print(cxy, cyx, cxz, czx, cyz, czy)
                    bs5=-0.5*b + bs1
                    bs6=-0.5*b - bs1
                    vrem[0]=(c+a*bs5)*0.25/bs1
                    vrem[1]=(c+a*bs6)*0.25/bs1
                    bs7=sqrt( bs5 )
                    bs8=sqrt( bs6 )
                    (vz[0])[0]= bs7 - vrem[0]
                    (vz[0])[1]= bs8 + vrem[1]
                    (vz[0])[2]=-bs7 - vrem[0]
                    (vz[0])[3]=-bs8 + vrem[1]

                    D21_min_vzsq[0] =  -( bs4 +vrem[0]*( -2*bs7 + vrem[0] ) )
                    D21_min_vzsq[1] =  - ( bs2 + vrem[1]*(2*bs8+vrem[1])  )
                    D21_min_vzsq[2] =  -( bs4 -vrem[0]*( -2*bs7 - vrem[0] ) )
                    D21_min_vzsq[3] =  - ( bs2 - vrem[1]*(2*bs8-vrem[1])  )
                    for i in range(4):
                        D34D43_min_vD33vD44[i]= D21_min_vzsq[i] -(vz[0])[i]*a +xld +bm
#                        print( D34D43_min_vD33vD44[i] )
#                        print( one_minus_vyvy_div_ezz*(MS[0]).D43 - ( (vz[0])[i]-D33 )*((vz[0])[i]-D44) )

            else:
                Delta0_bd=bquadr+12.*d
                Delta0_ac=-3.*a*c
                Delta0=Delta0_ac+Delta0_bd
                Delta0quadr=Delta0*Delta0

                Delta1_bd=2.*b*(bquadr-36.*d)
                Delta1_ac=-9.*a*b*c+27.*(aquadr*d+c*c)
                Delta1=Delta1_bd+Delta1_ac
                Deltaremainder = Delta0_ac*(3*Delta0_bd*Delta0 +Delta0_ac*Delta0_ac )
                bd_splitter=-4.*(bm*vsigmasquared+dm) + cquadr(bm-xld)
                ri_bd=-432.*d*cquadr(bd_splitter)
                ri_ac=Delta1_ac*(2.*Delta1_bd+Delta1_ac) -4.*Deltaremainder
                ri=ri_bd+ri_ac
                root=sqrt( ri )
                Qinternal=0.5*(Delta1+root)
                Q=pow(Qinternal, 0.3333333333333333)

                qu=0.5*a*xld
                qm=(0.5*a*(-bm +aquadr/4) + vy*(MS[0]).mixmag*(MS[0]).inverseezz)
                q=qu+qm

                pm=-3.*aquadr/8.0
                p=b+pm
                S=0.5*sqrt( -2.*p/3. + (Q+Delta0/Q)/3.0 )
                F1= 9.*b*bd_splitter +0.5*Delta1_ac+0.5*root +8.*pm*(3.*b*p +pm*pm)
                F2=9.*bd_splitter*( bquadr*bquadr - 24.*bquadr*d - 48.*dquadr ) \
                +4.*Delta1_bd*pm*(3*b*p +pm*pm ) + 4.*(Delta1_ac+root)*ccube(p) + Deltaremainder
                Qsquared=Q*Q
                M=(-Q*F1/(3.*( Qinternal-2*p*Qsquared+( Q*cquadr(2*p) ) )  )) - (F2/( ( Delta0quadr-2*p*Q*Delta0+cquadr(2*p*Q) )*3.*Q))



                q_div_S=q/S
                root_MqS1=sqrt(M-q_div_S)
                root_MqS2=sqrt(M+q_div_S)
                vrem[0]= min_a_div_4+0.5*root_MqS1
                vrem[1]= min_a_div_4-0.5*root_MqS1
                vrem[2]= min_a_div_4-0.5*root_MqS2
                vrem[3]= min_a_div_4+0.5*root_MqS2


                (vz[0])[0]= S+vrem[0]
                (vz[0])[1]= S+vrem[1]
                (vz[0])[2]=-S+vrem[2]
                (vz[0])[3]=-S+vrem[3]



                bs1=bm-xld
                bs2=bs1*( bs1*bs1-6*b*vsigmasquared )
                bs3=0.5*(Delta1_ac+root)
                bs4=4*vsigmasquared
                bs5=bs4*Q

                D21_min_SS=- (  -2.*(pm +bs1 ) \
                 + (  bs2 -36*d*bs1 +72*vsigmasquared*(xld*vsigmasquared+dm )     +bs3 )/(Q*Q +Q*bs4 + bs4*bs4 ) \
                 +  (  9.0*bd_splitter*( bquadr*bquadr - 24.0*bquadr*d - 48.0*dquadr ) + Deltaremainder -8*( bs2*Qinternal -b*bquadr*bs3 )  ) /(Q*( Delta0*Delta0 + bs5*Delta0 + cquadr(bs5)  )   )  )/12


                D21_min_vzsq[0] = D21_min_SS  - vrem[0]*(2*S+vrem[0]) -bm
                D21_min_vzsq[1] = D21_min_SS  - vrem[1]*(2*S+vrem[1]) -bm
                D21_min_vzsq[2] = D21_min_SS  - vrem[2]*(-2*S+vrem[2]) -bm
                D21_min_vzsq[3] = D21_min_SS  - vrem[3]*(-2*S+vrem[3]) -bm
                for i in range(4):
                    D34D43_min_vD33vD44[i]= D21_min_vzsq[i] -(vz[0])[i]*a +xld +bm
             #   print("hallo")


                #Force special cases:

                if( ( (D31==0) and ( (MS[0]).D41==0) ) or ( ( (MS[0]).D23==0 )  and (D24==0) ) ):
                    (vz[0])[0]=sqrt(vsigmasquared)
                    (vz[0])[1]=sqrt(vpisquared)
                    (vz[0])[2]=-(vz[0])[0]
                    (vz[0])[3]=-(vz[0])[1]
                    D21_min_vzsq[0]=0
                    D21_min_vzsq[1]=-xld
                    D21_min_vzsq[2]=0
                    D21_min_vzsq[3]=-xld





               # print( D34D43_min_vD33vD44[i] )

            if( ( (vz[0])[0] ).real < 0 ):
               # switch=1
                if( ( (vz[0])[2] ).real >0 ):
                    bs1=(vz[0])[0]
                    (vz[0])[0]=(vz[0])[2]
                    (vz[0])[2]=bs1
                    bs1=D21_min_vzsq[0]
                    D21_min_vzsq[0]=D21_min_vzsq[2]
                    D21_min_vzsq[2]=bs1
                    bs1=D34D43_min_vD33vD44[0]
                    D34D43_min_vD33vD44[0]=D34D43_min_vD33vD44[2]
                    D34D43_min_vD33vD44[2]=bs1
                #    print("switching 0 2")
                else:
                    bs1=(vz[0])[0]
                    (vz[0])[0]=(vz[0])[3]
                    (vz[0])[3]=bs1
                    bs1=D21_min_vzsq[0]
                    D21_min_vzsq[0]=D21_min_vzsq[3]
                    D21_min_vzsq[3]=bs1
                    bs1=D34D43_min_vD33vD44[0]
                    D34D43_min_vD33vD44[0]=D34D43_min_vD33vD44[3]
                    D34D43_min_vD33vD44[3]=bs1
                   # print("switching 0 3")
            if( ( (vz[0])[1] ).real < 0 ):
             #   switch=1
                if( ( (vz[0])[2] ).real >0 ):
                    bs1=(vz[0])[1]
                    (vz[0])[1]=(vz[0])[2]
                    (vz[0])[2]=bs1
                    bs1=D21_min_vzsq[1]
                    D21_min_vzsq[1]=D21_min_vzsq[2]
                    D21_min_vzsq[2]=bs1
                    bs1=D34D43_min_vD33vD44[1]
                    D34D43_min_vD33vD44[1]=D34D43_min_vD33vD44[2]
                    D34D43_min_vD33vD44[2]=bs1
                 #   print("switching 1 2")
                else:
                    bs1=(vz[0])[1]
                    (vz[0])[1]=(vz[0])[3]
                    (vz[0])[3]=bs1
                    bs1=D21_min_vzsq[1]
                    D21_min_vzsq[1]=D21_min_vzsq[3]
                    D21_min_vzsq[3]=bs1
                    bs1=D34D43_min_vD33vD44[1]
                    D34D43_min_vD33vD44[1]=D34D43_min_vD33vD44[3]
                    D34D43_min_vD33vD44[3]=bs1
                  #  print("switching 1 3")

           # i=0

            comp1=Cmaxnorm(D31)
            comp2=Cmaxnorm( (MS[0]).D41 )

       #     print(comp1, comp2, comp3)
       #     print((MS[0]).D41)
            #print("first evs")

            for i in range(0,4,2):
                phipre[0]=1
                phipre[1]=(vz[0])[i]
                comp3=Cmaxnorm( D21_min_vzsq[i] )
                if( (comp1==0) and (comp2==0) and (comp3==0) ):
                    phipre[2] = 0
                    phipre[3] = 0
                else:
                    if( comp1<=comp2 ):
                       # print("hi1")
                        if( (comp1==0) and (comp2==0) ):
                            bs1=0
                        else:
                            bs1=D31/(MS[0]).D41

                        denom=one_minus_vyvy_div_ezz-bs1*(D44- (vz[0])[i] )
                        if(denom==0):
                            denom=cx*2.0e-16
                        rel34u[i] =( (vz[0])[i] - D33 )/denom
                        rel34m[i] =bs1*(MS[0]).D43/denom
                        rel34[i]=rel34m[i]+rel34u[i]
                        bs3 = (MS[0]).D23 + rel34[i]*D24

                        if( (comp2>comp3) or ( (bs3==0) ) ):
                         #   print("use b")
                            bs2=bs1*(D44- (vz[0])[i] )
#                            print(bs1, bs2)
#                            print( (one_minus_vyvy_div_ezz-bs2  )  )
#                            print( (  + (D44-(vz[0])[i] )*rel34m[i] +( D34D43_min_vD33vD44[i]-(MS[0]).D43*bs2 )/(one_minus_vyvy_div_ezz-bs2  ) ) )
                            phipre[2] = -(MS[0]).D41/(  + (D44-(vz[0])[i] )*rel34m[i] +( D34D43_min_vD33vD44[i]-(MS[0]).D43*bs2 )/(one_minus_vyvy_div_ezz-bs2  ) )
                            phipre[3] = rel34[i]*phipre[2]
                        else:
                          #  print("use c")
                            phipre[2] = -D21_min_vzsq[i]/bs3
                            phipre[3] = rel34[i]*phipre[2]

                    elif( comp2<comp1 ):

                        bs1=(MS[0]).D41/D31
                        denom=bs1*one_minus_vyvy_div_ezz-(D44- (vz[0])[i] )
                        if(denom==0):
                            denom=cx*2.0e-16
                        rel34u[i] =(MS[0]).D43/denom
                        rel34m[i] = bs1*( (vz[0])[i] - D33 )/denom
                        rel34[i]=rel34m[i]+rel34u[i]
                        bs3 = (MS[0]).D23 + rel34[i]*D24

                        if( (comp1>comp3) or (bs3==0) ):
                            bs2=bs1*one_minus_vyvy_div_ezz
                        #    print( D33-(vz[0])[i] +one_minus_vyvy_div_ezz*rel34[i] )
                        #    print(  +one_minus_vyvy_div_ezz*rel34m[i] + (  (D33-(vz[0])[i])*bs2+ D34D43_min_vD33vD44[i] )/(bs2-(D44- (vz[0])[i] ))   )
                            phipre[2] = -D31/(  +one_minus_vyvy_div_ezz*rel34m[i] + (  (D33-(vz[0])[i])*bs2+ D34D43_min_vD33vD44[i] )/(bs2-(D44- (vz[0])[i] ))   )
                            phipre[3] = rel34[i]*phipre[2]
                        else:
                            phipre[2] = -D21_min_vzsq[i]/bs3
                            phipre[3] = rel34[i]*phipre[2]


                for j in range(4):
                    (PHI[0])[j][i]=phipre[j]
           # print("second evs")

            for i in range(1,4,2):
                phipre[2]=1
                if( (comp1==0) and (comp2==0) and (comp3==0) ):

                    phipre[0] = 0
                    phipre[1] = 0
                    phipre[3] = ( (vz[0])[i] - D33 )/one_minus_vyvy_div_ezz
                else:
                    if( comp1<=comp2 ):
                        if( (comp1==0) and (comp2==0) ):
                            bs1=0
                        else:
                            bs1=D31/(MS[0]).D41

                        denom=one_minus_vyvy_div_ezz-bs1*(D44- (vz[0])[i] )
                        if(denom==0):
                            denom=cx*2.0e-16
                        rel34u[i] =( (vz[0])[i] - D33 )/denom
                        rel34m[i] =bs1*(MS[0]).D43/denom
                        rel34[i]=rel34m[i]+rel34u[i]
                        phipre[3] = rel34[i]
                        bs3 = (MS[0]).D23 + rel34[i]*D24
                        if( (comp2>comp3) or (bs3==0)  ):
                            bs2=bs1*(D44- (vz[0])[i] )
                            phipre[0] = -(  + (D44-(vz[0])[i] )*rel34m[i] +( D34D43_min_vD33vD44[i]-(MS[0]).D43*bs2 )/(one_minus_vyvy_div_ezz-bs2  ) )/(MS[0]).D41
                            phipre[1] = (vz[0])[i]*phipre[0]
                        else:
                            phipre[0] = -bs3/D21_min_vzsq[i]
                            phipre[1] = (vz[0])[i]*phipre[0]
                    elif( comp2<comp1 ):
                        bs1=(MS[0]).D41/D31
                        denom=bs1*one_minus_vyvy_div_ezz-(D44- (vz[0])[i] )
                        if(denom==0):
                            denom=cx*2.0e-16
                        rel34u[i] =(MS[0]).D43/denom
                        rel34m[i] = bs1*( (vz[0])[i] - D33 )/denom
                        rel34[i]=rel34m[i]+rel34u[i]
                        phipre[3] = rel34[i]
                        bs3 = (MS[0]).D23 + rel34[i]*D24

                        if( (comp1>comp3) or (bs3==0) ):
                            bs2=bs1*one_minus_vyvy_div_ezz
                        #    print( D33-(vz[0])[i] +one_minus_vyvy_div_ezz*rel34[i] )
                        #    print(  +one_minus_vyvy_div_ezz*rel34m[i] + (  (D33-(vz[0])[i])*bs2+ D34D43_min_vD33vD44[i] )/(bs2-(D44- (vz[0])[i] ))   )
                            phipre[0] = -(  +one_minus_vyvy_div_ezz*rel34m[i] + (  (D33-(vz[0])[i])*bs2+ D34D43_min_vD33vD44[i] )/(bs2-(D44- (vz[0])[i] ))   )/D31
                            phipre[1] = (vz[0])[i]*phipre[0]
                        else:
                            phipre[0] = -bs3/D21_min_vzsq[i]
                            phipre[1] = (vz[0])[i]*phipre[0]
                for j in range(4):
                    (PHI[0])[j][i]=phipre[j]

            NormalizePHI(PHI)

            Matrix4Invert(PHI, PSI)





#        i=0
#        while(i<4):
#         #   print("check if Eigenvector", i)
#           # print( (PHI[0])[0][i], (PHI[0])[1][i], (PHI[0])[2][i], (PHI[0])[3][i] )
#            a=(PHI[0])[1][i]
#            b=(vz[0])[i]*(PHI[0])[0][i]
#
#           # print( a,b )
#          #  if(b!=0):
#           #     print( cabsvalue((a-b)/(a+b) ) )
#            if( a==0 and b==0 ):
#                evch[i][0]=0
#            else:
#                #evch[i][3]=cabsvalue((a-b)/(a+b) )
#                evch[i][0]=cabsvalue((a-b) )
#               # print(i, 0, evch[i][0])
#
#            a=(1.+cx-vyvy - czx*cxz/(1+cz))*(PHI[0])[0][i] + (cxy-cxz*czy/(1+cz) )*(PHI[0])[2][i] - (PHI[0])[3][i]*vy*cxz/(1+cz)
#            b=(vz[0])[i]*(PHI[0])[1][i]
#
#
#
#        #    print( a,b )
#        #    if(b!=0):
#        #        print( cabsvalue((a-b)/(a+b) ) )
#            if( a==0 and b==0 ):
#                evch[i][1]=0
#            else:
#                evch[i][3]=cabsvalue((a-b)/(a+b) )
#                evch[i][1]=cabsvalue((a-b) )
#          #      print(i, 1, evch[i][1])
#            a= -( vy*czx/(1.+cz) )*(PHI[0])[0][i] + -( vy*czy/(1.+cz) )*(PHI[0])[2][i] + (1.-vyvy/(1.+cz) )*(PHI[0])[3][i]
#            b=(vz[0])[i]*(PHI[0])[2][i]
#        #    print( a,b )
#        #    if(b!=0):
#         #       print( cabsvalue((a-b)/(a+b) ) )
#            if( a==0 and b==0 ):
#                evch[i][2]=0
#            else:
#                evch[i][3]=cabsvalue((a-b)/(a+b) )
#                evch[i][2]=cabsvalue((a-b) )
#            #    print(i, 2, evch[i][2])
#            a=(cyx-cyz*czx/(1+cz) )*(PHI[0])[0][i] + ( 1.+cy - cyz*czy/(1.+cz) )*(PHI[0])[2][i] +  ( -vy*cyz/(1.+cz) )*(PHI[0])[3][i]
#            b=(vz[0])[i]*(PHI[0])[3][i]
#         #   print( a,b )
#          #  if(b!=0):
#         #       print( cabsvalue((a-b)/(a+b) ) )
#            if( a==0 and b==0 ):
#                evch[i][3]=0
#            else:
#                evch[i][3]=cabsvalue((a-b)/(a+b) )
#                evch[i][3]=cabsvalue((a-b) )
#             #   print(i, 3, evch[i][3])
#
#            i+=1
#         #   print(i)
#        comp1=0
#        for i in range(4):
#            for j in range(4):
#                if( evch[i][j]>comp1):
#                    comp1=evch[i][j]
#        print("im here ", cx)
#        if(comp1>1.0e-13):
#
##                if( evch[i][j]>2.0e-14 ):
#            print("This diagonalization is bad")
#            print( evch[i][j], i, j, vy)
#            print("[", cx, cxy, cxz, cyx, cy, cyz, czx, czy, cz, "]", vy)
#         #####   check if properly inverted:
#
#        print("Inversion:")

#
#        Mult4x4_leftside(PHI, PSI)
#        for i in range(4):
#            for j in range(4):
#                if( i==j ):
#                    (PHI[0])[i][j]-=1
#                if( cabsvalue((PHI[0])[i][j])>1.0e-15 ):
#                    print("This inversion is bad")
#                    print( cabsvalue((PHI[0])[i][j]), i, j, vy )
#                    print("[", cx, cxy, cxz, cyx, cy, cyz, czx, czy, cz, "]", vy)
               ### print( (PHI[0])[i][j] )
    elif( L.type==3 ):
        #Transversal MOKE:
       # print("Moke x")
       # print("MOKE", L.magdir)
        if( L.magdir==1 ):
          #  print("hallo")
            (vz[0])[0]=CalculateVZsigma(vyvy, cx)
            (vz[0])[1]=CalculateVZpi_m(vyvy, cy, cz, cg)
            (vz[0])[2]=-(vz[0])[0]
            (vz[0])[3]=-(vz[0])[1]
            a = 1./(1.+cz)
            b=( (vz[0])[1]-vy*cg*a )
            c=1.0/(1.-vyvy*a )
            d=( (vz[0])[3]-vy*cg*a )
            S=2*(vz[0])[3]
            (PHI[0])[0][0]=1
            (PHI[0])[1][0]=(vz[0])[0]
            (PHI[0])[2][0]=0
            (PHI[0])[3][0]=0

            (PHI[0])[0][1]=0
            (PHI[0])[1][1]=0
            (PHI[0])[2][1]=1
            (PHI[0])[3][1]=b*c

            (PHI[0])[0][2]=1
            (PHI[0])[1][2]=(vz[0])[2]
            (PHI[0])[2][2]=0
            (PHI[0])[3][2]=0

            (PHI[0])[0][3]=0
            (PHI[0])[1][3]=0
            (PHI[0])[2][3]=1
            (PHI[0])[3][3]=d*c

            (PSI[0])[0][0]=0.5
            (PSI[0])[0][1]=0.5/(vz[0])[0]
            (PSI[0])[0][2]=0
            (PSI[0])[0][3]=0

            (PSI[0])[1][0]=0
            (PSI[0])[1][1]=0
            (PSI[0])[1][2]=d/S
            (PSI[0])[1][3]=-1/(  c*S  )

            (PSI[0])[2][0]=0.5
            (PSI[0])[2][1]=-(PSI[0])[0][1]
            (PSI[0])[2][2]=0
            (PSI[0])[2][3]=0

            (PSI[0])[3][0]=0
            (PSI[0])[3][1]=0
            (PSI[0])[3][2]=-b/S
            (PSI[0])[3][3]=-(PSI[0])[1][3]

#            i=0
#            print("check if Eigenvector", i)
#            print( (PHI[0])[1][i], (vz[0])[i]*(PHI[0])[0][i]  )
#            print( (1.+cx-vyvy)*(PHI[0])[0][i], (vz[0])[i]*(PHI[0])[1][i]  )
#            print(  ( vy*cg/(1.+cz) )*(PHI[0])[2][i] + (1.-vyvy/(1.+cz) )*(PHI[0])[3][i]  , (vz[0])[i]*(PHI[0])[2][i]  )
#            print( ( 1.+cy +cg*cg/(1.+cz) )*(PHI[0])[2][i] +  ( -vy*cg/(1.+cz) )*(PHI[0])[3][i], (vz[0])[i]*(PHI[0])[3][i] )
#            i=1
#            print("check if Eigenvector", i)
#            print( (PHI[0])[1][i], (vz[0])[i]*(PHI[0])[0][i]  )
#            print( (1.+cx-vyvy)*(PHI[0])[0][i], (vz[0])[i]*(PHI[0])[1][i]  )
#            print(  ( vy*cg/(1.+cz) )*(PHI[0])[2][i] + (1.-vyvy/(1.+cz) )*(PHI[0])[3][i]  , (vz[0])[i]*(PHI[0])[2][i]  )
#            print( ( 1.+cy +cg*cg/(1.+cz) )*(PHI[0])[2][i] +  ( -vy*cg/(1.+cz) )*(PHI[0])[3][i], (vz[0])[i]*(PHI[0])[3][i] )
#
#            i=0
#            print("check if Eigenvector", i)
#            print( (PHI[0])[1][i], (vz[0])[i]*(PHI[0])[0][i]  )
#            print( (1.+cx-vyvy)*(PHI[0])[0][i], (vz[0])[i]*(PHI[0])[1][i]  )
#            print(  ( vy*cg/(1.+cz) )*(PHI[0])[2][i] + (1.-vyvy/(1.+cz) )*(PHI[0])[3][i]  , (vz[0])[i]*(PHI[0])[2][i]  )
#            print( ( 1.+cy +cg*cg/(1.+cz) )*(PHI[0])[2][i] +  ( -vy*cg/(1.+cz) )*(PHI[0])[3][i], (vz[0])[i]*(PHI[0])[3][i] )
#            i=1
#            print("check if Eigenvector", i)
#            print( (PHI[0])[1][i], (vz[0])[i]*(PHI[0])[0][i]  )
#            print( (1.+cx-vyvy)*(PHI[0])[0][i], (vz[0])[i]*(PHI[0])[1][i]  )
#            print(  ( vy*cg/(1.+cz) )*(PHI[0])[2][i] + (1.-vyvy/(1.+cz) )*(PHI[0])[3][i]  , (vz[0])[i]*(PHI[0])[2][i]  )
#            print( ( 1.+cy +cg*cg/(1.+cz) )*(PHI[0])[2][i] +  ( -vy*cg/(1.+cz) )*(PHI[0])[3][i], (vz[0])[i]*(PHI[0])[3][i] )

#            ##check if properly inverted:
#            Mult4x4_leftside(PHI, PSI)
#            for i in range(4):
#                for j in range(4):
#                    print( (PHI[0])[i][j] )
            #####

        if( L.magdir==2 ):
         #   print("Moke y")


            a=1.0/(1.+cz)
            xld= (cx-cy) + vyvy*(cy-cz)*a
            c = cquadr( xld )
            du=cx+cy+cx*cy
            dm=1.0+a
            p=-vy*cg*a
            bm=-cg*cg*a
            bu=-cy-cx+vyvy*cy*a
            Delta0=bm+bu
            Delta1=-2.0+vyvy*dm
            Q=bm*(2.*bu + bm + 2.0*(  vyvy*dm + 2*cy   ) )
            root=sqrt( Q + c )
            b=Delta0+Delta1
            (vz[0])[0]=sqrt((-b-root)/2)
            (vz[0])[1]=sqrt((-b+root)/2)
            (vz[0])[2]=-(vz[0])[0]
            (vz[0])[3]=-(vz[0])[1]
            bs1 = xld-root
            bs2 = xld+root

            if( Cmaxnorm(bs1) < Cmaxnorm( bs2 ) ):
             #   print("case a")
                phipre[3]=(  bm -bs2  )/(2*p)
            else:
            #    print("case b")
                phipre[3]=( bm +Q/bs1 )/(2*p)
            phipre[2]=phipre[3]*(vz[0])[0]/(1.+cy)
            phipre[0]=1
            phipre[1]=(vz[0])[0]

            for j in range(4):
                (PHI[0])[j][0] = phipre[j]

            phipre[2]=1
            phipre[3]=(1+cy)/(vz[0])[1]
            phipre[1]=(PHI[0])[3][0]
            phipre[0]=phipre[1]/(vz[0])[1]
            #check if Eigenvector:
            for j in range(4):
                (PHI[0])[j][1] = phipre[j]

            (PHI[0])[0][2]=(PHI[0])[0][0]
            (PHI[0])[1][2]=-(PHI[0])[1][0]
            (PHI[0])[2][2]=-(PHI[0])[2][0]
            (PHI[0])[3][2]=(PHI[0])[3][0]

            (PHI[0])[0][3]=-(PHI[0])[0][1]
            (PHI[0])[1][3]=(PHI[0])[1][1]
            (PHI[0])[2][3]=(PHI[0])[2][1]
            (PHI[0])[3][3]=-(PHI[0])[3][1]
#            i=3
#            print("check if Eigenvector", i)
#            print( (PHI[0])[1][i], (vz[0])[i]*(PHI[0])[0][i]  )
#            print( (1.+cx-vyvy + cg*cg*a )*(PHI[0])[0][i] - ( vy*cg/(1.+cz) )*(PHI[0])[3][i] , (vz[0])[i]*(PHI[0])[1][i]  )
#            print(  ( vy*cg/(1.+cz) )*(PHI[0])[0][i] + (1.-vyvy/(1.+cz) )*(PHI[0])[3][i]  , (vz[0])[i]*(PHI[0])[2][i]  )
#            print( ( 1.+cy  )*(PHI[0])[2][i] , (vz[0])[i]*(PHI[0])[3][i] )

          #  (PHI[0])
            a=2*( (PHI[0])[3][0]*(PHI[0])[0][1] - (PHI[0])[3][1] )
            b=2*( (PHI[0])[1][0] - (PHI[0])[2][0]*(PHI[0])[1][1] )

            (PSI[0])[0][0] = - (PHI[0])[3][1]/a
            (PSI[0])[1][0] =  (PHI[0])[3][0]/a
            (PSI[0])[2][0] = (PSI[0])[0][0]
            (PSI[0])[3][0] = - (PSI[0])[1][0]

            (PSI[0])[0][1] = 1.0/b
            (PSI[0])[1][1] = -(PHI[0])[2][0]/b
            (PSI[0])[2][1] = -(PSI[0])[0][1]
            (PSI[0])[3][1] = (PSI[0])[1][1]

            (PSI[0])[0][2] = -(PHI[0])[1][1]/b
            (PSI[0])[1][2] = (PHI[0])[1][0]/b
            (PSI[0])[2][2] = -(PSI[0])[0][2]
            (PSI[0])[3][2] = (PSI[0])[1][2]

            (PSI[0])[0][3] = (PHI[0])[0][1]/a
            (PSI[0])[1][3] = - 1./a
            (PSI[0])[2][3] = (PSI[0])[0][3]
            (PSI[0])[3][3] = - (PSI[0])[1][3]

            ##check if properly inverted:
#            Mult4x4_leftside(PHI, PSI)
#            for i in range(4):
#                for j in range(4):
#                    print( (PHI[0])[i][j] )
        if(L.magdir == 3):


            a=1.0/(1.+cz)
            xld= (cx-cy) + vyvy*(cy-cz)*a
            c = cquadr(xld)

            du=(1-vyvy*a)
            dm=( cy+1.0 )*a

            bu=(  dm  +1)*(1.-vyvy)/2
            bm=(cz*dm+cx)/2
            root=sqrt(0.25*c-cg*cg*du )

            b=bu + bm
            (vz[0])[0]=sqrt(b+root)
            (vz[0])[1]=sqrt(b-root)
            (vz[0])[2]=-(vz[0])[0]
            (vz[0])[3]=-(vz[0])[1]

            bs1 = xld/2-root
            bs2 = xld/2+root

            phipre[0]=1
            phipre[1]=(vz[0])[0]
            if( Cmaxnorm(bs1) < Cmaxnorm( bs2 ) ):
             #   print("a")
                phipre[2]=  -cg*du/( bs2 )
            else:
              #  print("b")
                phipre[2]=(  -bs1 )/cg
           # phipre[2]=(  -bs1 )/cg
            phipre[3]=((vz[0])[0]/du)*phipre[2]
            for j in range(4):
                (PHI[0])[j][0] = phipre[j]



            (PHI[0])[2][1] = 1
            (PHI[0])[3][1] = (vz[0])[1]/du
            (PHI[0])[0][1] = (PHI[0])[2][0]/du
            (PHI[0])[1][1] = (vz[0])[1]*(PHI[0])[0][1]



            (PHI[0])[0][2]=(PHI[0])[0][0]
            (PHI[0])[1][2]=-(PHI[0])[1][0]
            (PHI[0])[2][2]=(PHI[0])[2][0]
            (PHI[0])[3][2]=-(PHI[0])[3][0]

            (PHI[0])[0][3]=(PHI[0])[0][1]
            (PHI[0])[1][3]=-(PHI[0])[1][1]
            (PHI[0])[2][3]=(PHI[0])[2][1]
            (PHI[0])[3][3]=-(PHI[0])[3][1]
#            i=0
#            print((PHI[0])[0][i], (PHI[0])[1][i], (PHI[0])[2][i], (PHI[0])[3][i])
#            print("check if Eigenvector", i)
#            print( (PHI[0])[1][i], (vz[0])[i]*(PHI[0])[0][i]  )
#            print( (1.+cx-vyvy )*(PHI[0])[0][i] +cg*(PHI[0])[2][i] , (vz[0])[i]*(PHI[0])[1][i]  )
#            print(   (1.-vyvy/(1.+cz) )*(PHI[0])[3][i]  , (vz[0])[i]*(PHI[0])[2][i]  )
#            print( -cg*(PHI[0])[0][i] + ( 1.+cy  )*(PHI[0])[2][i] , (vz[0])[i]*(PHI[0])[3][i] )
#            i=1
#            print((PHI[0])[0][i], (PHI[0])[1][i], (PHI[0])[2][i], (PHI[0])[3][i])
#            print("check if Eigenvector", i)
#            print( (PHI[0])[1][i], (vz[0])[i]*(PHI[0])[0][i]  )
#            print( (1.+cx-vyvy )*(PHI[0])[0][i] +cg*(PHI[0])[2][i] , (vz[0])[i]*(PHI[0])[1][i]  )
#            print(   (1.-vyvy/(1.+cz) )*(PHI[0])[3][i]  , (vz[0])[i]*(PHI[0])[2][i]  )
#            print( -cg*(PHI[0])[0][i] + ( 1.+cy  )*(PHI[0])[2][i] , (vz[0])[i]*(PHI[0])[3][i] )
#            i=2
#            print((PHI[0])[0][i], (PHI[0])[1][i], (PHI[0])[2][i], (PHI[0])[3][i])
#            print("check if Eigenvector", i)
#            print( (PHI[0])[1][i], (vz[0])[i]*(PHI[0])[0][i]  )
#            print( (1.+cx-vyvy )*(PHI[0])[0][i] +cg*(PHI[0])[2][i] , (vz[0])[i]*(PHI[0])[1][i]  )
#            print(   (1.-vyvy/(1.+cz) )*(PHI[0])[3][i]  , (vz[0])[i]*(PHI[0])[2][i]  )
#            print( -cg*(PHI[0])[0][i] + ( 1.+cy  )*(PHI[0])[2][i] , (vz[0])[i]*(PHI[0])[3][i] )
#            i=3
#            print((PHI[0])[0][i], (PHI[0])[1][i], (PHI[0])[2][i], (PHI[0])[3][i])
#            print("check if Eigenvector", i)
#            print( (PHI[0])[1][i], (vz[0])[i]*(PHI[0])[0][i]  )
#            print( (1.+cx-vyvy )*(PHI[0])[0][i] +cg*(PHI[0])[2][i] , (vz[0])[i]*(PHI[0])[1][i]  )
#            print(   (1.-vyvy/(1.+cz) )*(PHI[0])[3][i]  , (vz[0])[i]*(PHI[0])[2][i]  )
#            print( -cg*(PHI[0])[0][i] + ( 1.+cy  )*(PHI[0])[2][i] , (vz[0])[i]*(PHI[0])[3][i] )

            ###Norm:




            a=2*( (PHI[0])[2][0]*(PHI[0])[0][1] - 1 )
            b=2*( (PHI[0])[1][0]*(PHI[0])[3][1] - (PHI[0])[3][0]*(PHI[0])[1][1] )

            (PSI[0])[0][0] = -1./a
            (PSI[0])[1][0] = (PHI[0])[2][0]/a
            (PSI[0])[2][0] = (PSI[0])[0][0]
            (PSI[0])[3][0] = (PSI[0])[1][0]

            (PSI[0])[0][1] = (PHI[0])[3][1]/b
            (PSI[0])[1][1] = -(PHI[0])[3][0] /b
            (PSI[0])[2][1] = -(PSI[0])[0][1]
            (PSI[0])[3][1] = -(PSI[0])[1][1]

            (PSI[0])[0][2] = (PHI[0])[0][1]/a
            (PSI[0])[1][2] = -1 /a
            (PSI[0])[2][2] = (PSI[0])[0][2]
            (PSI[0])[3][2] = (PSI[0])[1][2]

            (PSI[0])[0][3] = -(PHI[0])[1][1]/b
            (PSI[0])[1][3] = (PHI[0])[1][0]/b
            (PSI[0])[2][3] = -(PSI[0])[0][3]
            (PSI[0])[3][3] = -(PSI[0])[1][3]
            #check if properly inverted:
#            NormalizePHI(PHI)
#            a=2*( (PHI[0])[2][0]*(PHI[0])[0][1] - (PHI[0])[2][1]*(PHI[0])[0][0] )
#            b=2*( (PHI[0])[1][0]*(PHI[0])[3][1] - (PHI[0])[3][0]*(PHI[0])[1][1] )
#
#            (PSI[0])[0][0] = -(PHI[0])[2][1]/a
#            (PSI[0])[1][0] = (PHI[0])[2][0]/a
#            (PSI[0])[2][0] = (PSI[0])[0][0]
#            (PSI[0])[3][0] = (PSI[0])[1][0]
#
#            (PSI[0])[0][1] = (PHI[0])[3][1]/b
#            (PSI[0])[1][1] = -(PHI[0])[3][0] /b
#            (PSI[0])[2][1] = -(PSI[0])[0][1]
#            (PSI[0])[3][1] = -(PSI[0])[1][1]
#
#            (PSI[0])[0][2] = (PHI[0])[0][1]/a
#            (PSI[0])[1][2] = -(PHI[0])[0][0]/a
#            (PSI[0])[2][2] = (PSI[0])[0][2]
#            (PSI[0])[3][2] = (PSI[0])[1][2]
#
#            (PSI[0])[0][3] = -(PHI[0])[1][1]/b
#            (PSI[0])[1][3] = (PHI[0])[1][0]/b
#            (PSI[0])[2][3] = -(PSI[0])[0][3]
#            (PSI[0])[3][3] = -(PSI[0])[1][3]
            #check if properly inverted:




 #           Mult4x4_leftside(PHI, PSI)
#            for i in range(4):
#                for j in range(4):
#                    print( (PHI[0])[i][j] )
    else:
        (vz[0])[0]=CalculateVZsigma(vyvy, cx)
        (vz[0])[1]=CalculateVZpi(vyvy, cy, cz)
        (vz[0])[2]=-(vz[0])[0]
        (vz[0])[3]=-(vz[0])[1]
        a=1.+cy
        (PHI[0])[0][0]=1
        (PHI[0])[1][0]=(vz[0])[0]
        (PHI[0])[2][0]=0
        (PHI[0])[3][0]=0

        (PHI[0])[0][1]=0
        (PHI[0])[1][1]=0
        (PHI[0])[2][1]=1
        (PHI[0])[3][1]=a/(vz[0])[1]

        (PHI[0])[0][2]=1
        (PHI[0])[1][2]=(vz[0])[2]
        (PHI[0])[2][2]=0
        (PHI[0])[3][2]=0

        (PHI[0])[0][3]=0
        (PHI[0])[1][3]=0
        (PHI[0])[2][3]=1
        (PHI[0])[3][3]=-(PHI[0])[3][1]

        (PSI[0])[0][0]=0.5
        (PSI[0])[0][1]=0.5/(vz[0])[0]
        (PSI[0])[0][2]=0
        (PSI[0])[0][3]=0

        (PSI[0])[1][0]=0
        (PSI[0])[1][1]=0
        (PSI[0])[1][2]=0.5
        (PSI[0])[1][3]=0.5*(vz[0])[1]/a

        (PSI[0])[2][0]=0.5
        (PSI[0])[2][1]=0.5/(vz[0])[2]
        (PSI[0])[2][2]=0
        (PSI[0])[2][3]=0

        (PSI[0])[3][0]=0
        (PSI[0])[3][1]=0
        (PSI[0])[3][2]=0.5
        (PSI[0])[3][3]=-(PSI[0])[1][3]




cdef void Full_Matrix(Heterostructure* HS, MatrixSafer *AllMS, int* Layer_type_to_Matrixsafe, double th, double wavelength, double complex (*rtot)[2][2]):
    cdef double k0=6.283185307179586/wavelength
    cdef double complex ik0=1j*k0

    cdef double vy=cos(two_pi_div_360()*th)




    cdef double vyvy=vy*vy

    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j,k,l

    cdef CLayer UpperLayer, LowerLayer

    cdef double complex PSI1[4][4]
    cdef double complex PSI2[4][4]
    cdef double complex PHI1[4][4]
    cdef double complex PHI2[4][4]

    cdef double complex (*PHIpointer_upper)[4][4]
    cdef double complex (*PSIpointer_upper)[4][4]
    cdef double complex (*PHIpointer_lower)[4][4]
    cdef double complex (*PSIpointer_lower)[4][4]


    cdef double complex W[4][4]
    cdef double complex P[4][4]
    cdef double complex Pcomp[4][4]
    cdef double complex PHIcomp_start[4][4]
    cdef double complex PSIcomp_start[4][4]




    cdef double complex vzstart[4]
    cdef double complex pstart[4]
    cdef double complex p[4]

    cdef double complex vz1[4]
    cdef double complex vz2[4]
    cdef int vzfilled
    cdef double complex (*vzpointer_lower)[4]
    cdef double complex (*vzpointer_upper)[4]
    cdef double vzvac = sin(two_pi_div_360()*th)
    cdef double vzvz=vzvac*vzvac
    cdef double complex PHIvac[4][4]

    cdef int Cap=NLAYERS-1
    cdef int Upper, Lower
    cdef double complex divide
    PHIvac[0][0]=1
    PHIvac[1][0]=vzvac
    PHIvac[2][0]=0
    PHIvac[3][0]=0

    PHIvac[0][1]=0
    PHIvac[1][1]=0
    PHIvac[2][1]=vzvac
    PHIvac[3][1]=1

    PHIvac[0][2]=1
    PHIvac[1][2]=-vzvac
    PHIvac[2][2]=0
    PHIvac[3][2]=0

    PHIvac[0][3]=0
    PHIvac[1][3]=0
    PHIvac[2][3]=-vzvac
    PHIvac[3][3]=1
    Lower=MLCOMP[0][0]
    LowerLayer=LR[Lower]


    Calculate_Phi_and_Psi(LowerLayer, &( AllMS[ Layer_type_to_Matrixsafe[Lower] ] ), vy, vzvz, vyvy, &vz1, &PHI1, &P)


#    print("initial PHI1")
#    for i in range(4):
#        for j in range(4):
#            print(i,j,P[i][j])

    if(NLAYERS==1):
        Mult4x4_leftside(&P, &PHIvac)
        vz2[0]=vzvac
        vz2[1]=vzvac
        vz2[2]=-vzvac
        vz2[3]=-vzvac
    else:
        Upper=MLCOMP[1][0]
        UpperLayer=LR[Upper]
        Calculate_Phi_and_Psi(UpperLayer, &( AllMS[ Layer_type_to_Matrixsafe[Upper] ] ), vy, vzvz, vyvy, &vz2, &PHI2, &PSI2)



        Mult4x4_leftside(&P, &PHI2)
        vzfilled=2
        for k in range(4):
            p[k] = exp(-ik0*UpperLayer.Thickness*vz2[k] )

#    print("vy ", vy)
#    print("Layer 1")
#    print("PSI-matrix")
#    for k in range(4):
#        for l in range(4):
#            print( PHI2[k][l] )
#    print("vz")
#    for k in range(4):
#        print( vz2[k] )

    ######################################################

    roughfac =-0.5*quadr( LowerLayer.Roughness*k0)

    for k in range(4):
        for l in range(4):

            W[k][l] =exp(roughfac*cquadr(vz1[k] - vz2[l] ) )
            P[k][l]*=W[k][l]
#    roughfac =-0.5*quadr( LowerLayer.Roughness*k0)
#
#    for k in range(2):
#        for l in range(2,4):
#
#            W[k][l] =exp(roughfac*4*vzvac*vzvac*(1. + 1.0e-15j) )
#            P[k][l]*=W[k][l]
#    for k in range(2,4):
#        for l in range(2):
#
#            W[k][l] =exp(roughfac*4*vzvac*vzvac*(1. + 1.0e-15j) )
#            P[k][l]*=W[k][l]

    ######################################################


    i=1
    while i<NLAYERS:
        if(MLLENGTH[i]==1):
            Lower=MLCOMP[i][0]
            LowerLayer=LR[ Lower ]

          #  print("lower is", Lower)

            if(vzfilled==2):
                vzpointer_lower=&vz2
                vzpointer_upper=&vz1
                PHIpointer_lower=&PHI2
                PHIpointer_upper=&PHI1
                PSIpointer_lower=&PSI2
                PSIpointer_upper=&PSI1
                vzfilled=1

            else:
                vzpointer_lower=&vz1
                vzpointer_upper=&vz2
                PHIpointer_lower=&PHI1
                PHIpointer_upper=&PHI2
                PSIpointer_lower=&PSI1
                PSIpointer_upper=&PSI2
                vzfilled=2


            Mult4x4_leftside_diag(&P, &p)
            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
                UpperLayer=LR[ Upper ]
                Calculate_Phi_and_Psi(UpperLayer, &( AllMS[ Layer_type_to_Matrixsafe[Upper] ] ), vy, vzvz, vyvy, vzpointer_upper, PHIpointer_upper, PSIpointer_upper)





                Mult4x4_leftside(PSIpointer_lower, PHIpointer_upper)

            else:

                Mult4x4_leftside(PSIpointer_lower, &PHIvac)
                (vzpointer_upper[0])[0]=vzvac
                (vzpointer_upper[0])[1]=vzvac
                (vzpointer_upper[0])[2]=-vzvac
                (vzpointer_upper[0])[3]=-vzvac
#                print("Cap vz")
#                for k in range(4):
#                    print( (vzpointer_lower[0])[k] )


            ###############################
            roughfac =-0.5*quadr( LowerLayer.Roughness*k0)

            for k in range(4):
                for l in range(4):
                    W[k][l] =exp(roughfac*cquadr( (vzpointer_lower[0])[k] - (vzpointer_upper[0])[l] ) )
                    (PSIpointer_lower[0])[k][l]*=W[k][l]
#            roughfac =-0.5*quadr( LowerLayer.Roughness*k0)
#            for k in range(4):
#                for l in range(4):
#                    W[k][l] =exp(roughfac*cquadr(vz1[k] - vz2[l] ) )
#                    (PSIpointer_lower[0])[k][l]*=W[k][l]
#            roughfac =-0.5*quadr( LowerLayer.Roughness*k0)
#
#            for k in range(2):
#                for l in range(2,4):
#
#                    W[k][l] =exp(roughfac*4*vzvac*vzvac*(1. + 1.0e-15j) )
#                    P[k][l]*=W[k][l]
#            for k in range(2,4):
#                for l in range(2):
#
#                    W[k][l] =exp(roughfac*4*vzvac*vzvac*(1. + 1.0e-15j) )
#                    P[k][l]*=W[k][l]
            ###############################

            Mult4x4_leftside(&P, PSIpointer_lower)




            for k in range(4):
                p[k] = exp(-ik0*UpperLayer.Thickness*(vzpointer_upper[0])[k] )
        else:
            ##Copy data from the recent layer


            if(vzfilled==1):
#                ###
#                NormalizePHI(&PHI1)
#                NormalizePSI(&PSI1)
#                ###
                for k in range(4):
                    for l in range(4):
                        PHIcomp_start[k][l] = PHI1[k][l]
                        PSIcomp_start[k][l] = PSI1[k][l]
                        Pcomp[k][l] = PSI1[k][l]
                    vzstart[k]=vz1[k]
                    pstart[k]=p[k]
            else:
#                ###
#                NormalizePHI(&PHI2)
#                NormalizePSI(&PSI2)
#                ###
                for k in range(4):
                    for l in range(4):
                        PHIcomp_start[k][l] = PHI2[k][l]
                        PSIcomp_start[k][l] = PSI2[k][l]
                        Pcomp[k][l] = PSI2[k][l]
                    vzstart[k]=vz2[k]
                    pstart[k]=p[k]

            Lower=MLCOMP[i][0]
            LowerLayer=LR[ Lower ]
            if(vzfilled==2):
                vzpointer_lower=&vz2
                vzpointer_upper=&vz1
                PHIpointer_lower=&PHI2
                PHIpointer_upper=&PHI1
                PSIpointer_lower=&PSI2
                PSIpointer_upper=&PSI1
                vzfilled=1

            else:
                vzpointer_lower=&vz1
                vzpointer_upper=&vz2
                PHIpointer_lower=&PHI1
                PHIpointer_upper=&PHI2
                PSIpointer_lower=&PSI1
                PSIpointer_upper=&PSI2
                vzfilled=2


            Upper=MLCOMP[i][1]
            UpperLayer=LR[Upper]
            Calculate_Phi_and_Psi(UpperLayer, &( AllMS[ Layer_type_to_Matrixsafe[Upper] ] ), vy, vzvz, vyvy, vzpointer_upper, PHIpointer_upper, PSIpointer_upper)
#            ###
#            NormalizePHI(PHIpointer_upper)
#            NormalizePSI(PSIpointer_upper)
#            ###

            for k in range(4):
                p[k] = exp(-ik0*UpperLayer.Thickness*(vzpointer_upper[0])[k] )

            Mult4x4_leftside(&Pcomp, PHIpointer_upper)


#            roughfac =-0.5*quadr( LowerLayer.Roughness*k0)
#            for k in range(4):
#                for l in range(4):
#                    W[k][l] =exp(roughfac*cquadr(vz1[k] - vz2[l] ) )
#                    Pcomp[k][l]*=W[k][l]
            roughfac =-0.5*quadr( LowerLayer.Roughness*k0)
            for k in range(4):
                for l in range(4):
                    W[k][l] =exp(roughfac*cquadr( vzstart[k] - (vzpointer_upper[0])[l] ) )
                    Pcomp[k][l]*=W[k][l]

            j=1
            while j<( MLLENGTH[i]-1 ):
                Lower=MLCOMP[i][j]
                Upper=MLCOMP[i][(j+1)]
                j+=1
                UpperLayer=LR[ Upper ]
                LowerLayer=LR[ Lower ]
                if(vzfilled==2):
                    vzpointer_lower=&vz2
                    vzpointer_upper=&vz1
                    PHIpointer_lower=&PHI2
                    PHIpointer_upper=&PHI1
                    PSIpointer_lower=&PSI2
                    PSIpointer_upper=&PSI1
                    vzfilled=1

                else:
                    vzpointer_lower=&vz1
                    vzpointer_upper=&vz2
                    PHIpointer_lower=&PHI1
                    PHIpointer_upper=&PHI2
                    PSIpointer_lower=&PSI1
                    PSIpointer_upper=&PSI2
                    vzfilled=2

                Mult4x4_leftside_diag(&Pcomp, &p)
                Calculate_Phi_and_Psi(UpperLayer, &( AllMS[ Layer_type_to_Matrixsafe[Upper] ] ), vy, vzvz, vyvy, vzpointer_upper, PHIpointer_upper, PSIpointer_upper)

#                ###
#                NormalizePHI(PHIpointer_upper)
#                NormalizePSI(PSIpointer_upper)
#                ###

                for k in range(4):
                    p[k] = exp(-ik0*UpperLayer.Thickness*(vzpointer_upper[0])[k] )

                Mult4x4_leftside(PSIpointer_lower, PHIpointer_upper)
                roughfac =-0.5*quadr( LowerLayer.Roughness*k0)
                for k in range(4):
                    for l in range(4):
                        W[k][l] =exp(roughfac*cquadr( (vzpointer_lower[0])[k] - (vzpointer_upper[0])[l] ) )
                        (PSIpointer_lower[0])[k][l]*=W[k][l]
#                for k in range(4):
#                    for l in range(4):
#                        W[k][l] =exp(roughfac*cquadr(vz1[k] - vz2[l] ) )
#                        (PSIpointer_lower[0])[k][l]*=W[k][l]

                Mult4x4_leftside(&Pcomp, PSIpointer_lower)

            Mult4x4_leftside_diag(&Pcomp, &p)

            Mult4x4_leftside(PSIpointer_upper, &PHIcomp_start)
            roughfac =-0.5*quadr( UpperLayer.Roughness*k0)
            for k in range(4):
                for l in range(4):
                    W[k][l] =exp(roughfac*cquadr( (vzpointer_upper[0])[k] - vzstart[l] ) )
                    (PSIpointer_upper[0])[k][l]*=W[k][l]

            Mult4x4_leftside(&Pcomp, PSIpointer_upper)

            Mult4x4_leftside_diag(&Pcomp, &pstart)

         #   NormalizePSI(&Pcomp)
#            for k in range(4):
#                Pcomp[k][0]=Pcomp[k][0]*(1.1)
#                Pcomp[k][1]=Pcomp[k][1]*(1.1)
#            for k in range(4):
#                Pcomp[k][2]=Pcomp[k][2]/(1.1)
#                Pcomp[k][3]=Pcomp[k][3]/(1.1)
            Matrixexp( &Pcomp, MLREP[i]-1 )



#            print("this time Pcomp is")
#            for k in range(4):
#                for l in range(4):
#                    print(Pcomp[k][l])

            Mult4x4_leftside_diag(&P, &pstart)

            Mult4x4_leftside(&P, &Pcomp)

            Lower=MLCOMP[i][0]
            Upper=MLCOMP[i][1]
            UpperLayer=LR[ Upper ]
            LowerLayer=LR[ Lower ]



#            Calculate_Phi_and_Psi(UpperLayer, &( AllMS[ Layer_type_to_Matrixsafe[Upper] ] ), vy, vzvz, vyvy, vzpointer_upper, PHIpointer_upper, PSIpointer_upper)
#            Mult4x4_leftside(&PSIcomp_start, PHIpointer_upper)
#            roughfac =-0.5*quadr( LowerLayer.Roughness*k0)
#            for k in range(4):
#                for l in range(4):
#                    W[k][l] =exp(roughfac*cquadr( (vzpointer_upper[0])[k] - vzstart[l] ) )
#                    PSIcomp_start[k][l]*=W[k][l]
#
#            Mult4x4_leftside(&P,&PSIcomp_start)


            Calculate_Phi_and_Psi(UpperLayer, &( AllMS[ Layer_type_to_Matrixsafe[Upper] ] ), vy, vzvz, vyvy, &vz2, &PHI2, &PSI2)
            Mult4x4_leftside(&PSIcomp_start, &PHI2)
            roughfac =-0.5*quadr( LowerLayer.Roughness*k0)
#            for k in range(4):
#                for l in range(4):
#                    W[k][l] =exp(roughfac*cquadr( vz2[k] - vzstart[l] ) )
#                    PSIcomp_start[k][l]*=W[k][l]

            for k in range(4):
                for l in range(4):
                    W[k][l] =exp(roughfac*cquadr( vzstart[k] - vz2[l] ) )
                    PSIcomp_start[k][l]*=W[k][l]


#            print("hallo1")
#            for k in range(4):
#                for l in range(4):
#                    print( (PSIpointer_lower[0])[k][l]  )
            Mult4x4_leftside(&P,&PSIcomp_start)
            vzfilled=2


            for k in range(4):
                p[k] = exp(-ik0*UpperLayer.Thickness*vz2[k] )

            j=1
            while j<( MLLENGTH[i]-1 ):
                Lower=MLCOMP[i][j]
                Upper=MLCOMP[i][(j+1)]
                j+=1
                UpperLayer=LR[ Upper ]
                LowerLayer=LR[ Lower ]
                if(vzfilled==2):
                    vzpointer_lower=&vz2
                    vzpointer_upper=&vz1
                    PHIpointer_lower=&PHI2
                    PHIpointer_upper=&PHI1
                    PSIpointer_lower=&PSI2
                    PSIpointer_upper=&PSI1
                    vzfilled=1

                else:
                    vzpointer_lower=&vz1
                    vzpointer_upper=&vz2
                    PHIpointer_lower=&PHI1
                    PHIpointer_upper=&PHI2
                    PSIpointer_lower=&PSI1
                    PSIpointer_upper=&PSI2
                    vzfilled=2
                Mult4x4_leftside_diag(&P, &p)
                Calculate_Phi_and_Psi(UpperLayer, &( AllMS[ Layer_type_to_Matrixsafe[Upper] ] ), vy, vzvz, vyvy, vzpointer_upper, PHIpointer_upper, PSIpointer_upper)
                for k in range(4):
                    p[k] = exp(-ik0*UpperLayer.Thickness*(vzpointer_upper[0])[k] )

                Mult4x4_leftside(PSIpointer_lower, PHIpointer_upper)
                roughfac =-0.5*quadr( LowerLayer.Roughness*k0)
#                for k in range(4):
#                    for l in range(4):
#                        W[k][l] =exp(roughfac*cquadr(vz1[k] - vz2[l] ) )
#                        (PSIpointer_lower[0])[k][l]*=W[k][l]
                for k in range(4):
                    for l in range(4):
                        W[k][l] =exp(roughfac*cquadr( (vzpointer_lower[0])[k] - (vzpointer_upper[0])[l] ) )
                        (PSIpointer_lower[0])[k][l]*=W[k][l]
                Mult4x4_leftside(&P, PSIpointer_lower)

            Lower=MLCOMP[i][  MLLENGTH[i]-1 ]
            LowerLayer=LR[ Lower ]
            if(vzfilled==2):
                vzpointer_lower=&vz2
                vzpointer_upper=&vz1
                PHIpointer_lower=&PHI2
                PHIpointer_upper=&PHI1
                PSIpointer_lower=&PSI2
                PSIpointer_upper=&PSI1
                vzfilled=1

            else:
                vzpointer_lower=&vz1
                vzpointer_upper=&vz2
                PHIpointer_lower=&PHI1
                PHIpointer_upper=&PHI2
                PSIpointer_lower=&PSI1
                PSIpointer_upper=&PSI2
                vzfilled=2

            Mult4x4_leftside_diag(&P, &p)
            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
                UpperLayer=LR[ Upper ]
                Calculate_Phi_and_Psi(UpperLayer, &( AllMS[ Layer_type_to_Matrixsafe[Upper] ] ), vy, vzvz, vyvy, vzpointer_upper, PHIpointer_upper, PSIpointer_upper)
                Mult4x4_leftside(PSIpointer_lower, PHIpointer_upper)
            else:

                Mult4x4_leftside(PSIpointer_lower, &PHIvac)
                (vzpointer_upper[0])[0]=vzvac
                (vzpointer_upper[0])[1]=vzvac
                (vzpointer_upper[0])[2]=-vzvac
                (vzpointer_upper[0])[3]=-vzvac
            roughfac =-0.5*quadr( LowerLayer.Roughness*k0)
#            for k in range(4):
#                for l in range(4):
#                    W[k][l] =exp(roughfac*cquadr(vz1[k] - vz2[l] ) )
#                    (PSIpointer_lower[0])[k][l]*=W[k][l]
            for k in range(4):
                for l in range(4):
                    W[k][l] =exp(roughfac*cquadr( (vzpointer_lower[0])[k] - (vzpointer_upper[0])[l] ) )
                    (PSIpointer_lower[0])[k][l]*=W[k][l]

#            print("hallo2")
#            for k in range(4):
#                for l in range(4):
#                    print( (PSIpointer_lower[0])[k][l]  )
            Mult4x4_leftside(&P, PSIpointer_lower)

            for k in range(4):
                p[k] = exp(-ik0*UpperLayer.Thickness*(vzpointer_upper[0])[k] )


        i+=1

    divide=1.0/(P[0][1]*P[1][0]-P[1][1]*P[0][0])
 #   print( P[1][1], P[0][2], P[0][1], P[1][2], divide )
    (rtot[0])[0][0]=(P[1][1]*P[0][2]-P[0][1]*P[1][2])*divide # Incoming 1 reflected 1
    (rtot[0])[0][1]=-(P[1][1]*P[0][3]-P[0][1]*P[1][3])*divide # Incoming 2 reflected 1
    (rtot[0])[1][0]=-(P[0][0]*P[1][2]-P[1][0]*P[0][2])*divide # Incoming 1 reflected 2
    (rtot[0])[1][1]=(P[0][0]*P[1][3]-P[1][0]*P[0][3])*divide # Incoming 2 reflected 2
