import cython
from libc.stdlib cimport malloc, free
cimport numpy as np
import numpy as np
np.import_array()
import ctypes
from numpy cimport ndarray

from libc.math cimport  sin, cos
cdef extern from "complexobject.h":
    pass

cdef extern from "math.h":
    double complex exp(double complex) nogil
    double complex sqrt(double complex)  nogil


cdef struct CLayer:
    double Thickness, Roughness
    double complex ex, ey, ez, eg
    int type, magdir


cdef double two_times_pi=6.283185307179586
cdef double deg_to_rad=0.017453292522222223

cdef struct Heterostructure:
    int NLayers
    int NLayers_types
    int *MLLENGTH
    int *MLREP
    int **MLCOMP
    CLayer *LR



cdef class Lowestlayer:
    cdef CLayer Layercontent
    cdef Heterostructure Mother
    cdef int islowestlayer
    def __cinit__(self, double d, double sigma, MLstructure, int NL_types):
      #  print "init is called"
        self.Layercontent.Thickness = d
        self.Layercontent.Roughness = sigma
        self.Layercontent.type=1
        self.Layercontent.magdir = 0
        self.Layercontent.ex=1
        self.Layercontent.ey=1
        self.Layercontent.ez=1
        self.Layercontent.eg=0
        self.islowestlayer=1
        cdef int NLayers
        cdef int *MLLENGTH
        cdef int *MLREP
        cdef int **MLCOMP
        cdef CLayer *LR
        NLayers=FindLayerNumber(MLstructure, NL_types)
        MLLENGTH = <int*>malloc( NLayers* sizeof(int))
        MLREP = <int*>malloc( NLayers* sizeof(int))
        MLCOMP = <int**>malloc( NLayers* sizeof(int*))
        LR=<CLayer*>malloc( NL_types* sizeof(CLayer) )

        if(MLstructure=="default"):
            for i in range(NL_types):
                MLREP[i]=1
                MLLENGTH[i]=1
                MLCOMP[i]= <int*>malloc(sizeof(int))
                MLCOMP[i][0]=i
        else:
            SaveNumbers=[]

            a=MLstructure.split(",")

            N=len(a)
            j=0
            IsInBracket=False
            k=0
            for i in range(N):
                if (a[i].find('(')!=-1):
                    b=a[i].split('*(')
                    MLREP[j]=int(b[0])
                    SaveNumbers=SaveNumbers+[int(b[1])]
                    IsInBracket=True
                    k=1
                elif (a[i].find(')')!=-1):
                    b=a[i].split(')')
                    SaveNumbers=SaveNumbers+[int(b[0])]
                    IsInBracket=False
                    k=k+1
                    MLLENGTH[j]=k
                    MLCOMP[j]= <int*>malloc(int(MLLENGTH[j])*sizeof(int))
                    for l in range(MLLENGTH[j]):
                        MLCOMP[j][l]=int(SaveNumbers[l])
                    SaveNumbers=[]
                    j=j+1
                else:
                    if(IsInBracket):
                        k=k+1
                        SaveNumbers=SaveNumbers+[int(a[i])]
                    else:
                        MLREP[j]=1
                        MLLENGTH[j]=1
                        MLCOMP[j]= <int*>malloc(sizeof(int))
                        MLCOMP[j][0]=int(a[i])
                        j=j+1




        self.Mother.MLLENGTH=MLLENGTH
        self.Mother.MLREP=MLREP
        self.Mother.MLCOMP=MLCOMP
        self.Mother.NLayers=NLayers
        self.Mother.NLayers_types=NL_types
        self.Mother.LR=LR
        MakeConsistencyCheck(MLstructure, &(self.Mother), NL_types)

    def __dealloc__(self):
    #    print "dealloc is called"
        for i in range(self.Mother.NLayers):
            free(self.Mother.MLCOMP[i])
        free(self.Mother.MLCOMP)
        free(self.Mother.MLLENGTH)
        free(self.Mother.MLREP)
        free(self.Mother.LR)
    cpdef int isthisthelowestlayer(self):
        return self.islowestlayer
    cpdef long layercontentpointer(self):
        return <long>(&(self.Layercontent))
    cpdef long motherpointer(self):
        return <long>(&(self.Mother))
    cpdef double d(self):
        return self.Layercontent.Thickness
    cpdef setd(self, double d):
        self.Layercontent.Thickness=float(d)
    cpdef double sigma(self):
        return self.Layercontent.Roughness
    cpdef setsigma(self, double sigma):
        self.Layercontent.Roughness=float(sigma)
    def seteps(self, epsilon):
        if(hasattr(epsilon, "__len__")):
            if(len(epsilon) in (1,3,4)):

                if(len(epsilon)==1):
                    self.Layercontent.type=1
                    #self.Layercontent.eps[0]=epsilon[0]
                    self.Layercontent.ex=epsilon[0]
                    self.Layercontent.ey=epsilon[0]
                    self.Layercontent.ez=epsilon[0]
                    self.Layercontent.eg=0
                    self.Layercontent.magdir=0
                elif(len(epsilon)==3):
                    self.Layercontent.type=2
                    self.Layercontent.ex=epsilon[0]
                    self.Layercontent.ey=epsilon[1]
                    self.Layercontent.ez=epsilon[2]
                    self.Layercontent.eg=0
                    self.Layercontent.magdir=0
                elif(len(epsilon)==4):
                    self.Layercontent.type=3
                    self.Layercontent.ex=epsilon[0]
                    self.Layercontent.ey=epsilon[1]
                    self.Layercontent.ez=epsilon[2]
                    if(epsilon[3]==0):
                        self.Layercontent.magdir=0
                        self.Layercontent.type=2
                    elif(self.Layercontent.magdir==0):
                        raise Exception("Please use setmag to set a magnetization direction for this layer before you set a gyrotropy")
                    self.Layercontent.eg=epsilon[3]



            else:
                raise Exception('Epsilon must be a number or an array of length 1, 3 or 4')
        else:
            self.Layercontent.type=1
            self.Layercontent.ex=epsilon
            self.Layercontent.ey=epsilon
            self.Layercontent.ez=epsilon
    def eps(self, wl=-1):
        if(self.Layercontent.type==1):
            return self.Layercontent.ex
        elif(self.Layercontent.type==2):
            return np.array([ self.Layercontent.ex, self.Layercontent.ey, self.Layercontent.ez ])
        else:
            return np.array([ self.Layercontent.ex, self.Layercontent.ey, self.Layercontent.ez , self.Layercontent.eg ])
    cpdef double complex epsxx(self):
        return self.Layercontent.ex
    cpdef double complex epsyy(self):
        return self.Layercontent.ey
    cpdef double complex epszz(self):
        return self.Layercontent.ez
    cpdef double complex eg(self):
        return self.Layercontent.eg




    def setmag(self, dir):
        if(dir=='x'):
            self.Layercontent.magdir=1
        elif(dir=='y'):
            self.Layercontent.magdir=2
        elif(dir=='z'):
            self.Layercontent.magdir=3
        elif(dir=='0'):
            self.Layercontent.magdir=0
        else:
            raise Exception("Allowed input for setmag is 'x', 'y', 'z' or '0'")
    def mag(self):
        if(self.Layercontent.type==3):
            if(self.Layercontent.magdir==1):
                return "Magnetized along the x direction"
            elif(self.Layercontent.magdir==2):
                return "Magnetized along the y direction"
            elif(self.Layercontent.magdir==3):
                return "Magnetized along the z direction"
            else:
                return "No magnetization has been set"
        else:
            return "This layer is not magnetic"
    cdef int GetType(self):
        return self.Layercontent.type
    cdef int dir(self):
        return self.Layercontent.magdir
    def setepsxx(self, x):
        self.Layercontent.ex=x
    def setepsyy(self, x):
        self.Layercontent.ey=x
    def setepszz(self, x):
        self.Layercontent.ez=x
    def setepsg(self,  x):
        if(x==0):
            self.Layercontent.magdir=0
            self.Layercontent.type=2
        elif(self.Layercontent.magdir==0):
            raise Exception("Please use setmag to set a magnetization direction for this layer before you set a gyrotropy")
        self.Layercontent.eg=x
#    def __getstate__(self):
#        return list([self.Layercontent.Thickness, \
#        self.Layercontent.Roughness, \
#        self.Layercontent.type, \
#        self.Layercontent.magdir,\
#        self.Layercontent.ex,\
#        self.Layercontent.ey,\
#        self.Layercontent.ez,\
#        self.Layercontent.eg,\
#        self.islowestlayer ])
#    def __setstate__(self, x):
#        self.Layercontent.Thickness=x[0]
#        self.Layercontent.Roughness=x[1]
#        self.Layercontent.type=x[2]
#        self.Layercontent.magdir=x[3]
#        self.Layercontent.ex=x[4]
#        self.Layercontent.ey=x[5]
#        self.Layercontent.ez=x[6]
#        self.Layercontent.eg=x[7]
#        self.islowestlayer==x[8]

cdef class Layer:
    cdef CLayer Layercontent
    cdef int islowestlayer
    def __cinit__(self, double d, double sigma=0):
        self.Layercontent.Thickness = d
        self.Layercontent.Roughness = sigma
        self.Layercontent.type=1
        self.Layercontent.magdir = 0
        self.Layercontent.ex=1
        self.Layercontent.ey=1
        self.Layercontent.ez=1
        self.Layercontent.eg=0
        self.islowestlayer=0
    cpdef int isthisthelowestlayer(self):
        return self.islowestlayer
    cpdef long layercontentpointer(self):
        return <long>(&(self.Layercontent))
    cpdef double d(self):
        return self.Layercontent.Thickness
    cpdef setd(self, double d):
        self.Layercontent.Thickness=float(d)
    cpdef double sigma(self):
        return self.Layercontent.Roughness
    cpdef setsigma(self, double sigma):
        self.Layercontent.Roughness=float(sigma)
    def seteps(self, epsilon):
        if(hasattr(epsilon, "__len__")):
            if(len(epsilon) in (1,3,4)):
                if(len(epsilon)==1):
                    self.Layercontent.type=1
                    #self.Layercontent.eps[0]=epsilon[0]
                    self.Layercontent.ex=epsilon[0]
                    self.Layercontent.ey=epsilon[0]
                    self.Layercontent.ez=epsilon[0]
                    self.Layercontent.eg=0
                    self.Layercontent.magdir=0
                elif(len(epsilon)==3):
                    self.Layercontent.type=2
                    self.Layercontent.ex=epsilon[0]
                    self.Layercontent.ey=epsilon[1]
                    self.Layercontent.ez=epsilon[2]
                    self.Layercontent.eg=0
                    self.Layercontent.magdir=0
                elif(len(epsilon)==4):
                    self.Layercontent.type=3
                    self.Layercontent.ex=epsilon[0]
                    self.Layercontent.ey=epsilon[1]
                    self.Layercontent.ez=epsilon[2]
                    self.Layercontent.eg=epsilon[3]
                    if(epsilon[3]==0):
                        self.Layercontent.magdir=0
                        self.Layercontent.type=2
                    elif(self.Layercontent.magdir==0):
                        raise Exception("Please use setmag to set a magnetization direction for this layer before you set a gyrotropy")
            else:
                raise Exception('Epsilon must be a number or an array of length 1, 3 or 4')
        else:
            self.Layercontent.type=1
            self.Layercontent.ex=epsilon
            self.Layercontent.ey=epsilon
            self.Layercontent.ez=epsilon
    def eps(self):
        if(self.Layercontent.type==1):
            return self.Layercontent.ex
        elif(self.Layercontent.type==2):
            return np.array([ self.Layercontent.ex, self.Layercontent.ey, self.Layercontent.ez ])
        else:
            return np.array([ self.Layercontent.ex, self.Layercontent.ey, self.Layercontent.ez , self.Layercontent.eg ])
    cpdef double complex epsxx(self):
        return self.Layercontent.ex
    cpdef double complex epsyy(self):
        return self.Layercontent.ey
    cpdef double complex epszz(self):
        return self.Layercontent.ez
    cpdef double complex epsg(self):
        return self.Layercontent.eg

    def setmag(self, dir):
        if(dir=='x'):
            self.Layercontent.magdir=1
        elif(dir=='y'):
            self.Layercontent.magdir=2
        elif(dir=='z'):
            self.Layercontent.magdir=3
        elif(dir=='0'):
            self.Layercontent.magdir=0
        else:
            raise Exception("Allowed input for setmag is 'x', 'y', 'z' or '0'")
    def mag(self):
        if(self.Layercontent.type==3):
            if(self.Layercontent.magdir==1):
                return "Magnetized along the x direction"
            elif(self.Layercontent.magdir==2):
                return "Magnetized along the y direction"
            elif(self.Layercontent.magdir==3):
                return "Magnetized along the z direction"
            else:
                return "No magnetization has been set"
        else:
            return "This layer is not magnetic"
    cdef int GetType(self):
        return self.Layercontent.type
    cdef int dir(self):
        return self.Layercontent.magdir
    def setepsxx(self, x):
        self.Layercontent.ex=x
    def setepsyy(self, x):
        self.Layercontent.ey=x
    def setepszz(self, x):
        self.Layercontent.ez=x
    def setepsg(self,  x):
        if(x==0):
            self.Layercontent.magdir=0
            self.Layercontent.type=2
        elif(self.Layercontent.magdir==0):
            raise Exception("Please use setmag to set a magnetization direction for this layer before you set a gyrotropy")
        self.Layercontent.eg=x









cdef inline double complex cquadr(double complex x):
    return x*x
cdef inline double quadr(double x):
    return x*x
cdef inline double cabsquadr(double complex x):
    return quadr(x.real)+quadr(x.imag)

cdef double complex cconj(double complex x):
    return x.real -1j*x.imag

cdef inline double complex CalculateVZsigma(double vyvy, double complex ex):
    return sqrt(ex-vyvy)

cdef inline double complex CalculateVZpi(double vyvy, double complex ey, double complex ez):
    return sqrt((1.-vyvy/ez)*ey)

cdef inline double complex CalculateVZpi_m(double vyvy, double complex ey, double complex ez, double complex eg):
    return sqrt((1.-vyvy/ez)*ey+eg*eg/ez)


cdef void Calculate_Multilayer(double complex *t_comp1_up, double complex *t_comp2_up, double complex *t_comp1_do, double complex *t_comp2_do, double complex *r_ML_in1, double complex *r_ML_in2, double complex *r_ML_ba1, double complex *r_ML_ba2, double N):

    cdef double complex rres1, rres2, tres_up, tres_do, MLfac

    rres1=r_ML_in1[0]
    rres2=r_ML_ba1[0]
    tres_up=t_comp1_up[0]
    tres_do=t_comp1_do[0]
    if(N==0):
        t_comp2_up[0]=1
        t_comp2_do[0]=1
        r_ML_in2[0]=0
        r_ML_ba2[0]=0
        return
    N=N-1

    while(N):
        if(N%2): #Adding one Compound Layer to the heterostructure
            MLfac=1.0/(1-rres1*r_ML_ba1[0])
            rres1=r_ML_in1[0]+t_comp1_up[0]*t_comp1_do[0]*rres1*MLfac
            rres2=rres2+tres_up*tres_do*r_ML_ba1[0]*MLfac

            tres_up=t_comp1_up[0]*tres_up*MLfac
            tres_do=t_comp1_do[0]*tres_do*MLfac

        N=N/2


        #Doubling the heterostructure:
        MLfac=1.0/(1-r_ML_in1[0]*r_ML_ba1[0])
        r_ML_in1[0]=r_ML_in1[0]+t_comp1_up[0]*t_comp1_do[0]*r_ML_in1[0]*MLfac
        r_ML_ba1[0]=r_ML_ba1[0]+t_comp1_up[0]*t_comp1_do[0]*r_ML_ba1[0]*MLfac
        t_comp1_up[0]=cquadr(t_comp1_up[0])*MLfac
        t_comp1_do[0]=cquadr(t_comp1_do[0])*MLfac

    t_comp2_up[0]=tres_up
    t_comp2_do[0]=tres_do
    r_ML_in2[0]=rres1
    r_ML_ba2[0]=rres2



cdef double complex ipow(double complex base, double exp):
    cdef double complex result = 1.
    while (exp):
        if (exp%2): #If exp is Uneven
            result *= base
        exp=exp/2 #exp halbieren
        base=base*base
    return result



cdef double complex LinDicParatt_Sigma(Heterostructure* HS, double th, double wavelength) except *:
    if((th<=0)|(th>90)):
        raise Exception("Theta must be in the range 0<theta<=90")


    cdef double k0=two_times_pi/wavelength
    cdef double sintheta=sin(deg_to_rad*th)
    cdef double vyvy=quadr(cos(deg_to_rad*th))


    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR



    cdef int i,j

    cdef int Cap=NLAYERS-1
    cdef double complex rtot, rprime, vz, p,t,MSfac,pquad, rsl, vzupper, vzlower
    #cdef double complex vzlower

    cdef double complex r_ML_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_2 # Upper reflectivity from the Multilayer
    cdef double complex ptot
    cdef int Lower, Upper
    cdef double complex rough, groundrough
    cdef double roughfak=-2.*quadr(k0)



    vzlower=CalculateVZsigma(vyvy, LR[ MLCOMP[0][0] ].ex )



  #  print NLAYERS

    if(NLAYERS==1):
        vzupper=sintheta
    else:
        vzupper=CalculateVZsigma(vyvy,   LR[ MLCOMP[1][0] ].ex  )


    rough=exp(vzlower*vzupper*quadr(   LR[ MLCOMP[0][0] ].Roughness )*roughfak)
  #  print "all:", MLCOMP[0][0], MLCOMP[1][0], LR[ MLCOMP[0][0] ].ex, LR[ MLCOMP[1][0] ].ex, vzupper, vzlower
    rtot=(vzupper-vzlower)/(vzupper+vzlower)*rough
  #  print rtot
    i=1
    while i<NLAYERS:

        if(MLLENGTH[i]==1):
            vzlower=vzupper
            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
                vzupper=CalculateVZsigma(vyvy,  LR[ Upper ].ex  )
            else:
                vzupper=sintheta
            Lower=MLCOMP[i][0]

            rough=exp(roughfak*vzlower*vzupper*quadr(   LR[ Lower ].Roughness  ))

            rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough

            pquad=exp(2j*k0*( LR[ Lower ] ).Thickness*vzlower)
          #  print pquad, rprime
           # No Multiple scattering:
            rtot=(rprime+rtot*pquad)

        else:

            rough=exp(roughfak*vzlower*vzupper*quadr(  ( LR[ MLCOMP[i][0] ] ).Roughness ))

            vzlower=vzupper
            vzupper=CalculateVZsigma(vyvy, LR[ MLCOMP[i][1] ].ex)
            rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough


            r_ML_1=rprime
           # print "2"

            j=1
            p=1.

           # print MLLENGTH[i]
            while j<MLLENGTH[i]:
              #  print "3 ", j
                vzlower=vzupper
                pquad=exp(2j*k0* LR[MLCOMP[i][j]].Thickness*vzlower)
                p*=pquad
             #   print pquad, p
                if(j==(MLLENGTH[i]-1)):
                    Upper=MLCOMP[i][0]
                else:
                    Upper=MLCOMP[i][j+1]
                vzupper=CalculateVZsigma(vyvy,  LR[ Upper ].ex )

                rough=exp(roughfak*vzlower*vzupper*quadr(  LR[ MLCOMP[i][j] ].Roughness))

                rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough

                r_ML_1=rprime+r_ML_1*pquad
              #  print "5 ", j
             #   print "new rml", r_ML_1
                j+=1
            pquad=exp(2j*k0* LR[Upper].Thickness*vzupper)
            p*=pquad
          #  print pquad, p
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
                Upper=MLCOMP[i][j+1]
                vzupper=CalculateVZsigma(vyvy,  LR[ Upper ].ex )
                rough=exp(roughfak*vzlower*vzupper*quadr(  LR[ Lower ].Roughness )  )
                pquad=exp(2j*k0*( LR[ Lower ] ).Thickness*vzlower)
                rprime=(vzupper-vzlower)/(vzlower+vzupper)*rough
                rtot=(rprime+rtot*pquad)
                j+=1
            #Now the next layer begins:
            vzlower=vzupper
            Lower=MLCOMP[i][j]
            if(i!=Cap):
               # print "hier1"

                Upper=MLCOMP[i+1][0]
                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)

                vzupper=CalculateVZsigma(vyvy,LR[ Upper ].ex  )
                rough=exp(roughfak*vzlower*vzupper*quadr(  LR[ Lower ].Roughness )  )
                rprime=(vzupper-vzlower)/(vzlower+vzupper)*rough
                rtot=rprime+rtot*pquad
            else:

                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)
                rough=exp(roughfak*sintheta*vzupper*quadr(  LR[ Lower ].Roughness )  )
                rprime=(sintheta-vzlower)/(vzlower+sintheta)*rough
                rtot=rprime+rtot*pquad
        i=i+1
    return rtot


cdef double complex LinDicParatt_Pi(Heterostructure* HS, double th, double wavelength) except *:
    if((th<=0)|(th>90)):
        raise Exception("Theta must be in the range 0<theta<=90")

    cdef double k0=two_times_pi/wavelength
    cdef double sintheta=sin(deg_to_rad*th)
    cdef double vyvy=quadr(cos(deg_to_rad*th))

    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j

    Cap=NLAYERS-1
    cdef double complex rtot, rprime, vz, p,t,MSfac,pquad, rsl, vzupper, vzlower
    #cdef double complex vzlower

    cdef double complex r_ML_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_2 # Upper reflectivity from the Multilayer
    cdef double complex ptot
    cdef int Lower, Upper
    cdef double complex rough, groundrough
    cdef double roughfak=-2.*quadr(k0)

    #rough=1.

    vzlower=CalculateVZpi(vyvy, LR[MLCOMP[0][0]].ey, LR[MLCOMP[0][0]].ez)


    if(NLAYERS==1):
        vzupper=sintheta
        rough=exp(vzlower*vzupper*quadr(LR[MLCOMP[0][0]].Roughness)*roughfak)
        rtot=(vzupper*LR[MLCOMP[0][0]].ey-vzlower)/(vzupper*LR[MLCOMP[0][0]].ey+vzlower)*rough
    else:
        vzupper=CalculateVZpi(vyvy, LR[MLCOMP[1][0]].ey, LR[MLCOMP[1][0]].ez)
        rough=exp(vzlower*vzupper*quadr(LR[MLCOMP[0][0]].Roughness )*roughfak)
        rtot=(vzupper*LR[MLCOMP[0][0]].ey-vzlower*LR[MLCOMP[1][0]].ey)/(vzupper*LR[MLCOMP[0][0]].ey+vzlower*LR[MLCOMP[1][0]].ey)*rough



    i=1
    while i<NLAYERS:

        if(MLLENGTH[i]==1):
            vzlower=vzupper
            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
                vzupper=CalculateVZpi(vyvy, LR[Upper].ey, LR[Upper].ez)
                Lower=MLCOMP[i][0]
                rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[i][0]].Roughness))
                rprime=(vzupper*LR[Lower].ey-vzlower*LR[Upper].ey)/(vzupper*LR[Lower].ey+vzlower*LR[Upper].ey)*rough
            else:
                Lower=MLCOMP[i][0]
                rough=exp(roughfak*vzlower*sintheta*quadr(LR[MLCOMP[i][0]].Roughness))
                rprime=(sintheta*LR[Lower].ey-vzlower)/(sintheta*LR[Lower].ey+vzlower)*rough




            pquad=exp(2j*k0*LR[Lower].Thickness*vzlower)

           # No Multiple scattering:
            rtot=(rprime+rtot*pquad)

        else:
            rough=exp(roughfak*vzlower*vzupper*quadr(  ( LR[ MLCOMP[i][0] ] ).Roughness ))

            vzlower=vzupper
            vzupper=CalculateVZpi(vyvy, LR[ MLCOMP[i][1] ].ey, LR[ MLCOMP[i][1] ].ez)
           # rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough
            r_ML_1=(vzupper*LR[MLCOMP[i][0]].ey-vzlower*LR[MLCOMP[i][1]].ey)/(vzupper*LR[MLCOMP[i][0]].ey+vzlower*LR[MLCOMP[i][1]].ey)*rough

           # print "2"

            j=1
            p=1.


            while j<MLLENGTH[i]:
              #  print "3 ", j
                vzlower=vzupper
                pquad=exp(2j*k0* LR[MLCOMP[i][j]].Thickness*vzlower)

                p*=pquad

                if(j==(MLLENGTH[i]-1)):
                    Upper=MLCOMP[i][0]
                else:
                    Upper=MLCOMP[i][j+1]
                vzupper=CalculateVZpi(vyvy, LR[ Upper ].ey, LR[ Upper ].ez)

                rough=exp(roughfak*vzlower*vzupper*quadr(  LR[ MLCOMP[i][j] ].Roughness))

              #  rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough
                rprime=(vzupper*LR[MLCOMP[i][j]].ey-vzlower*LR[Upper].ey)/(vzupper*LR[MLCOMP[i][j]].ey+vzlower*LR[Upper].ey)*rough
                r_ML_1=rprime+r_ML_1*pquad
              #  print "5 ", j
             #   print "new rml", r_ML_1
                j+=1
            pquad=exp(2j*k0* LR[Upper].Thickness*vzupper)
            p*=pquad
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
                Upper=MLCOMP[i][j+1]
                vzupper=CalculateVZpi(vyvy, LR[ Upper ].ey, LR[ Upper ].ez)
                rough=exp(roughfak*vzlower*vzupper*quadr(  LR[ Lower ].Roughness )  )
                pquad=exp(2j*k0*( LR[ Lower ] ).Thickness*vzlower)
                rprime=(vzupper*LR[Lower].ey-vzlower*LR[Upper].ey)/(vzupper*LR[Lower].ey+vzlower*LR[Upper].ey)*rough
                rtot=(rprime+rtot*pquad)
                j+=1
            #Now the next layer begins:
            vzlower=vzupper
            Lower=MLCOMP[i][j]
            if(i!=Cap):
               # print "hier1"

                Upper=MLCOMP[i+1][0]
                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)

                vzupper=CalculateVZpi(vyvy, LR[ Upper ].ey, LR[ Upper ].ez)
                rough=exp(roughfak*vzlower*vzupper*quadr(  LR[ Lower ].Roughness )  )
                rprime=(vzupper*LR[Lower].ey-vzlower*LR[Upper].ey)/(vzupper*LR[Lower].ey+vzlower*LR[Upper].ey)*rough
                rtot=rprime+rtot*pquad
            else:

                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)
                rough=exp(roughfak*sintheta*vzupper*quadr(  LR[ Lower ].Roughness )  )
                rprime=(sintheta*LR[Lower].ey-vzlower)/(sintheta*LR[Lower].ey+vzlower)*rough
                rtot=rprime+rtot*pquad

        i=i+1
    return rtot

cdef double complex LinDicParatt_Sigma_MS(Heterostructure* HS, double th, double wavelength) except *:
    if((th<0)|(th>90)):
        raise Exception("Theta must be in the range 0<=theta<=90")
    if(th==0):
        return 1.0
    cdef double k0=two_times_pi/wavelength
    cdef double sintheta=sin(deg_to_rad*th)
    cdef double vyvy=quadr(cos(deg_to_rad*th))

    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
    cdef double roughfak=-2.*quadr(k0)
    cdef double complex rough, rough2
    Cap=NLAYERS-1
    cdef double complex rtot, rprime, p,t,MSfac,pquad, rsl, vzupper, vzlower
    #cdef double complex vzlower
    cdef double complex t_comp1_down #Transmission through one Compound layer, down
    cdef double complex t_comp1_up #Transmission through one Compound layer, up
    cdef double complex t_comp2_down #Transmission through one Compound layer, down
    cdef double complex t_comp2_up #Transmission through one Compound layer, up
    cdef double complex r_ML_in_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_in_2 # Upper reflectivity from the Multilayer
    cdef double complex r_ML_back_1 # Back reflectivity from one Compound layer
    cdef double complex r_ML_back_2 # Back reflectivity from the Multilayer
    cdef int Lower, Upper

    vzlower=CalculateVZsigma(vyvy, LR[MLCOMP[0][0]].ex)

    if(NLAYERS==1):
        vzupper=sintheta
    else:
        vzupper=CalculateVZsigma(vyvy, LR[MLCOMP[1][0]].ex)

    rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[0][0]].Roughness))


    rtot=(vzupper-vzlower)/(vzupper+vzlower)*rough

    i=1
    while i<NLAYERS:
        if(MLLENGTH[i]==1):
            vzlower=vzupper
            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
                vzupper=CalculateVZsigma(vyvy, LR[Upper].ex)
            else:
                vzupper=sintheta
            Lower=MLCOMP[i][0]
            rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[i][0]].Roughness))

            rprime=(vzupper-vzlower)/(vzupper+vzlower)*rough
            pquad=exp(2j*k0*LR[Lower].Thickness*vzlower)
           # if(MultipleScattering):
            rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
           # else:
            #    rtot=rprime+rtot*pquad

        else:
            vzlower=vzupper
            vzupper=CalculateVZsigma(vyvy, LR[MLCOMP[i][1]].ex)


            r_ML_in_1=(vzupper-vzlower)/(vzupper+vzlower)
            rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[i][0]].Roughness))
            r_ML_in_1*=rough

            t_comp1_up=1-r_ML_in_1

            j=1
            while j<MLLENGTH[i]:
                if(j+1<MLLENGTH[i]):
                    Upper=MLCOMP[i][j+1]
                else:
                    Upper=MLCOMP[i][0]
                Lower=MLCOMP[i][j]
                vzlower=vzupper
                vzupper=CalculateVZsigma(vyvy, LR[Upper].ex)
                rough=exp(roughfak*vzupper*vzlower*quadr(LR[Lower].Roughness))
                rprime=rough*(vzupper-vzlower)/(vzupper+vzlower)

                p=exp(1j*k0*LR[Lower].Thickness*vzlower)
                pquad=cquadr(p)
                MSfac=1.0/(1+rprime*r_ML_in_1*pquad)
                t=1-rprime
                t_comp1_up*=p*t*MSfac
                r_ML_in_1=(rprime+r_ML_in_1*pquad)*MSfac
                j+=1
            p=exp(1j*k0*LR[Upper].Thickness*vzupper)
            t_comp1_up*=p

            r_ML_back_1=-rprime
            t_comp1_down=1-r_ML_back_1
            j=MLLENGTH[i]-2
            while j>=0:
                vzupper=vzlower
                Upper=MLCOMP[i][j+1]
                Lower=MLCOMP[i][j]
                vzlower=CalculateVZsigma(vyvy, LR[Lower].ex)
                rough=exp(roughfak*vzupper*vzlower*quadr(LR[Lower].Roughness))

                rprime=rough*(vzlower-vzupper)/(vzupper+vzlower)
                t=1-rprime
                p=exp(1j*k0*LR[Upper].Thickness*vzupper)
                pquad=cquadr(p)

                MSfac=1.0/(1+rprime*r_ML_back_1*pquad)
                t_comp1_down*=p*t*MSfac

                r_ML_back_1=(rprime+r_ML_back_1*pquad)*MSfac
                j-=1
            p=exp(1j*k0*LR[Lower].Thickness*vzlower)
            t_comp1_down*=p
            r_ML_back_1*=cquadr(p)


            Calculate_Multilayer(&t_comp1_up, &t_comp2_up,&t_comp1_down, &t_comp2_down, &r_ML_in_1, &r_ML_in_2, &r_ML_back_1, &r_ML_back_2, MLREP[i]-1)



            rtot=r_ML_in_2+t_comp2_up*t_comp2_down*rtot/(1.-r_ML_back_2*rtot)
            #Now the next layer begins:
            vzupper=vzlower #CalculateVZsigma(vyvy, LR[MLCOMP[i][0]].ex)
            j=0
            while(j<MLLENGTH[i]-1):
                vzlower=vzupper
                Lower=MLCOMP[i][j]
                Upper=MLCOMP[i][j+1]
                vzupper=CalculateVZsigma(vyvy,  LR[ Upper ].ex )
                rough=exp(roughfak*vzlower*vzupper*quadr(  LR[ Lower ].Roughness )  )
                pquad=exp(2j*k0*( LR[ Lower ] ).Thickness*vzlower)
                rprime=(vzupper-vzlower)/(vzlower+vzupper)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
                j+=1
            #Now the next layer begins:
            vzlower=vzupper
            Lower=MLCOMP[i][j]
            if(i!=Cap):
               # print "hier1"

                Upper=MLCOMP[i+1][0]
                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)

                vzupper=CalculateVZsigma(vyvy,LR[ Upper ].ex  )
                rough=exp(roughfak*vzlower*vzupper*quadr(  LR[ Lower ].Roughness )  )
                rprime=(vzupper-vzlower)/(vzlower+vzupper)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
            else:

                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)
                rough=exp(roughfak*sintheta*vzlower*quadr(  LR[ Lower ].Roughness )  )
                rprime=(sintheta-vzlower)/(vzlower+sintheta)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)


        i=i+1


    return rtot


cdef double complex LinDicParatt_Pi_MS(Heterostructure* HS, double th, double wavelength) except *:
    if((th<0)|(th>90)):
        raise Exception("Theta must be in the range 0<=theta<=90")
    if(th==0):
        return 1.0

    cdef double k0=two_times_pi/wavelength
    cdef double sintheta=sin(deg_to_rad*th)
    cdef double vyvy=quadr(cos(deg_to_rad*th))


    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
    cdef double roughfak=-2.*quadr(k0)
    cdef double complex rough, rough2
    Cap=NLAYERS-1
    cdef double complex rtot, rprime, vz, p,t,MSfac,pquad, rsl, vzupper, vzlower
    #cdef double complex vzlower
    cdef double complex t_comp1_down #Transmission through one Compound layer, down
    cdef double complex t_comp1_up #Transmission through one Compound layer, up
    cdef double complex t_comp2_down #Transmission through one Compound layer, down
    cdef double complex t_comp2_up #Transmission through one Compound layer, up
    cdef double complex r_ML_in_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_in_2 # Upper reflectivity from the Multilayer
    cdef double complex r_ML_back_1 # Back reflectivity from one Compound layer
    cdef double complex r_ML_back_2 # Back reflectivity from the Multilayer
    cdef int Lower, Upper

    vzlower=CalculateVZpi(vyvy, LR[MLCOMP[0][0]].ey, LR[MLCOMP[0][0]].ez)

    if(NLAYERS==1):
        vzupper=sintheta
        rough=exp(vzlower*vzupper*quadr(LR[MLCOMP[0][0]].Roughness)*roughfak)
        rtot=(vzupper*LR[MLCOMP[0][0]].ey-vzlower)/(vzupper*LR[MLCOMP[0][0]].ey+vzlower)*rough
    else:
        vzupper=CalculateVZpi(vyvy, LR[MLCOMP[1][0]].ey, LR[MLCOMP[1][0]].ez)
        rough=exp(vzlower*vzupper*quadr(LR[MLCOMP[0][0]].Roughness)*roughfak)
        rtot=(vzupper*LR[MLCOMP[0][0]].ey-vzlower*LR[MLCOMP[1][0]].ey)/(vzupper*LR[MLCOMP[0][0]].ey+vzlower*LR[MLCOMP[1][0]].ey)*rough



    i=1
    while i<NLAYERS:
        if(MLLENGTH[i]==1):
            vzlower=vzupper
            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
                Lower=MLCOMP[i][0]
                vzupper=CalculateVZpi(vyvy, LR[Upper].ey, LR[Upper].ez)
                rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[i][0]].Roughness))
                rprime=(vzupper*LR[Lower].ey-vzlower*LR[Upper].ey)/(vzupper*LR[Lower].ey+vzlower*LR[Upper].ey)*rough
            else:
                vzupper=sintheta
                Lower=MLCOMP[i][0]
                rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[i][0]].Roughness))
                rprime=(vzupper*LR[Lower].ey-vzlower)/(vzupper*LR[Lower].ey+vzlower)*rough




            pquad=exp(2j*k0*LR[Lower].Thickness*vzlower)
           # if(MultipleScattering):
            rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
           # else:
            #    rtot=rprime+rtot*pquad

        else:

            vzlower=vzupper
            vzupper=CalculateVZpi(vyvy, LR[MLCOMP[i][1]].ey, LR[MLCOMP[i][1]].ez )


            r_ML_in_1=(vzupper*LR[MLCOMP[i][0]].ey-LR[MLCOMP[i][1]].ey*vzlower)/(vzupper*LR[MLCOMP[i][0]].ey+LR[MLCOMP[i][1]].ey*vzlower)
            rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[i][0]].Roughness))
            r_ML_in_1*=rough

            t_comp1_up=1-r_ML_in_1
            j=1
            while j<MLLENGTH[i]:
                if(j+1<MLLENGTH[i]):
                    Upper=MLCOMP[i][j+1]
                else:
                    Upper=MLCOMP[i][0]
                Lower=MLCOMP[i][j]
                vzlower=vzupper
                vzupper=CalculateVZpi(vyvy, LR[Upper].ey , LR[Upper].ez)
                rough=exp(roughfak*vzupper*vzlower*quadr(LR[Lower].Roughness))
                rprime=rough*(vzupper*LR[Lower].ey-vzlower*LR[Upper].ey)/(vzupper*LR[Lower].ey+vzlower*LR[Upper].ey)

                p=exp(1j*k0*LR[Lower].Thickness*vzlower)
                pquad=cquadr(p)
                MSfac=1.0/(1+rprime*r_ML_in_1*pquad)
                t=1-rprime
                t_comp1_up*=p*t*MSfac
                r_ML_in_1=(rprime+r_ML_in_1*pquad)*MSfac
                j+=1
            p=exp(1j*k0*LR[Upper].Thickness*vzupper)
            t_comp1_up*=p

            r_ML_back_1=-rprime
            t_comp1_down=1-r_ML_back_1

            j=MLLENGTH[i]-2
            while j>=0:
                vzupper=vzlower
                Upper=MLCOMP[i][j+1]
                Lower=MLCOMP[i][j]
                vzlower=CalculateVZpi(vyvy, LR[Lower].ey , LR[Lower].ez)
                rough=exp(roughfak*vzupper*vzlower*quadr(LR[Lower].Roughness))

                rprime=rough*(vzlower*LR[Upper].ey-vzupper*LR[Lower].ey)/(vzupper*LR[Lower].ey+vzlower*LR[Upper].ey)
                t=1-rprime
                p=exp(1j*k0*LR[Upper].Thickness*vzupper)
                pquad=cquadr(p)

                MSfac=1.0/(1+rprime*r_ML_back_1*pquad)
                t_comp1_down*=p*t*MSfac

                r_ML_back_1=(rprime+r_ML_back_1*pquad)*MSfac
                j-=1
            p=exp(1j*k0*LR[Lower].Thickness*vzlower)
            t_comp1_down*=p
            r_ML_back_1*=cquadr(p)


            Calculate_Multilayer(&t_comp1_up, &t_comp2_up,&t_comp1_down, &t_comp2_down, &r_ML_in_1, &r_ML_in_2, &r_ML_back_1, &r_ML_back_2, MLREP[i]-1)



            rtot=r_ML_in_2+t_comp2_up*t_comp2_down*rtot/(1.-r_ML_back_2*rtot)
            #Now the next layer begins:
            vzupper=vzlower #CalculateVZsigma(vyvy, LR[MLCOMP[i][0]].ex)
            j=0
            while(j<MLLENGTH[i]-1):
                vzlower=vzupper
                Lower=MLCOMP[i][j]
                Upper=MLCOMP[i][j+1]
                vzupper=CalculateVZpi(vyvy, LR[Upper].ey , LR[Upper].ez)
                rough=exp(roughfak*vzlower*vzupper*quadr(  LR[ Lower ].Roughness )  )
                pquad=exp(2j*k0*( LR[ Lower ] ).Thickness*vzlower)
                rprime=(vzupper*LR[Lower].ey-vzlower*LR[Upper].ey)/(vzupper*LR[Lower].ey+vzlower*LR[Upper].ey)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
                j+=1
            #Now the next layer begins:
            vzlower=vzupper
            Lower=MLCOMP[i][j]
            if(i!=Cap):
               # print "hier1"

                Upper=MLCOMP[i+1][0]
                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)

                vzupper=CalculateVZpi(vyvy, LR[Upper].ey , LR[Upper].ez)
                rough=exp(roughfak*vzlower*vzupper*quadr(  LR[ Lower ].Roughness )  )
                rprime=(vzupper*LR[Lower].ey-vzlower*LR[Upper].ey)/(vzupper*LR[Lower].ey+vzlower*LR[Upper].ey)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)
            else:

                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)
                rough=exp(roughfak*sintheta*vzupper*quadr(  LR[ Lower ].Roughness )  )
                rprime=(sintheta*LR[Lower].ey-vzlower)/(sintheta*LR[Lower].ey+vzlower)*rough
                rtot=(rprime+rtot*pquad)/(1+rprime*rtot*pquad)



        i=i+1


    return rtot



cdef void Relevant_Stuff_for_xmag(double complex ey1, double complex ey2, double complex ez1, double complex ez2,\
                                  double complex eg1, double complex eg2, double complex vz1, double complex vz2, \
                                  double vy, double k0, double sigma, \
                                  double complex *r, double complex *rp, double complex *t, double complex *tp):
    cdef double vyvy=vy*vy
    cdef double complex div1
    cdef double complex a,b,c,d,e
    div1=-vy*eg1/ez1
    a=vz1+div1
    b=vz1-div1
    c=1.-vyvy/ez1
    div1=ez2-vyvy
    d=(ez2*vz2+eg2*vy)/div1
    e=(-ez2*vz2+eg2*vy)/div1
    div1=1./2.*vz1
    cdef double complex J11,J12,J21,J22
    cdef double roughfac=-0.5*quadr(sigma)*quadr(k0)
    cdef double complex roughplus=exp(roughfac*cquadr(vz1+vz2))
    cdef double complex roughminus=exp(roughfac*cquadr(vz1-vz2))

    J11=roughminus*(a+c*d)*div1
    J12=roughplus*(a+c*e)*div1
    J21=roughplus*(b-c*d)*div1
    J22=roughminus*(b-c*e)*div1

    r[0]=(-J12/J11)
    t[0]=J21*r[0]+J22
    rp[0]=J21/J11
    tp[0]=1./J11



cdef double complex LinDicParatt_Pi_xmag(Heterostructure* HS, double th, double wavelength) except *:
    if((th<=0)|(th>90)):
        raise Exception("Theta must be in the range 0<theta<=90")


    cdef double k0=two_times_pi/wavelength
    cdef double sintheta=sin(deg_to_rad*th)
    cdef double vy=cos(deg_to_rad*th)
    cdef double vyvy=quadr(vy)


    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j

    Cap=NLAYERS-1
    cdef double complex rtot, r, rprime, vz, p,t,tp,ttot,MSfac,pquad, rsl, vzupper, vzlower
    #cdef double complex vzlower

    cdef double complex r_ML_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_2 # Upper reflectivity from the Multilayer
    cdef double complex ptot
    cdef int Lower, Upper
   # cdef double complex rough, groundrough
  #  cdef double roughfak=-2.*quadr(k0)

    #rough=1.

    vzlower=CalculateVZpi_m(vyvy, LR[MLCOMP[0][0]].ey, LR[MLCOMP[0][0]].ez, LR[MLCOMP[0][0]].eg)


    if(NLAYERS==1):
        vzupper=sintheta
       # rough=exp(vzlower*vzupper*quadr(LR[MLCOMP[0][0]].Roughness)*roughfak)
        Relevant_Stuff_for_xmag(LR[MLCOMP[0][0]].ey, 1., LR[MLCOMP[0][0]].ez, 1, \
                                LR[MLCOMP[0][0]].eg, 0., vzlower, sintheta, \
                                vy, k0, LR[MLCOMP[0][0]].Roughness, \
                                &rtot, &rprime, &t, &tp)
    else:
        vzupper=CalculateVZpi(vyvy, LR[MLCOMP[1][0]].ey, LR[MLCOMP[1][0]].ez)
     #   rough=exp(vzlower*vzupper*quadr(LR[MLCOMP[0][0]].Roughness)*roughfak)
        Relevant_Stuff_for_xmag(LR[MLCOMP[0][0]].ey, LR[MLCOMP[1][0]].ey, LR[MLCOMP[0][0]].ez, LR[MLCOMP[1][0]].ez, \
                                LR[MLCOMP[0][0]].eg, LR[MLCOMP[1][0]].eg, vzlower, vzupper, \
                                vy, k0, LR[MLCOMP[0][0]].Roughness, \
                                &rtot, &rprime, &t, &tp)
    #rtot=rtot*rough


    i=1
    while i<NLAYERS:

        if(MLLENGTH[i]==1):
            vzlower=vzupper
            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
              #  print Upper
                vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
                Lower=MLCOMP[i][0]
              #  rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[i][0]].Roughness))
                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[MLCOMP[i][0]].Roughness, \
                                &r, &rprime, &t, &tp)
              #  print "endif"
            else:
                vzupper=sintheta
                Lower=MLCOMP[i][0]
             #   rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[i][0]].Roughness))
                Relevant_Stuff_for_xmag(LR[Lower].ey, 1, LR[Lower].ez, 1, \
                                LR[Lower].eg, 0, vzlower, sintheta, \
                                vy, k0, LR[MLCOMP[i][0]].Roughness, \
                                &r, &rprime, &t, &tp)
              #  print "endelse"
           # print "r and rprime: ", r, rprime



            pquad=exp(2j*k0*LR[Lower].Thickness*vzlower)
           # print "Layer ", Upper, " - ", Lower
         #   print "r", r
          #  print "t and tp", t, tp
       #     print "t*tp", t*tp
          #  print "pquad", pquad
          #  print "rtot", rtot
          #  print "together", t*tp*pquad*rtot
           # No Multiple scattering:
           # rtot=(r*rough+t*tp*pquad*rtot)
            rtot=r+t*tp*pquad*rtot
        else:

            vzlower=vzupper
            Upper=MLCOMP[i][1]
            Lower=MLCOMP[i][0]
            vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)

            Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper,  \
                                vy, k0, LR[Lower].Roughness, \
                                &r_ML_1, &rprime, &t, &tp)

            j=1
            p=1.
            ttot=t*tp

            while j<MLLENGTH[i]:
              #  print "3 ", j
                vzlower=vzupper
                pquad=exp(2j*k0*LR[MLCOMP[i][j]].Thickness*vzlower)
                Lower=MLCOMP[i][j]
                p*=pquad
                if(j==(MLLENGTH[i]-1)):

                    Upper=MLCOMP[i][0]
                    vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)

                    Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[Lower].Roughness, \
                                &r, &rprime, &t, &tp)
                else:
                    Upper=MLCOMP[i][j+1]
                    vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)

                    Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[Lower].Roughness, \
                                &r, &rprime, &t, &tp)
              #  r=r*rough
                ttot*=t*tp
                r_ML_1=r+t*tp*r_ML_1*pquad

                j=j+1
            pquad=exp(2j*k0*LR[Upper].Thickness*vzupper)
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
            rtot=r_ML_2+ptot*rtot

            j=0
            while(j<MLLENGTH[i]-1):
                vzlower=vzupper
                Lower=MLCOMP[i][j]
                Upper=MLCOMP[i][j+1]
                vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[Lower].Roughness, \
                                &r, &rprime, &t, &tp)
                pquad=exp(2j*k0*( LR[ Lower ] ).Thickness*vzlower)
                rtot=(r+rtot*pquad)
                j+=1
            #Now the next layer begins:
            vzlower=vzupper
            Lower=MLCOMP[i][j]
            if(i!=Cap):
               # print "hier1"

                Upper=MLCOMP[i+1][0]
                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)

                vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[Lower].Roughness, \
                                &r, &rprime, &t, &tp)
                rtot=r+rtot*pquad
            else:

                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)
                Relevant_Stuff_for_xmag(LR[Lower].ey, 1, LR[Lower].ez, 1, \
                                LR[Lower].eg, 0, vzlower, sintheta, \
                                vy, k0, LR[Lower].Roughness, \
                                &r, &rprime, &t, &tp)
                rtot=r+rtot*pquad


        i=i+1
   # print "8 "
    return rtot


cdef double complex LinDicParatt_Pi_xmag_MS(Heterostructure* HS, double th, double wavelength) except *:
    if((th<=0)|(th>=90)):
        raise Exception("Theta must be in the range 0<theta<90")


    cdef double k0=two_times_pi/wavelength
    cdef double sintheta=sin(deg_to_rad*th)
    cdef double vy=cos(deg_to_rad*th)
    cdef double vyvy=quadr(vy)


    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
   # cdef double roughfak=-2.*quadr(k0)
    cdef double complex rough, rough2
    Cap=NLAYERS-1
    cdef double complex rtot, rprime, vz,r, p,t,tp,MSfac,pquad, rsl, vzupper, vzlower
    #cdef double complex vzlower
    cdef double complex t_ML_in_1 #Transmission through one Compound layer, down
    cdef double complex t_ML_back_1 #Transmission through one Compound layer, up
    cdef double complex t_ML_in_2 #Transmission through all Compound layer, down
    cdef double complex t_ML_back_2 #Transmission through all Compound layer, up
    cdef double complex r_ML_in_1 # Upper reflectivity from one Compound layer
    cdef double complex r_ML_in_2 # Upper reflectivity from the Multilayer
    cdef double complex r_ML_back_1 # Back reflectivity from one Compound layer
    cdef double complex r_ML_back_2 # Back reflectivity from the Multilayer
    cdef int Lower, Upper
    #rough=1.
    Lower=MLCOMP[0][0]
    Upper=MLCOMP[1][0]
    vzlower=CalculateVZpi_m(vyvy, LR[Lower].ey, LR[Lower].ez, LR[Lower].eg)

    if(NLAYERS==1):
        vzupper=sintheta
       # rough=exp(vzlower*vzupper*quadr(LR[Lower].Roughness)*roughfak)
        Relevant_Stuff_for_xmag(LR[Lower].ey, 1., LR[Lower].ez, 1, \
                                LR[Lower].eg, 0., vzlower, sintheta, \
                                vy, k0, LR[Lower].Roughness, \
                                &rtot, &rprime, &t, &tp)
    else:
        vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
     #   rough=exp(vzlower*vzupper*quadr(LR[Lower].Roughness)*roughfak)
        Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[Lower].Roughness, \
                                &rtot, &rprime, &t, &tp)
    #rtot=rtot*rough
   # print rtot

    i=1
    while i<NLAYERS:

        if(MLLENGTH[i]==1):
            vzlower=vzupper
            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
              #  print Upper
                vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
                Lower=MLCOMP[i][0]
              #  rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[i][0]].Roughness))
                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[MLCOMP[i][0]].Roughness, \
                                &r, &rprime, &t, &tp)
              #  print "endif"
            else:
                vzupper=sintheta
                Lower=MLCOMP[i][0]
             #   rough=exp(roughfak*vzlower*vzupper*quadr(LR[MLCOMP[i][0]].Roughness))
                Relevant_Stuff_for_xmag(LR[Lower].ey, 1, LR[Lower].ez, 1, \
                                LR[Lower].eg, 0, vzlower, sintheta, \
                                vy, k0, LR[MLCOMP[i][0]].Roughness, \
                                &r, &rprime, &t, &tp)
              #  print "endelse"
           # print "r and rprime: ", r, rprime



            pquad=exp(2j*k0*LR[Lower].Thickness*vzlower)
            rtot=r+t*tp*pquad*rtot/(1.-rtot*rprime*pquad)
          #  print rtot, vzlower, vzupper
        else:

            vzlower=vzupper

            Lower=MLCOMP[i][0]
            Upper=MLCOMP[i][1]
            vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
            Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[Lower].Roughness, \
                                &r_ML_in_1, &rprime, &t, &t_ML_back_1)
            j=1
            while j<MLLENGTH[i]:
                if(j+1<MLLENGTH[i]):
                    Upper=MLCOMP[i][j+1]
                else:
                    Upper=MLCOMP[i][0]
                Lower=MLCOMP[i][j]

                vzlower=vzupper
                vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[Lower].Roughness, \
                                &r, &rprime, &t, &tp)
                p=exp(1j*k0*LR[Lower].Thickness*vzlower)
                pquad=cquadr(p)


                MSfac=1.0/(1.-rprime*r_ML_in_1*pquad)

                t_ML_back_1=t_ML_back_1*p*tp*MSfac
                r_ML_in_1=r+t*tp*pquad*r_ML_in_1*MSfac

                j+=1
            p=exp(1j*k0*LR[Upper].Thickness*vzupper)
            t_ML_back_1*=p


            r_ML_back_1=rprime
            t_ML_in_1=t




            j=MLLENGTH[i]-2
            while j>=0:
                vzupper=vzlower
                Upper=MLCOMP[i][j+1]
                Lower=MLCOMP[i][j]
                vzlower=CalculateVZpi_m(vyvy, LR[Lower].ey, LR[Lower].ez, LR[Lower].eg)
                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[Lower].Roughness, \
                                &r, &rprime, &t, &tp)
                p=exp(1j*k0*LR[Upper].Thickness*vzupper)
                pquad=cquadr(p)

                MSfac=1.0/(1.-r*r_ML_back_1*pquad)

                t_ML_in_1*=t*p*MSfac
                r_ML_back_1=rprime+t*tp*pquad*r_ML_back_1*MSfac
                j-=1

            p=exp(1j*k0*LR[Lower].Thickness*vzlower)
            t_ML_in_1*=p
            r_ML_back_1*=cquadr(p)


            Calculate_Multilayer(&t_ML_back_1, &t_ML_back_2,&t_ML_in_1, &t_ML_in_2, &r_ML_in_1, &r_ML_in_2, &r_ML_back_1, &r_ML_back_2, MLREP[i]-1)



            rtot=r_ML_in_2+t_ML_back_2*t_ML_in_2*rtot/(1.-r_ML_back_2*rtot)
            #Now the next layer begins:

            Lower=MLCOMP[i][0]
            vzlower=CalculateVZpi_m(vyvy, LR[Lower].ey, LR[Lower].ez, LR[Lower].eg)

            vzupper=vzlower #CalculateVZsigma(vyvy, LR[MLCOMP[i][0]].ex)
            j=0
            while(j<MLLENGTH[i]-1):
                vzlower=vzupper
                Upper=MLCOMP[i][j+1]
                Lower=MLCOMP[i][j]
                vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[Lower].Roughness, \
                                &r, &rprime, &t, &tp)
                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)

                rtot=r+t*tp*pquad*rtot/(1.-rtot*rprime*pquad)
             #   print rtot, vzlower, vzupper
                j+=1
            #Now the next layer begins:
            vzlower=vzupper
            Lower=MLCOMP[i][j]
            if(i!=Cap):
               # print "hier1"

                Upper=MLCOMP[i+1][0]
                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)

                vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
                Relevant_Stuff_for_xmag(LR[Lower].ey, LR[Upper].ey, LR[Lower].ez, LR[Upper].ez, \
                                LR[Lower].eg, LR[Upper].eg, vzlower, vzupper, \
                                vy, k0, LR[Lower].Roughness, \
                                &r, &rprime, &t, &tp)
                rtot=r+t*tp*pquad*rtot/(1.-rtot*rprime*pquad)
            else:

                pquad=exp(2j*k0*LR[ Lower ].Thickness*vzlower)
                vzupper=CalculateVZpi_m(vyvy, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg)
                Relevant_Stuff_for_xmag(LR[Lower].ey, 1, LR[Lower].ez, 1, \
                                LR[Lower].eg, 0, vzlower, sintheta, \
                                vy, k0, LR[Lower].Roughness, \
                                &r, &rprime, &t, &tp)
                rtot=r+t*tp*pquad*rtot/(1.-rtot*rprime*pquad)
          #  print rtot, vzlower, vzupper
        i=i+1

    return rtot



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


cdef void Calculate_rt(double complex (*PSI1)[4], double complex (*PHI1)[4], double complex (*PSI2)[4], double complex (*PHI2)[4], double complex (*r)[2][2], double complex (*rprime)[2][2], double complex (*t)[2][2], double complex (*tprime)[2][2], int Magnetic, \
                       double complex vz1, double complex vz2, double complex vz3, double complex vz4, double sigma, double k0):

    cdef double complex J[2][4]
    cdef double complex b,d, div;
    cdef double roughfac=-0.5*quadr(sigma*  k0)
  #  cdef double complex roughplus=exp(roughfac*cquadr(vz1+vz2))
  #  cdef double complex roughminus=exp(roughfac*cquadr(vz1-vz2))

    if(Magnetic):



        b=(PHI1[0])[3]*(PSI1[0])[3]+(PSI1[0])[0]
        d=(PHI1[0])[2]*(PSI1[0])[2]+(PSI1[0])[1]*(PHI1[0])[1]
        J[0][0]=(b+d)*exp(roughfac*cquadr(vz1-vz3))
        J[0][2]=(b-d)*exp(roughfac*cquadr(vz1+vz3))



        b=(PHI2[0])[0]*(PSI2[0])[0]+(PSI2[0])[3]
        d=(PHI2[0])[2]*(PSI2[0])[2]+(PSI2[0])[1]*(PHI2[0])[1]
        J[1][1]=(b+d)*exp(roughfac*cquadr(vz2-vz4))
        J[1][3]=(b-d)*exp(roughfac*cquadr(vz2+vz4))


        b=(PHI2[0])[0]*(PSI1[0])[0]+(PSI1[0])[3]
        d=(PHI2[0])[1]*(PSI1[0])[1]+(PHI2[0])[2]*(PSI1[0])[2]
        J[0][1]=(b+d)*exp(roughfac*cquadr(vz1-vz4))
        J[0][3]=(b-d)*exp(roughfac*cquadr(vz1+vz4))

        b=(PHI1[0])[3]*(PSI2[0])[3]+(PSI2[0])[0]
        d=(PHI1[0])[1]*(PSI2[0])[1]+(PHI1[0])[2]*(PSI2[0])[2]
        J[1][0]=(b+d)*exp(roughfac*cquadr(vz2-vz3))
        J[1][2]=(b-d)*exp(roughfac*cquadr(vz2+vz3))

        div=(J[0][1]*J[1][0]-J[1][1]*J[0][0])

        (r[0])[0][0]=(J[1][1]*J[0][2]-J[0][1]*J[1][2])/div # Incoming 1 reflected 1
        (r[0])[0][1]=(J[1][1]*J[0][3]-J[0][1]*J[1][3])/div # Incoming 2 reflected 1
        (r[0])[1][0]=(J[0][0]*J[1][2]-J[1][0]*J[0][2])/div # Incoming 1 reflected 2
        (r[0])[1][1]=(J[0][0]*J[1][3]-J[1][0]*J[0][3])/div # Incoming 2 reflected 2



        (t[0])[0][0]=J[0][2]*(r[0])[0][0]+J[0][3]*(r[0])[1][0]+J[0][0] # Incoming 1 transmitted 1
        (t[0])[0][1]=J[0][2]*(r[0])[0][1]+J[0][3]*(r[0])[1][1]+J[0][1] # Incoming 2 transmitted 1
        (t[0])[1][0]=J[1][2]*(r[0])[0][0]+J[1][3]*(r[0])[1][0]+J[1][0] # Incoming 1 transmitted 2
        (t[0])[1][1]=J[1][2]*(r[0])[0][1]+J[1][3]*(r[0])[1][1]+J[1][1] # Incoming 2 transmitted 2

        (tprime[0])[0][0]=-J[1][1]/div
        (tprime[0])[0][1]=J[0][1]/div
        (tprime[0])[1][0]=J[1][0]/div
        (tprime[0])[1][1]=-J[0][0]/div

        (rprime[0])[0][0]=J[0][2]*((tprime[0])[0][0])+J[0][3]*((tprime[0])[1][0])
        (rprime[0])[0][1]=J[0][2]*((tprime[0])[0][1])+J[0][3]*((tprime[0])[1][1])
        (rprime[0])[1][0]=J[1][2]*((tprime[0])[0][0])+J[1][3]*((tprime[0])[1][0])
        (rprime[0])[1][1]=J[1][2]*((tprime[0])[0][1])+J[1][3]*((tprime[0])[1][1])


    else:

        b=(PSI1[0])[0]
        d=(PSI1[0])[1]*(PHI1[0])[1]

        J[0][0]=(b+d)*exp(roughfac*cquadr(vz1-vz3))
        J[0][2]=(b-d)*exp(roughfac*cquadr(vz1+vz3))


        b=(PSI2[0])[3]
        d=(PHI2[0])[2]*(PSI2[0])[2]
        J[1][1]=(b+d)*exp(roughfac*cquadr(vz2-vz4))
        J[1][3]=(b-d)*exp(roughfac*cquadr(vz2+vz4))


        (r[0])[0][0]=-J[0][2]/J[0][0] # Incoming 1 reflected 1
        (r[0])[0][1]=0 # Incoming 2 reflected 1
        (r[0])[1][0]=0 #/ Incoming 1 reflected 2
        (r[0])[1][1]=-J[1][3]/J[1][1] # Incoming 2 reflected 2


        (t[0])[0][0]=+J[0][2]*(r[0])[0][0]+J[0][0] # Incoming 1 transmitted 1
        (t[0])[0][1]=0 # Incoming 2 transmitted 1
        (t[0])[1][0]=0 # Incoming 1 transmitted 2
        (t[0])[1][1]=+J[1][3]*(r[0])[1][1]+J[1][1] # Incoming 2 transmitted 2

        (tprime[0])[0][0]=+1.0/J[0][0]
        (tprime[0])[0][1]=0
        (tprime[0])[1][0]=0
        (tprime[0])[1][1]=+1.0/J[1][1]


#        (rprime[0])[0][0]=J[0][2]*((tprime[0])[0][0])
#        (rprime[0])[0][1]=0
#        (rprime[0])[1][0]=0
#        (rprime[0])[1][1]=J[1][3]*((tprime[0])[1][1])
        (rprime[0])[0][0]=J[0][2]*((tprime[0])[0][0])
        (rprime[0])[0][1]=0
        (rprime[0])[1][0]=0
        (rprime[0])[1][1]=J[1][3]*((tprime[0])[1][1])

cdef  inline double complex rootfunc(double complex res, double complex D21, double complex D31, double complex eyy):
    return sqrt(  cquadr(res-D21)/4.+cquadr(D31)*eyy  )


cdef void PHI_to_PSI(double complex (*PSI1)[4], double complex (*PHI1)[4], double complex (*PSI2)[4], double complex (*PHI2)[4], int previously_magnetic):
    cdef double complex b,d
    if(previously_magnetic):
        b=2*((PHI1[0])[3]*(PHI2[0])[0]-1)
        d=2*((PHI1[0])[2]*(PHI2[0])[1]-(PHI1[0])[1]*(PHI2[0])[2])

        (PSI1[0])[0]=-1.0/b
        (PSI1[0])[1]=-(PHI2[0])[2]/d
        (PSI1[0])[2]=(PHI2[0])[1]/d
        (PSI1[0])[3]=(PHI2[0])[0]/b

        (PSI2[0])[0]=(PHI1[0])[3]/b
        (PSI2[0])[1]=(PHI1[0])[2]/d
        (PSI2[0])[2]=-(PHI1[0])[1]/d
        (PSI2[0])[3]=-1.0/b
    else:
        (PSI1[0])[0]=0.5
        (PSI1[0])[1]=0.5/(PHI1[0])[1]
        (PSI1[0])[2]=0
        (PSI1[0])[3]=0
        (PSI2[0])[2]=0.5/(PHI2[0])[2]
        (PSI2[0])[3]=0.5
        (PSI2[0])[1]=0
        (PSI2[0])[0]=0


cdef void MagneticPhi(double complex epsxx, double complex epsyy, double complex epszz, double complex epsg, double complex *vz3, double complex *vz4, double complex (*PHI1)[4], double complex (*PHI2)[4], double vy, double vyvy):
    cdef double complex D34, D21, D31, exzexz, b, d, root
    exzexz=cquadr(epsg)
    D34=1.-vyvy/epszz
    D21=epsxx-vyvy+exzexz/epszz
    D31=-vy*epsg/epszz

    b=-exzexz/epszz-epsyy-epsxx+vyvy*(1+epsyy/epszz)
    d=epsyy*epsxx-vyvy*(epsxx*epsyy/epszz+epsyy)+(exzexz+quadr(vyvy))*epsyy/epszz
    root=sqrt(b*b-4*d)

    vz3[0]=sqrt((-b-root)/2)
    vz4[0]=sqrt((-b+root)/2)

    (PHI1[0])[0]=1. #Eigenvectors
    (PHI1[0])[1]=vz3[0]
    (PHI1[0])[3]=-(cquadr(vz3[0])-D21)/D31
    (PHI1[0])[2]=(PHI1[0])[3]*vz3[0]/epsyy

    (PHI2[0])[3]=1.
    (PHI2[0])[2]=vz4[0]/epsyy
    (PHI2[0])[0]=(vz4[0]*(PHI2[0])[2]-D34)/D31
    (PHI2[0])[1]=(PHI2[0])[0]*vz4[0]


cdef void NormalPhi(double complex epsxx, double complex epsyy, double complex epszz, double complex *vz3, double complex *vz4, double complex (*PHI1)[4], double complex (*PHI2)[4], double vyvy):
 #   print "1"
    vz3[0]=CalculateVZsigma(vyvy, epsxx)
    vz4[0]=CalculateVZpi(vyvy, epsyy, epszz)
 #   print "2"
    (PHI1[0])[0]=1. #Eigenvectors
    (PHI1[0])[1]=vz3[0]
    (PHI1[0])[2]=0
    (PHI1[0])[3]=0
    (PHI2[0])[3]=1.
    (PHI2[0])[2]=vz4[0]/epsyy
    (PHI2[0])[0]=0
    (PHI2[0])[1]=0


cdef void Calculate_ANXBN(double complex (*A)[2][2], double complex (*B)[2][2], double complex (*X)[2][2], double N):
    cdef double expite;
    cdef int i,j;
    cdef double complex  resA[2][2];
    cdef double complex  resB[2][2];
    expite=N;


    for i in range(2):
        for j in range(2):
            if(i==j):
                resA[i][j]=1;
                resB[i][j]=1;
            else:
                resA[i][j]=0;
                resB[i][j]=0;




    while(expite):
        if (expite%2):
            Mult2x2_leftside(&resA, A);
            Mult2x2_leftside(&resB, B);

        expite = expite/2
        Mult2x2_leftside(A, A);
        Mult2x2_leftside(B, B);


    for i in range(2):
        for j in range(2):
            (A[0])[i][j]=resA[i][j];
            (B[0])[i][j]=resB[i][j];


    Mult2x2_leftside(X, B)
    Mult2x2_rightside(A, X)




cdef void Calculate_Multilayer_equation(double complex  (*A)[2][2], double complex  (*B)[2][2], double complex (*X)[2][2], double complex  (*result)[2][2], double N):
    # This function calculates efficiently
    # X + AXB + A^2 X B^2 + ... + A^(N-1) X B^(N-1)
    # for complex 2x2-Matrices

    cdef double complex EliB[4];

    cdef double complex Mat[4][4];


    Mat[0][0]=1-(A[0])[0][0]*(B[0])[0][0];
    Mat[0][1]= -(A[0])[0][0]*(B[0])[1][0];
    Mat[0][2]= -(A[0])[0][1]*(B[0])[0][0];
    Mat[0][3]= -(A[0])[0][1]*(B[0])[1][0];

    Mat[1][0]= -(A[0])[0][0]*(B[0])[0][1];
    Mat[1][1]=1-(A[0])[0][0]*(B[0])[1][1];
    Mat[1][2]= -(A[0])[0][1]*(B[0])[0][1];
    Mat[1][3]= -(A[0])[0][1]*(B[0])[1][1];

    Mat[2][0]= -(A[0])[1][0]*(B[0])[0][0];
    Mat[2][1]= -(A[0])[1][0]*(B[0])[1][0];
    Mat[2][2]=1-(A[0])[1][1]*(B[0])[0][0];
    Mat[2][3]= -(A[0])[1][1]*(B[0])[1][0];

    Mat[3][0]= -(A[0])[1][0]*(B[0])[0][1];
    Mat[3][1]= -(A[0])[1][0]*(B[0])[1][1];
    Mat[3][2]= -(A[0])[1][1]*(B[0])[0][1];
    Mat[3][3]=1-(A[0])[1][1]*(B[0])[1][1];



    EliB[0]=(X[0])[0][0];
    EliB[1]=(X[0])[0][1];
    EliB[2]=(X[0])[1][0];
    EliB[3]=(X[0])[1][1];
    Calculate_ANXBN(A, B, X, N);




    EliB[0]-=(X[0])[0][0];
    EliB[1]-=(X[0])[0][1];
    EliB[2]-=(X[0])[1][0];
    EliB[3]-=(X[0])[1][1];
   # printf("EliB: \n");
#    PrintComplex((*X)[0][0]);
#    PrintComplex((*X)[0][1]);
#    PrintComplex((*X)[1][0]);
#    PrintComplex((*X)[1][1]);
    Elimination_4x4(&Mat, &EliB);
    (result[0])[0][0]=EliB[0];
    (result[0])[0][1]=EliB[1];
    (result[0])[1][0]=EliB[2];
    (result[0])[1][1]=EliB[3];




cdef void Paratt_magnetic_y(Heterostructure* HS, double th, double wavelength, double complex (*rtot)[2][2]) except *:


    if((th<=0)|(th>=90)):
        raise Exception("Theta must be in the range 0<theta<90")

    cdef double k0=two_times_pi/wavelength
    cdef double sintheta=sin(deg_to_rad*th)

    cdef double vy=cos(deg_to_rad*th)
    cdef double vyvy=quadr(vy)
    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
    cdef int Lower, Upper
    cdef int ML_is_diagonal=1
    cdef double complex PHI1[4]
    cdef double complex PSI1[4]
    cdef double complex PHI2[4]
    cdef double complex PSI2[4]
    cdef double complex vz1, vz2, vz3, vz4
    cdef double complex D34,D21,D31
    cdef double complex res, root, b, d, exzexz
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


    p[0][1]=0
    p[1][0]=0
   # print "0"
  #  cdef double complex test=LR[MLCOMP[0][0]].ey

    cdef int Cap=NLAYERS-1
    cdef int previously_magnetic=0
   # print "1"
    if(NLAYERS==1):
        Upper=MLCOMP[0][0]
        if(LR[Upper].magdir==2): # magnetic

            MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz1, &vz2, &PHI1, &PHI2, vy, vyvy)

            b=2*(PHI1[3]*PHI2[0]-1)
            d=2*(PHI1[2]*PHI2[1]-PHI1[1]*PHI2[2])

            PSI1[0]=-1.0/b #Inversion of Eigenvectors
            PSI1[1]=-PHI2[2]/d
            PSI1[2]=PHI2[1]/d
            PSI1[3]=PHI2[0]/b

            PSI2[0]=PHI1[3]/b
            PSI2[1]=PHI1[2]/d
            PSI2[2]=-PHI1[1]/d
            PSI2[3]=-1.0/b
            PHI1[0]=1
            PHI1[1]=sintheta
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1
            PHI2[2]=sintheta
            PHI2[1]=0
            PHI2[0]=0
            #Vacuum:
            Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime, 1, vz1, vz2, sintheta,sintheta,LR[MLCOMP[0][0]].Roughness, k0)
        else: # Not magnetic

            NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz1, &vz2, &PHI1, &PHI2, vyvy)


            PSI1[0]=0.5
            PSI1[1]=0.5/vz1
            PSI1[2]=0
            PSI1[3]=0
            PSI2[2]=0.5/PHI2[2]
            PSI2[3]=0.5
            PSI2[1]=0
            PSI2[0]=0
        #Vacuum:
            PHI1[0]=1
            PHI1[1]=sintheta
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1
            PHI2[2]=sintheta
            PHI2[1]=0
            PHI2[0]=0
            Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime, 0, vz1, vz2, sintheta,sintheta,LR[MLCOMP[0][0]].Roughness, k0)
    else:
        Lower=MLCOMP[0][0]
        Upper=MLCOMP[1][0]
        if((LR[Lower].magdir==2)|(LR[Upper].magdir==2)):
            if(LR[Lower].magdir==2): # magnetic

                MagneticPhi(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, LR[Lower].eg, &vz1, &vz2, &PHI1, &PHI2, vy, vyvy)
                previously_magnetic=1
            else:

                NormalPhi(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, &vz1, &vz2, &PHI1, &PHI2, vyvy)
                previously_magnetic=0
            PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)

            if(LR[Upper].magdir==2):

                MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                previously_magnetic=1
            else:
              #  print "here are the phi"
                NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                previously_magnetic=0
           # p=tbtmatrix(exp(2j*k0*LR[1].Thickness*vz3),0,0,exp(2j*k0*LR[1].Thickness*vz4))
            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
#            print "here p is calculated for the normal structure", p[0][0], p[1][1]
#            print "Components:"
#            print k0, Upper, LR[1].Thickness, vz3, vz4
            Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime,1 , vz1, vz2, vz3,vz4,LR[MLCOMP[0][0]].Roughness, k0)
        #    print "here p is calculated for the normal structure", p[0][0], p[1][1]
       #     print "Components:"
         #   print k0, Upper, LR[Upper].Thickness, vz3, vz4

        else:

          #  print "lowest interface is not magnetic"
            vz1=CalculateVZsigma(vyvy, LR[MLCOMP[0][0]].ex)
            vz2=CalculateVZpi(vyvy, LR[MLCOMP[0][0]].ey, LR[MLCOMP[0][0]].ez)
            vz3=CalculateVZsigma(vyvy, LR[MLCOMP[1][0]].ex)
            vz4=CalculateVZpi(vyvy, LR[MLCOMP[1][0]].ey, LR[MLCOMP[1][0]].ez)
            PHI1[0]=1.
            PHI1[1]=vz3
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1.
            PHI2[2]=vz4/LR[MLCOMP[1][0]].ey
            PHI2[0]=0
            PHI2[1]=0

            p[0][0]=exp(1j*k0*LR[MLCOMP[1][0]].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[MLCOMP[1][0]].Thickness*vz4)
         #   print "here p is calculated for the normal structure", p[0][0], p[1][1]
         #   print "Components:"
         #   print k0, "1", LR[1].Thickness, vz3, vz4

            (rtot[0])[0][0]=((vz3-vz1)/(vz3+vz1))*exp(-2.0*quadr(k0)*vz1*vz3*quadr(LR[MLCOMP[0][0]].Roughness))
            (rtot[0])[1][0]=0
            (rtot[0])[0][1]=0
            (rtot[0])[1][1]=((vz4*LR[MLCOMP[0][0]].ey-vz2*LR[MLCOMP[1][0]].ey)/(vz4*LR[MLCOMP[0][0]].ey+vz2*LR[MLCOMP[1][0]].ey))*exp(-2.0*quadr(k0)*vz2*vz4*quadr(LR[MLCOMP[0][0]].Roughness))



    i=1
    while i<NLAYERS:
       # print "loop start", i
        if(MLLENGTH[i]==1):

            vz1=vz3
            vz2=vz4
             #Inversion of Eigenvectors
            PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)

            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
               # print "i, Upper", i, Upper
                if(LR[Upper].magdir==2):
                  #  print "hallo"

                    MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                 #   print "x"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)
                 #   print "y"
                    if(previously_magnetic):
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    previously_magnetic=0
            else:
              #  print "vacuum"
                PHI1[0]=1. #Eigenvectors
                PHI1[1]=sintheta
                PHI1[2]=0
                PHI1[3]=0
                PHI2[3]=1.
                PHI2[2]=sintheta
                PHI2[0]=0
                PHI2[1]=0
                if(previously_magnetic):
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][0]].Roughness, k0)
                else:
                  #  print "else"
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][0]].Roughness, k0)


            Mult2x2_rightside(&p, rtot)

            Mult2x2_leftside(rtot, &p)

            Mult2x2_rightside(&tprime, rtot)

            Mult2x2_leftside(rtot, &t)

            (rtot[0])[0][0]+=r[0][0]
            (rtot[0])[1][0]+=r[1][0]
            (rtot[0])[0][1]+=r[0][1]
            (rtot[0])[1][1]+=r[1][1]

#                Mult2x2_rightside(&p, &r_ML_in_1) #p_B r p_B and so on
#                Mult2x2_leftside(&r_ML_in_1, &p)
#                Mult2x2_rightside(&tprime, &r_ML_in_1) #t' p_B r p_b t and so on
#                Mult2x2_leftside(&r_ML_in_1, &t)
#                r_ML_in_1[0][0]+=r[0][0]
#                r_ML_in_1[1][0]+=r[1][0]
#                r_ML_in_1[0][1]+=r[0][1]
#                r_ML_in_1[1][1]+=r[1][1]

#            if(i==1):
#                print "new rtot"
#                print (rtot[0])[0][0]
#                print (rtot[0])[0][1]
#                print (rtot[0])[1][0]
#                print (rtot[0])[1][1]


            if(i!=Cap):
                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)




        else: #A Multilayer
          #  print "1"
            PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
          #  print "2"
            Upper=MLCOMP[i][1]
            vz1=vz3
            vz2=vz4
            if(LR[Upper].magdir==2):
                ML_is_diagonal=0
              #  print "x"
                MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                previously_magnetic=1

            else:
           #     print "Hallo"
                NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                if(previously_magnetic):
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                else:
               #     print "else"
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                previously_magnetic=0
         #   print "3"
#            print "t BA"
#            print t_ML_in_1[0][0]
#            print t_ML_in_1[0][1]
#            print t_ML_in_1[1][0]
#            print t_ML_in_1[1][1]
#
#            print "p A"
#
#            print p[0][0]
#            print p[1][1]
            Mult2x2_leftside(&t_ML_back_1, &p) # t'(AB)*p(A)
            Mult2x2_rightside(&p, &t_ML_in_1) # p(A) * t(AB)
#            print "r_ML_in_1 ", 0
#            print r_ML_in_1[0][0]
#            print r_ML_in_1[0][1]
#            print r_ML_in_1[1][0]
#            print r_ML_in_1[1][1]


            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)

            j=1
            while j<MLLENGTH[i]:
                Upper=MLCOMP[i][(j+1)%MLLENGTH[i]]
                vz1=vz3
                vz2=vz4
                PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "5"
              #  Upper=MLCOMP[i][1]
                if(LR[Upper].magdir==2):
                  #  print "hallo"
                    MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)
                    ML_is_diagonal=0
                #    print "y"
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    previously_magnetic=0
#                if(j==1):
#                    print "t CB"
#                elif(j==2):
#                    print "t AC"
#                print t[0][0]
#                print t[0][1]
#                print t[1][0]
#                print t[1][1]
                Mult2x2_rightside(&p, &r_ML_in_1) #p_B r p_B and so on

                Mult2x2_leftside(&r_ML_in_1, &p)

                Mult2x2_rightside(&tprime, &r_ML_in_1) #t' p_B r p_b t and so on

                Mult2x2_leftside(&r_ML_in_1, &t)


                r_ML_in_1[0][0]+=r[0][0]
                r_ML_in_1[1][0]+=r[1][0]
                r_ML_in_1[0][1]+=r[0][1]
                r_ML_in_1[1][1]+=r[1][1]
#                print "r_ML_in_1 ", j
#                print r_ML_in_1[0][0]
#                print r_ML_in_1[0][1]
#                print r_ML_in_1[1][0]
#                print r_ML_in_1[1][1]
#                print "new rml"
#                print r_ML_in_1[0][0]
#                print r_ML_in_1[0][1]
#                print r_ML_in_1[1][0]
#                print r_ML_in_1[1][1]

#                if(j==1):
#                    print "p B"
#                elif(j==2):
#                    print "p C"
              #  print p[0][0]
            #    print p[1][1]
                Mult2x2_leftside(&t_ML_in_1, &p)
                Mult2x2_rightside(&p, &t_ML_back_1)

                Mult2x2_leftside(&t_ML_in_1, &t) # p(A) * t(AB) p_B t(BC) and so on
                Mult2x2_rightside(&tprime, &t_ML_back_1)# t'(BC) p_B t'(AB)*p(A) and so on


                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
           #     if(j==1):
            #        print "p C components"
            #        print Upper, LR[Upper].Thickness, vz3, vz4

                j=j+1



            if(ML_is_diagonal):
                r_ML_in_2[0][0]=r_ML_in_1[0][0]*(1-ipow(t_ML_in_1[0][0]*t_ML_back_1[0][0], MLREP[i]-1))/(1-t_ML_in_1[0][0]*t_ML_back_1[0][0])
                r_ML_in_2[1][1]=r_ML_in_1[1][1]*(1-ipow(t_ML_in_1[1][1]*t_ML_back_1[1][1], MLREP[i]-1))/(1-t_ML_in_1[1][1]*t_ML_back_1[1][1])
#                print "diagonal"
#                print t_ML_in_1[0][0]
#                print t_ML_in_1[1][1]
#                print t_ML_back_1[0][0]
#                print t_ML_back_1[1][1]
                t_ML_in_1[0][0]=ipow(t_ML_in_1[0][0],  MLREP[i]-1)
                t_ML_in_1[1][1]=ipow(t_ML_in_1[1][1],  MLREP[i]-1)
                t_ML_back_1[0][0]=ipow(t_ML_back_1[0][0],  MLREP[i]-1)
                t_ML_back_1[1][1]=ipow(t_ML_back_1[1][1],  MLREP[i]-1)

            else:

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
                vz1=vz3
                vz2=vz4
                Upper=MLCOMP[i][j]
                PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "10"
                if(LR[Upper].magdir==2):
                  #  print "hallo"
                    MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)


                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    previously_magnetic=0
              #  print "11"
                Mult2x2_rightside(&p, rtot)
                Mult2x2_leftside(rtot, &p)
                Mult2x2_rightside(&tprime, rtot)
                Mult2x2_leftside(rtot, &t)
                (rtot[0])[0][0]+=r[0][0]
                (rtot[0])[1][0]+=r[1][0]
                (rtot[0])[0][1]+=r[0][1]
                (rtot[0])[1][1]+=r[1][1]


                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
                j=j+1
            vz1=vz3
            vz2=vz4
            if(i==(NLAYERS-1)):
                PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "13"
              #  vz3=CalculateVZsigma(vyvy, LR[Upper].ex)
              #  vz4=CalculateVZpi(vyvy, LR[Upper].ey, LR[Upper].ez)
                PHI1[0]=1. #Eigenvectors
                PHI1[1]=sintheta
                PHI1[2]=0
                PHI1[3]=0
                PHI2[3]=1.
                PHI2[2]=sintheta
                PHI2[0]=0
                PHI2[1]=0
                if(previously_magnetic):
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                else:
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
              #  print "14"
            else:
                Upper=MLCOMP[i+1][0]
                PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
            #    print "15"
                if(LR[Upper].magdir==2):
                  #  print "hallo"
                    MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    previously_magnetic=0
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

            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
          #  print "18"
        i=i+1


cdef void Invert2x2(double complex (*M)[2][2]):
    cdef double complex dconj=cconj( (M[0])[0][0]*(M[0])[1][1]-(M[0])[0][1]*(M[0])[1][0] )
    cdef double complex safe=(M[0])[0][0]
    dconj=dconj/cabsquadr(dconj)
    (M[0])[0][0]=(M[0])[1][1]*dconj
    (M[0])[0][1]=-(M[0])[0][1]*dconj
    (M[0])[1][0]=-(M[0])[1][0]*dconj
    (M[0])[1][1]=safe*dconj

cdef void FillC0(double complex (*C0)[2][2], double complex  (*rprime)[2][2], double complex (*rtot)[2][2],double complex  (*p)[2][2]):

    (C0[0])[0][0]=(p[0])[0][0]
    (C0[0])[0][1]=0
    (C0[0])[1][0]=0
    (C0[0])[1][1]=(p[0])[1][1]

    Mult2x2_rightside(rtot, C0)
    Mult2x2_rightside(p,C0)
    Mult2x2_rightside(rprime,C0)
    (C0[0])[0][0]=1.-(C0[0])[0][0]
    (C0[0])[0][1]=-(C0[0])[0][1]
    (C0[0])[1][0]=-(C0[0])[1][0]
    (C0[0])[1][1]=1.-(C0[0])[1][1]
    Invert2x2(C0)


cdef void Calculate_Multilayer_with_Matrices(double complex (*t_comp1_up)[2][2], double complex (*t_comp2_up)[2][2], double complex (*t_comp1_do)[2][2], double complex (*t_comp2_do)[2][2], double complex (*r_ML_in1)[2][2], double complex (*r_ML_in2)[2][2], double complex (*r_ML_ba1)[2][2], double complex (*r_ML_ba2)[2][2], double N):

    cdef double complex MLfac1[2][2]
    cdef double complex MLfac2[2][2]
    cdef double complex safe[2][2]

    if(N==0):
        (r_ML_in2[0])[0][0]=0
        (r_ML_in2[0])[0][1]=0
        (r_ML_in2[0])[1][0]=0
        (r_ML_in2[0])[1][1]=0

        (r_ML_ba2[0])[0][0]=0
        (r_ML_ba2[0])[0][1]=0
        (r_ML_ba2[0])[1][0]=0
        (r_ML_ba2[0])[1][1]=0

        (t_comp2_up[0])[0][0]=1
        (t_comp2_up[0])[0][1]=0
        (t_comp2_up[0])[1][0]=0
        (t_comp2_up[0])[1][1]=1

        (t_comp2_do[0])[0][0]=1
        (t_comp2_do[0])[0][1]=0
        (t_comp2_do[0])[1][0]=0
        (t_comp2_do[0])[1][1]=1
        return

    (r_ML_in2[0])[0][0]=(r_ML_in1[0])[0][0]
    (r_ML_in2[0])[0][1]=(r_ML_in1[0])[0][1]
    (r_ML_in2[0])[1][0]=(r_ML_in1[0])[1][0]
    (r_ML_in2[0])[1][1]=(r_ML_in1[0])[1][1]

    (r_ML_ba2[0])[0][0]=(r_ML_ba1[0])[0][0]
    (r_ML_ba2[0])[0][1]=(r_ML_ba1[0])[0][1]
    (r_ML_ba2[0])[1][0]=(r_ML_ba1[0])[1][0]
    (r_ML_ba2[0])[1][1]=(r_ML_ba1[0])[1][1]

    (t_comp2_up[0])[0][0]=(t_comp1_up[0])[0][0]
    (t_comp2_up[0])[0][1]=(t_comp1_up[0])[0][1]
    (t_comp2_up[0])[1][0]=(t_comp1_up[0])[1][0]
    (t_comp2_up[0])[1][1]=(t_comp1_up[0])[1][1]

    (t_comp2_do[0])[0][0]=(t_comp1_do[0])[0][0]
    (t_comp2_do[0])[0][1]=(t_comp1_do[0])[0][1]
    (t_comp2_do[0])[1][0]=(t_comp1_do[0])[1][0]
    (t_comp2_do[0])[1][1]=(t_comp1_do[0])[1][1]

    N=N-1
  #  print N
    while(N):
        if(N%2): #Adding one Compound Layer to the heterostructure
           # MLfac=1.0/(1-rres1*r_ML_ba1[0])
            MLfac1[0][0]=(r_ML_ba1[0])[0][0]
            MLfac1[0][1]=(r_ML_ba1[0])[0][1]
            MLfac1[1][0]=(r_ML_ba1[0])[1][0]
            MLfac1[1][1]=(r_ML_ba1[0])[1][1]
            Mult2x2_rightside(r_ML_in2, &MLfac1)
            MLfac1[0][0]=1.-MLfac1[0][0]
            MLfac1[0][1]=-MLfac1[0][1]
            MLfac1[1][0]=-MLfac1[1][0]
            MLfac1[1][1]=1.-MLfac1[1][1]
            Invert2x2(&MLfac1)

            MLfac2[0][0]=(r_ML_ba1[0])[0][0]
            MLfac2[0][1]=(r_ML_ba1[0])[0][1]
            MLfac2[1][0]=(r_ML_ba1[0])[1][0]
            MLfac2[1][1]=(r_ML_ba1[0])[1][1]
            Mult2x2_leftside(&MLfac2, r_ML_in2)
            MLfac2[0][0]=1.-MLfac2[0][0]
            MLfac2[0][1]=-MLfac2[0][1]
            MLfac2[1][0]=-MLfac2[1][0]
            MLfac2[1][1]=1.-MLfac2[1][1]
            Invert2x2(&MLfac2)

           # rres1=r_ML_in1[0]+t_comp1_up[0]*t_comp1_do[0]*rres1*MLfac

            safe[0][0]=(t_comp1_do[0])[0][0]
            safe[0][1]=(t_comp1_do[0])[0][1]
            safe[1][0]=(t_comp1_do[0])[1][0]
            safe[1][1]=(t_comp1_do[0])[1][1]

            Mult2x2_rightside(&MLfac2, &safe)
            Mult2x2_rightside(r_ML_in2, &safe)
            Mult2x2_rightside(t_comp1_up, &safe)

            (r_ML_in2[0])[0][0]=(r_ML_in1[0])[0][0]+safe[0][0]
            (r_ML_in2[0])[0][1]=(r_ML_in1[0])[0][1]+safe[0][1]
            (r_ML_in2[0])[1][0]=(r_ML_in1[0])[1][0]+safe[1][0]
            (r_ML_in2[0])[1][1]=(r_ML_in1[0])[1][1]+safe[1][1]

          #  rres2=rres2+tres_up*tres_do*r_ML_ba1[0]*MLfac

            safe[0][0]=(t_comp2_up[0])[0][0]
            safe[0][1]=(t_comp2_up[0])[0][1]
            safe[1][0]=(t_comp2_up[0])[1][0]
            safe[1][1]=(t_comp2_up[0])[1][1]

            Mult2x2_rightside(&MLfac1, &safe)
            Mult2x2_rightside(r_ML_ba1, &safe)
            Mult2x2_rightside(t_comp2_do, &safe)

            (r_ML_ba2[0])[0][0]=(r_ML_ba2[0])[0][0]+safe[0][0]
            (r_ML_ba2[0])[0][1]=(r_ML_ba2[0])[0][1]+safe[0][1]
            (r_ML_ba2[0])[1][0]=(r_ML_ba2[0])[1][0]+safe[1][0]
            (r_ML_ba2[0])[1][1]=(r_ML_ba2[0])[1][1]+safe[1][1]

          #  tres_up=t_comp1_up[0]*tres_up*MLfac

            Mult2x2_rightside(&MLfac1, t_comp2_up)
            Mult2x2_rightside(t_comp1_up, t_comp2_up)


          #  tres_do=t_comp1_do[0]*tres_do*MLfac

            Mult2x2_leftside(t_comp2_do, &MLfac2)
            Mult2x2_leftside(t_comp2_do, t_comp1_do)
        #endif
        N=N/2


        #Doubling the heterostructure:
     #   MLfac=1.0/(1-r_ML_in1[0]*r_ML_ba1[0])

        MLfac1[0][0]=(r_ML_ba1[0])[0][0]
        MLfac1[0][1]=(r_ML_ba1[0])[0][1]
        MLfac1[1][0]=(r_ML_ba1[0])[1][0]
        MLfac1[1][1]=(r_ML_ba1[0])[1][1]
        Mult2x2_rightside(r_ML_in1, &MLfac1)
        MLfac1[0][0]=1.-MLfac1[0][0]
        MLfac1[0][1]=-MLfac1[0][1]
        MLfac1[1][0]=-MLfac1[1][0]
        MLfac1[1][1]=1.-MLfac1[1][1]
        Invert2x2(&MLfac1)

        MLfac2[0][0]=(r_ML_ba1[0])[0][0]
        MLfac2[0][1]=(r_ML_ba1[0])[0][1]
        MLfac2[1][0]=(r_ML_ba1[0])[1][0]
        MLfac2[1][1]=(r_ML_ba1[0])[1][1]
        Mult2x2_leftside(&MLfac2, r_ML_in1)
        MLfac2[0][0]=1.-MLfac2[0][0]
        MLfac2[0][1]=-MLfac2[0][1]
        MLfac2[1][0]=-MLfac2[1][0]
        MLfac2[1][1]=1.-MLfac2[1][1]
        Invert2x2(&MLfac2)

     #   r_ML_in1[0]=r_ML_in1[0]+t_comp1_up[0]*t_comp1_do[0]*r_ML_in1[0]*MLfac

        safe[0][0]=(t_comp1_do[0])[0][0]
        safe[0][1]=(t_comp1_do[0])[0][1]
        safe[1][0]=(t_comp1_do[0])[1][0]
        safe[1][1]=(t_comp1_do[0])[1][1]
        Mult2x2_rightside(&MLfac2, &safe)
        Mult2x2_rightside(r_ML_in1, &safe)
        Mult2x2_rightside(t_comp1_up, &safe)
        (r_ML_in1[0])[0][0]=(r_ML_in1[0])[0][0]+safe[0][0]
        (r_ML_in1[0])[0][1]=(r_ML_in1[0])[0][1]+safe[0][1]
        (r_ML_in1[0])[1][0]=(r_ML_in1[0])[1][0]+safe[1][0]
        (r_ML_in1[0])[1][1]=(r_ML_in1[0])[1][1]+safe[1][1]


     #   r_ML_ba1[0]=r_ML_ba1[0]+t_comp1_up[0]*t_comp1_do[0]*r_ML_ba1[0]*MLfac

        safe[0][0]=(t_comp1_up[0])[0][0]
        safe[0][1]=(t_comp1_up[0])[0][1]
        safe[1][0]=(t_comp1_up[0])[1][0]
        safe[1][1]=(t_comp1_up[0])[1][1]
        Mult2x2_rightside(&MLfac1, &safe)
        Mult2x2_rightside(r_ML_ba1, &safe)
        Mult2x2_rightside(t_comp1_do, &safe)
        (r_ML_ba1[0])[0][0]=(r_ML_ba1[0])[0][0]+safe[0][0]
        (r_ML_ba1[0])[0][1]=(r_ML_ba1[0])[0][1]+safe[0][1]
        (r_ML_ba1[0])[1][0]=(r_ML_ba1[0])[1][0]+safe[1][0]
        (r_ML_ba1[0])[1][1]=(r_ML_ba1[0])[1][1]+safe[1][1]
    #    t_comp1_up[0]=cquadr(t_comp1_up[0])*MLfac

        safe[0][0]=(t_comp1_up[0])[0][0]
        safe[0][1]=(t_comp1_up[0])[0][1]
        safe[1][0]=(t_comp1_up[0])[1][0]
        safe[1][1]=(t_comp1_up[0])[1][1]

        Mult2x2_rightside(&MLfac1, &safe)
        Mult2x2_leftside(t_comp1_up, &safe)

        safe[0][0]=(t_comp1_do[0])[0][0]
        safe[0][1]=(t_comp1_do[0])[0][1]
        safe[1][0]=(t_comp1_do[0])[1][0]
        safe[1][1]=(t_comp1_do[0])[1][1]
     #   t_comp1_do[0]=cquadr(t_comp1_do[0])*MLfac

        Mult2x2_rightside(&MLfac2, &safe)
        Mult2x2_leftside(t_comp1_do, &safe)



cdef void Paratt_magnetic_y_MS(Heterostructure* HS, double th, double wavelength, double complex (*rtot)[2][2]) except *:


    if((th<=0)|(th>=90)):
        raise Exception("Theta must be in the range 0<theta<90")

    cdef double k0=two_times_pi/wavelength
    cdef double sintheta=sin(deg_to_rad*th)

    cdef double vy=cos(deg_to_rad*th)
    cdef double vyvy=quadr(vy)
    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
    cdef int Lower, Upper
    cdef int ML_is_diagonal=1
    cdef double complex PHI1[4]
    cdef double complex PSI1[4]
    cdef double complex PHI2[4]
    cdef double complex PSI2[4]
    cdef double complex vz1, vz2, vz3, vz4
    cdef double complex D34,D21,D31
    cdef double complex res, root, b, d, exzexz
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

    p[0][1]=0
    p[1][0]=0
   # print "0"
   # cdef double complex test=LR[MLCOMP[0][0]].ey

    cdef int Cap=NLAYERS-1
    cdef int previously_magnetic=0
   # print "1"
    if(NLAYERS==1):
        Upper=MLCOMP[0][0]
        if(LR[Upper].magdir==2): # magnetic
          #  print "x"
            MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

            b=2*(PHI1[3]*PHI2[0]-1)
            d=2*(PHI1[2]*PHI2[1]-PHI1[1]*PHI2[2])

            PSI1[0]=-1.0/b #Inversion of Eigenvectors
            PSI1[1]=-PHI2[2]/d
            PSI1[2]=PHI2[1]/d
            PSI1[3]=PHI2[0]/b

            PSI2[0]=PHI1[3]/b
            PSI2[1]=PHI1[2]/d
            PSI2[2]=-PHI1[1]/d
            PSI2[3]=-1.0/b
            PHI1[0]=1
            PHI1[1]=sintheta
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1
            PHI2[2]=sintheta
            PHI2[1]=0
            PHI2[0]=0
            #Vacuum:
            Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime, 1, vz3, vz4, sintheta,sintheta,LR[MLCOMP[0][0]].Roughness, k0)
        else: # Not magnetic
       #     print "y"
            NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

            PSI1[0]=0.5
            PSI1[1]=0.5/vz3
            PSI1[2]=0
            PSI1[3]=0
            PSI2[2]=0.5/PHI2[2]
            PSI2[3]=0.5
            PSI2[1]=0
            PSI2[0]=0
        #Vacuum:
            PHI1[0]=1
            PHI1[1]=sintheta
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1
            PHI2[2]=sintheta
            PHI2[1]=0
            PHI2[0]=0
            Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime, 0, vz3, vz4, sintheta,sintheta,LR[MLCOMP[0][0]].Roughness, k0)
    else:
        Lower=MLCOMP[0][0]
        Upper=MLCOMP[1][0]
        if((LR[Lower].magdir==2)|(LR[Upper].magdir==2)):
            if(LR[Lower].magdir==2): # magnetic

                MagneticPhi(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, LR[Lower].eg, &vz1, &vz2, &PHI1, &PHI2, vy, vyvy)
                previously_magnetic=1
            else:

                NormalPhi(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, &vz1, &vz2, &PHI1, &PHI2, vyvy)
                previously_magnetic=0
            PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)

            if(LR[Upper].magdir==2):

                MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                previously_magnetic=1
            else:
              #  print "here are the phi"
                NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                previously_magnetic=0
           # p=tbtmatrix(exp(2j*k0*LR[1].Thickness*vz3),0,0,exp(2j*k0*LR[1].Thickness*vz4))
            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
#            print "here p is calculated for the normal structure", p[0][0], p[1][1]
#            print "Components:"
#            print k0, Upper, LR[1].Thickness, vz3, vz4
            Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime, 1 , vz1, vz2, vz3,vz4,LR[MLCOMP[0][0]].Roughness, k0)
        #    print "here p is calculated for the normal structure", p[0][0], p[1][1]
       #     print "Components:"
         #   print k0, Upper, LR[Upper].Thickness, vz3, vz4

        else:

          #  print "lowest interface is not magnetic"
            vz1=CalculateVZsigma(vyvy, LR[MLCOMP[0][0]].ex)
            vz2=CalculateVZpi(vyvy, LR[MLCOMP[0][0]].ey, LR[MLCOMP[0][0]].ez)
            vz3=CalculateVZsigma(vyvy, LR[MLCOMP[1][0]].ex)
            vz4=CalculateVZpi(vyvy, LR[MLCOMP[1][0]].ey, LR[MLCOMP[1][0]].ez)
            PHI1[0]=1.
            PHI1[1]=vz3
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1.
            PHI2[2]=vz4/LR[MLCOMP[1][0]].ey
            PHI2[0]=0
            PHI2[1]=0

            p[0][0]=exp(1j*k0*LR[MLCOMP[1][0]].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[MLCOMP[1][0]].Thickness*vz4)

            (rtot[0])[0][0]=((vz3-vz1)/(vz3+vz1))*exp(-2.0*quadr(k0)*vz1*vz3*quadr(LR[MLCOMP[0][0]].Roughness))
            (rtot[0])[1][0]=0
            (rtot[0])[0][1]=0
            (rtot[0])[1][1]=((vz4*LR[MLCOMP[0][0]].ey-vz2*LR[MLCOMP[1][0]].ey)/(vz4*LR[MLCOMP[0][0]].ey+vz2*LR[MLCOMP[1][0]].ey))*exp(-2.0*quadr(k0)*vz2*vz4*quadr(LR[MLCOMP[0][0]].Roughness))


    i=1
    while i<NLAYERS:
       # print "loop start", i
        if(MLLENGTH[i]==1):

            vz1=vz3
            vz2=vz4
             #Inversion of Eigenvectors
            PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)

            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
               # print "i, Upper", i, Upper
                if(LR[Upper].magdir==2):
                  #  print "hallo"

                    MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                 #   print "x"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)
                 #   print "y"
                    if(previously_magnetic):
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    previously_magnetic=0
            else:
              #  print "vacuum"
                PHI1[0]=1. #Eigenvectors
                PHI1[1]=sintheta
                PHI1[2]=0
                PHI1[3]=0
                PHI2[3]=1.
                PHI2[2]=sintheta
                PHI2[0]=0
                PHI2[1]=0
                if(previously_magnetic):
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][0]].Roughness, k0)
                else:
                  #  print "else"
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][0]].Roughness, k0)

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
                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
        else:
            PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
          #  print "2"
            Upper=MLCOMP[i][1]
            vz1=vz3
            vz2=vz4

            if(LR[Upper].magdir==2):
                ML_is_diagonal=0
              #  print "x"
                MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)

                previously_magnetic=1

            else:
           #     print "Hallo"
                NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                if(previously_magnetic):
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                else:
               #     print "else"
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)

                previously_magnetic=0


            Mult2x2_leftside(&t_ML_back_1, &p) # t'(AB)*p(A)


            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)




            j=1
            while j<MLLENGTH[i]:
                Upper=MLCOMP[i][(j+1)%MLLENGTH[i]]
              #  print "j and Upper", j, (j+1)%MLLENGTH[i], Upper
                vz1=vz3
                vz2=vz4
                PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "5"
              #  Upper=MLCOMP[i][1]

                if(LR[Upper].magdir==2):
                  #  print "hallo"
                    MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)
                    ML_is_diagonal=0
                #    print "y"
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    previously_magnetic=0
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
                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
                j=j+1


            p[0][0]=exp(1j*k0*LR[MLCOMP[i][MLLENGTH[i]-1]].Thickness*vz1)
            p[1][1]=exp(1j*k0*LR[MLCOMP[i][MLLENGTH[i]-1]].Thickness*vz2)
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
                Lower=MLCOMP[i][j-1]
                Upper=MLCOMP[i][j]
                vz3=vz1
                vz4=vz2

                if(LR[Lower].magdir==2):
                    MagneticPhi(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, LR[Lower].eg, &vz1, &vz2, &PHI1, &PHI2, vy, vyvy)
                    ML_is_diagonal=0
                    previously_magnetic=1
                else:
                    NormalPhi(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, &vz1, &vz2, &PHI1, &PHI2, vyvy)
                    previously_magnetic=0


                PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "5"
              #  Upper=MLCOMP[i][1]
                if(LR[Upper].magdir==2):
                  #  print "hallo"
                    MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)
                    ML_is_diagonal=0
                #    print "y"
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[Lower].Roughness, k0)
                  #  previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[Lower].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[Lower].Roughness, k0)



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
                p[0][0]=exp(1j*k0*LR[Lower].Thickness*vz1)
                p[1][1]=exp(1j*k0*LR[Lower].Thickness*vz2)

                Mult2x2_rightside(&p, &t_ML_in_1)

                j=j-1


            Mult2x2_rightside(&p, &r_ML_back_1)
            Mult2x2_leftside(&r_ML_back_1, &p)


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


            Mult2x2_rightside(&t_ML_back_2, rtot) #Here the ML part is missing
            Mult2x2_leftside(rtot, &C0)
            Mult2x2_leftside(rtot, &t_ML_in_2) # (t'(CA) p_C t'(BC) p_B t'(AB)*p(A))^N rtot (p(A) * t(AB) p_B t(BC) p_C t(CA))^N
            (rtot[0])[0][0]+=r_ML_in_2[0][0]
            (rtot[0])[1][0]+=r_ML_in_2[1][0]
            (rtot[0])[0][1]+=r_ML_in_2[0][1]
            (rtot[0])[1][1]+=r_ML_in_2[1][1]


            Upper=MLCOMP[i][0]
            if(LR[Upper].magdir==2):
              #  print "hallo"
                MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

         #       Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                previously_magnetic=1
            else:
           #     print "Hallo"
                NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                previously_magnetic=0
            j=1

            while j<MLLENGTH[i]:
                Upper=MLCOMP[i][j]
                vz1=vz3
                vz2=vz4

                PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "10"
                if(LR[Upper].magdir==2):
                  #  print "hallo"
                    MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)


                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    previously_magnetic=0
           #     print "4", vz1, vz2, vz3, vz4
              #  print "11"
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

                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
                j=j+1
              #  print "12"
            if(i==(NLAYERS-1)):
                vz1=vz3
                vz2=vz4
                PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "13"
              #  vz3=CalculateVZsigma(vyvy, LR[Upper].ex)
              #  vz4=CalculateVZpi(vyvy, LR[Upper].ey, LR[Upper].ez)
                PHI1[0]=1. #Eigenvectors
                PHI1[1]=sintheta
                PHI1[2]=0
                PHI1[3]=0
                PHI2[3]=1.
                PHI2[2]=sintheta
                PHI2[0]=0
                PHI2[1]=0
                if(previously_magnetic):
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                else:
                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
              #  print "14"

            else:
                Upper=MLCOMP[i+1][0]
                vz1=vz3
                vz2=vz4
                PHI_to_PSI(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
            #    print "15"
                if(LR[Upper].magdir==2):
                  #  print "hallo"
                    MagneticPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                    Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    previously_magnetic=0
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


            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)

        i=i+1


cdef  inline double complex rootfunc2(double complex D21, double complex b, double complex exyexy, double complex D34):
    return sqrt( cquadr(0.5*(D21-b))-exyexy*D34 )

cdef void MagneticPhi_z(double complex epsxx, double complex epsyy, double complex epszz, double complex epsg, double complex *vz3, double complex *vz4, double complex (*PHI1)[4], double complex (*PHI2)[4], double vy, double vyvy):
    cdef double complex D34, D21, D23, exyexy, b, root

    D34=1.-vyvy/epszz
    D21=epsxx-vyvy
    D23=-epsg
    exyexy=cquadr(epsg)

    b=D34*epsyy
    root=rootfunc2(D21, b, exyexy, D34)


    vz3[0]=sqrt( 0.5*(D21+b)+root)
    vz4[0]=sqrt( 0.5*(D21+b)-root)

    (PHI1[0])[0]=1. #Eigenvectors
    (PHI1[0])[1]=vz3[0]
    (PHI1[0])[2]=(cquadr(vz3[0])-D21)/D23

    (PHI1[0])[3]=(PHI1[0])[2]*(vz3[0]/D34)

    (PHI2[0])[3]=1.
    (PHI2[0])[2]=D34/vz4[0]
    (PHI2[0])[0]=-(vz4[0]-epsyy*(PHI2[0])[2])/D23
    (PHI2[0])[1]=vz4[0]*(PHI2[0])[0]


cdef void PHI_to_PSI_z(double complex (*PSI1)[4], double complex (*PHI1)[4], double complex (*PSI2)[4], double complex (*PHI2)[4], int previously_magnetic):
    cdef double complex b,d
    if(previously_magnetic):
        b=2*((PHI2[0])[2]-(PHI2[0])[0]*(PHI1[0])[2])
        d=2*((PHI1[0])[1]-(PHI2[0])[1]*(PHI1[0])[3])

        (PSI1[0])[0]=(PHI2[0])[2]/b
        (PSI1[0])[1]=1.0/d
        (PSI1[0])[2]=-(PHI2[0])[0]/b
        (PSI1[0])[3]=-(PHI2[0])[1]/d

        (PSI2[0])[0]=-(PHI1[0])[2]/b
        (PSI2[0])[1]=-(PHI1[0])[3]/d
        (PSI2[0])[2]=1.0/b
        (PSI2[0])[3]=(PHI1[0])[1]/d
    else:
        (PSI1[0])[0]=0.5
        (PSI1[0])[1]=0.5/(PHI1[0])[1]
        (PSI1[0])[2]=0
        (PSI1[0])[3]=0
        (PSI2[0])[2]=0.5/(PHI2[0])[2]
        (PSI2[0])[3]=0.5
        (PSI2[0])[1]=0
        (PSI2[0])[0]=0

cdef void Calculate_rt_z(double complex (*PSI1)[4], double complex (*PHI1)[4], double complex (*PSI2)[4], double complex (*PHI2)[4], double complex (*r)[2][2], double complex (*rprime)[2][2], double complex (*t)[2][2], double complex (*tprime)[2][2], int Magnetic, \
                       double complex vz1, double complex vz2, double complex vz3, double complex vz4, double sigma, double k0):

    cdef double complex J[2][4]
    cdef double complex b,d, div;
    cdef double roughfac=-0.5*quadr(sigma)*quadr(k0)
    if(Magnetic):




        b=(PHI1[0])[2]*(PSI1[0])[2]+(PSI1[0])[0]
        d=(PHI1[0])[3]*(PSI1[0])[3]+(PSI1[0])[1]*(PHI1[0])[1]
        J[0][0]=(b+d)*exp(roughfac*cquadr(vz1-vz3))
        J[0][2]=(b-d)*exp(roughfac*cquadr(vz1+vz3))



        b=(PHI2[0])[1]*(PSI2[0])[1]+(PSI2[0])[3]
        d=(PHI2[0])[2]*(PSI2[0])[2]+(PSI2[0])[0]*(PHI2[0])[0]
        J[1][1]=(b+d)*exp(roughfac*cquadr(vz2-vz4))
        J[1][3]=(b-d)*exp(roughfac*cquadr(vz2+vz4))


        b=(PHI2[0])[1]*(PSI1[0])[1]+(PSI1[0])[3]
        d=(PHI2[0])[2]*(PSI1[0])[2]+(PHI2[0])[0]*(PSI1[0])[0]
        J[0][1]=(b+d)*exp(roughfac*cquadr(vz1-vz4))
        J[0][3]=(b-d)*exp(roughfac*cquadr(vz1+vz4))

        b=(PHI1[0])[2]*(PSI2[0])[2]+(PSI2[0])[0]
        d=(PHI1[0])[3]*(PSI2[0])[3]+(PHI1[0])[1]*(PSI2[0])[1]
        J[1][0]=(b+d)*exp(roughfac*cquadr(vz2-vz3))
        J[1][2]=(b-d)*exp(roughfac*cquadr(vz2+vz3))


        div=(J[0][1]*J[1][0]-J[1][1]*J[0][0])

#        (r[0])[0][0]=(J[1][1]*J[0][2]-J[0][1]*J[1][2])/div # Incoming 1 reflected 1
#        (r[0])[0][1]=(J[1][1]*J[0][3]-J[0][1]*J[1][3])/div # Incoming 2 reflected 1
#        (r[0])[1][0]=(J[0][0]*J[1][2]-J[1][0]*J[0][2])/div # Incoming 1 reflected 2
#        (r[0])[1][1]=(J[0][0]*J[1][3]-J[1][0]*J[0][3])/div # Incoming 2 reflected 2

        (r[0])[0][0]=(J[1][1]*J[0][2]-J[0][1]*J[1][2])/div # Incoming 1 reflected 1
        (r[0])[0][1]=-(J[1][1]*J[0][3]-J[0][1]*J[1][3])/div # Incoming 2 reflected 1
        (r[0])[1][0]=-(J[0][0]*J[1][2]-J[1][0]*J[0][2])/div # Incoming 1 reflected 2
        (r[0])[1][1]=(J[0][0]*J[1][3]-J[1][0]*J[0][3])/div # Incoming 2 reflected 2



        (t[0])[0][0]=J[0][2]*(r[0])[0][0]+J[0][3]*(r[0])[1][0]+J[0][0] # Incoming 1 transmitted 1
        (t[0])[0][1]=J[0][2]*(r[0])[0][1]+J[0][3]*(r[0])[1][1]+J[0][1] # Incoming 2 transmitted 1
        (t[0])[1][0]=J[1][2]*(r[0])[0][0]+J[1][3]*(r[0])[1][0]+J[1][0] # Incoming 1 transmitted 2
        (t[0])[1][1]=J[1][2]*(r[0])[0][1]+J[1][3]*(r[0])[1][1]+J[1][1] # Incoming 2 transmitted 2

        div=(J[0][1]*J[1][0]-J[1][1]*J[0][0])

#        (tprime[0])[0][0]=-J[1][1]/div
#        (tprime[0])[0][1]=J[0][1]/div
#        (tprime[0])[1][0]=J[1][0]/div
#        (tprime[0])[1][1]=-J[0][0]/div
        (tprime[0])[0][0]=-J[1][1]/div
        (tprime[0])[0][1]=-J[0][1]/div
        (tprime[0])[1][0]=-J[1][0]/div
        (tprime[0])[1][1]=-J[0][0]/div

        (rprime[0])[0][0]=J[0][2]*((tprime[0])[0][0])+J[0][3]*((tprime[0])[1][0])
        (rprime[0])[0][1]=J[0][2]*((tprime[0])[0][1])+J[0][3]*((tprime[0])[1][1])
        (rprime[0])[1][0]=J[1][2]*((tprime[0])[0][0])+J[1][3]*((tprime[0])[1][0])
        (rprime[0])[1][1]=J[1][2]*((tprime[0])[0][1])+J[1][3]*((tprime[0])[1][1])


    else:

        b=(PSI1[0])[0]
        d=(PSI1[0])[1]*(PHI1[0])[1]

        J[0][0]=(b+d)*exp(roughfac*cquadr(vz1-vz3))
        J[0][2]=(b-d)*exp(roughfac*cquadr(vz1+vz3))


        b=(PSI2[0])[3]
        d=(PHI2[0])[2]*(PSI2[0])[2]
        J[1][1]=(b+d)*exp(roughfac*cquadr(vz2-vz4))
        J[1][3]=(b-d)*exp(roughfac*cquadr(vz2+vz4))


        (r[0])[0][0]=-J[0][2]/J[0][0] # Incoming 1 reflected 1
        (r[0])[0][1]=0 # Incoming 2 reflected 1
        (r[0])[1][0]=0 #/ Incoming 1 reflected 2
        (r[0])[1][1]=-J[1][3]/J[1][1] # Incoming 2 reflected 2


        (t[0])[0][0]=+J[0][2]*(r[0])[0][0]+J[0][0] # Incoming 1 transmitted 1
        (t[0])[0][1]=0 # Incoming 2 transmitted 1
        (t[0])[1][0]=0 # Incoming 1 transmitted 2
        (t[0])[1][1]=+J[1][3]*(r[0])[1][1]+J[1][1] # Incoming 2 transmitted 2

        (tprime[0])[0][0]=+1.0/J[0][0]
        (tprime[0])[0][1]=0
        (tprime[0])[1][0]=0
        (tprime[0])[1][1]=+1.0/J[1][1]


#        (rprime[0])[0][0]=J[0][2]*((tprime[0])[0][0])
#        (rprime[0])[0][1]=0
#        (rprime[0])[1][0]=0
#        (rprime[0])[1][1]=J[1][3]*((tprime[0])[1][1])
        (rprime[0])[0][0]=J[0][2]*((tprime[0])[0][0])
        (rprime[0])[0][1]=0
        (rprime[0])[1][0]=0
        (rprime[0])[1][1]=J[1][3]*((tprime[0])[1][1])


cdef void Paratt_magnetic_z(Heterostructure* HS, double th, double wavelength, double complex (*rtot)[2][2]) except *:


    if((th<=0)|(th>=90)):
        raise Exception("Theta must be in the range 0<theta<90")

    cdef double k0=two_times_pi/wavelength
    cdef double sintheta=sin(deg_to_rad*th)

    cdef double vy=cos(deg_to_rad*th)
    cdef double vyvy=quadr(vy)
    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
    cdef int Lower, Upper
    cdef int ML_is_diagonal=1
    cdef double complex PHI1[4]
    cdef double complex PSI1[4]
    cdef double complex PHI2[4]
    cdef double complex PSI2[4]
    cdef double complex vz1, vz2, vz3, vz4
    cdef double complex D34,D21,D31
    cdef double complex res, root, b, d, exzexz
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


    p[0][1]=0
    p[1][0]=0
   # print "0"
  #  cdef double complex test=LR[MLCOMP[0][0]].ey
  #  cdef int debug1
   # cdef double complex debug2


    cdef int Cap=NLAYERS-1
    cdef int previously_magnetic=0
   # print "1"
    if(NLAYERS==1):
        Upper=MLCOMP[0][0]
        if(LR[Upper].magdir==3): # magnetic

            MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz1, &vz2, &PHI1, &PHI2, vy, vyvy)

            PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, 1)


            PHI1[0]=1
            PHI1[1]=sintheta
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1
            PHI2[2]=sintheta
            PHI2[1]=0
            PHI2[0]=0
            #Vacuum:
            Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime, 1, vz1, vz2, sintheta,sintheta,LR[MLCOMP[0][0]].Roughness, k0)
        else: # Not magnetic

            NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz1, &vz2, &PHI1, &PHI2, vyvy)
            PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, 0)
        #Vacuum:
            PHI1[0]=1
            PHI1[1]=sintheta
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1
            PHI2[2]=sintheta
            PHI2[1]=0
            PHI2[0]=0
            Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime, 0, vz1, vz2, sintheta,sintheta,LR[MLCOMP[0][0]].Roughness, k0)
    #    print vz1, vz2
    else:
        Lower=MLCOMP[0][0]
        Upper=MLCOMP[1][0]
        if((LR[Lower].magdir==3)|(LR[Upper].magdir==3)):
            if(LR[Lower].magdir==3): # magnetic

                MagneticPhi_z(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, LR[Lower].eg, &vz1, &vz2, &PHI1, &PHI2, vy, vyvy)
                previously_magnetic=1
            else:

                NormalPhi(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, &vz1, &vz2, &PHI1, &PHI2, vyvy)
                previously_magnetic=0
            PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)

            if(LR[Upper].magdir==3):

                MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                previously_magnetic=1
            else:
              #  print "here are the phi"
                NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                previously_magnetic=0
           # p=tbtmatrix(exp(2j*k0*LR[1].Thickness*vz3),0,0,exp(2j*k0*LR[1].Thickness*vz4))
            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
#            print "here p is calculated for the normal structure", p[0][0], p[1][1]
#            print "Components:"
#            print k0, Upper, LR[1].Thickness, vz3, vz4
            Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime,1 , vz1, vz2, vz3,vz4,LR[MLCOMP[0][0]].Roughness, k0)
        #    print "here p is calculated for the normal structure", p[0][0], p[1][1]
       #     print "Components:"
         #   print k0, Upper, LR[Upper].Thickness, vz3, vz4

        else:

          #  print "lowest interface is not magnetic"
            vz1=CalculateVZsigma(vyvy, LR[MLCOMP[0][0]].ex)
            vz2=CalculateVZpi(vyvy, LR[MLCOMP[0][0]].ey, LR[MLCOMP[0][0]].ez)
            vz3=CalculateVZsigma(vyvy, LR[MLCOMP[1][0]].ex)
            vz4=CalculateVZpi(vyvy, LR[MLCOMP[1][0]].ey, LR[MLCOMP[1][0]].ez)
            PHI1[0]=1.
            PHI1[1]=vz3
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1.
            PHI2[2]=vz4/LR[MLCOMP[1][0]].ey
            PHI2[0]=0
            PHI2[1]=0

            p[0][0]=exp(1j*k0*LR[MLCOMP[1][0]].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[MLCOMP[1][0]].Thickness*vz4)
         #   print "here p is calculated for the normal structure", p[0][0], p[1][1]
         #   print "Components:"
         #   print k0, "1", LR[1].Thickness, vz3, vz4

            (rtot[0])[0][0]=((vz3-vz1)/(vz3+vz1))*exp(-2.0*quadr(k0)*vz1*vz3*quadr(LR[MLCOMP[0][0]].Roughness))
            (rtot[0])[1][0]=0
            (rtot[0])[0][1]=0
            (rtot[0])[1][1]=((vz4*LR[MLCOMP[0][0]].ey-vz2*LR[MLCOMP[1][0]].ey)/(vz4*LR[MLCOMP[0][0]].ey+vz2*LR[MLCOMP[1][0]].ey))*exp(-2.0*quadr(k0)*vz2*vz4*quadr(LR[MLCOMP[0][0]].Roughness))


    i=1
    while i<NLAYERS:
       # print "loop start", i
        if(MLLENGTH[i]==1):

            vz1=vz3
            vz2=vz4
             #Inversion of Eigenvectors
            PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)

            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
               # print "i, Upper", i, Upper
                if(LR[Upper].magdir==3):
                  #  print "hallo"

                    MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                 #   print "x"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)
                 #   print "y"
                    if(previously_magnetic):
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    previously_magnetic=0
            else:
              #  print "vacuum"
                PHI1[0]=1. #Eigenvectors
                PHI1[1]=sintheta
                PHI1[2]=0
                PHI1[3]=0
                PHI2[3]=1.
                PHI2[2]=sintheta
                PHI2[0]=0
                PHI2[1]=0
                if(previously_magnetic):
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][0]].Roughness, k0)
                else:
                  #  print "else"
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][0]].Roughness, k0)


            Mult2x2_rightside(&p, rtot)

            Mult2x2_leftside(rtot, &p)

            Mult2x2_rightside(&tprime, rtot)

            Mult2x2_leftside(rtot, &t)

            (rtot[0])[0][0]+=r[0][0]
            (rtot[0])[1][0]+=r[1][0]
            (rtot[0])[0][1]+=r[0][1]
            (rtot[0])[1][1]+=r[1][1]


            if(i!=Cap):
                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)




        else: #A Multilayer

          #  print "1"
            PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
          #  print "2"
            Upper=MLCOMP[i][1]
            vz1=vz3
            vz2=vz4
            if(LR[Upper].magdir==3):
                ML_is_diagonal=0
              #  print "x"
                MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                previously_magnetic=1

            else:
           #     print "Hallo"
                NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                if(previously_magnetic):
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                else:
               #     print "else"
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                previously_magnetic=0

            Mult2x2_leftside(&t_ML_back_1, &p) # t'(AB)*p(A)
            Mult2x2_rightside(&p, &t_ML_in_1) # p(A) * t(AB)


            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)

            j=1
            while j<MLLENGTH[i]:
                Upper=MLCOMP[i][(j+1)%MLLENGTH[i]]
                vz1=vz3
                vz2=vz4
                PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "5"
              #  Upper=MLCOMP[i][1]
                if(LR[Upper].magdir==3):
                  #  print "hallo"
                    MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)
                    ML_is_diagonal=0
                #    print "y"
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    previously_magnetic=0

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


                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
           #     if(j==1):
            #        print "p C components"
            #        print Upper, LR[Upper].Thickness, vz3, vz4

                j=j+1


            if(ML_is_diagonal):
                r_ML_in_2[0][0]=r_ML_in_1[0][0]*(1-ipow(t_ML_in_1[0][0]*t_ML_back_1[0][0], MLREP[i]-1))/(1-t_ML_in_1[0][0]*t_ML_back_1[0][0])
                r_ML_in_2[1][1]=r_ML_in_1[1][1]*(1-ipow(t_ML_in_1[1][1]*t_ML_back_1[1][1], MLREP[i]-1))/(1-t_ML_in_1[1][1]*t_ML_back_1[1][1])

                t_ML_in_1[0][0]=ipow(t_ML_in_1[0][0],  MLREP[i]-1)
                t_ML_in_1[1][1]=ipow(t_ML_in_1[1][1],  MLREP[i]-1)
                t_ML_back_1[0][0]=ipow(t_ML_back_1[0][0],  MLREP[i]-1)
                t_ML_back_1[1][1]=ipow(t_ML_back_1[1][1],  MLREP[i]-1)

            else:

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
                vz1=vz3
                vz2=vz4
                Upper=MLCOMP[i][j]
                PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "10"
                if(LR[Upper].magdir==3):
                  #  print "hallo"
                    MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)


                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    previously_magnetic=0
              #  print "11"
                Mult2x2_rightside(&p, rtot)
                Mult2x2_leftside(rtot, &p)
                Mult2x2_rightside(&tprime, rtot)
                Mult2x2_leftside(rtot, &t)
                (rtot[0])[0][0]+=r[0][0]
                (rtot[0])[1][0]+=r[1][0]
                (rtot[0])[0][1]+=r[0][1]
                (rtot[0])[1][1]+=r[1][1]


                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
                j=j+1
            vz1=vz3
            vz2=vz4
            if(i==(NLAYERS-1)):
                PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "13"
              #  vz3=CalculateVZsigma(vyvy, LR[Upper].ex)
              #  vz4=CalculateVZpi(vyvy, LR[Upper].ey, LR[Upper].ez)
                PHI1[0]=1. #Eigenvectors
                PHI1[1]=sintheta
                PHI1[2]=0
                PHI1[3]=0
                PHI2[3]=1.
                PHI2[2]=sintheta
                PHI2[0]=0
                PHI2[1]=0
                if(previously_magnetic):
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                else:
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
              #  print "14"
            else:
                Upper=MLCOMP[i+1][0]
                PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
            #    print "15"
                if(LR[Upper].magdir==3):
                  #  print "hallo"
                    MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    previously_magnetic=0
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


            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
         #   print "18"
        i=i+1

cdef void Paratt_magnetic_z_MS(Heterostructure* HS, double th, double wavelength, double complex (*rtot)[2][2]) except *:


    if((th<=0)|(th>=90)):
        raise Exception("Theta must be in the range 0<theta<90")

    cdef double k0=two_times_pi/wavelength
    cdef double sintheta=sin(deg_to_rad*th)

    cdef double vy=cos(deg_to_rad*th)
    cdef double vyvy=quadr(vy)
    cdef int NLAYERS=(HS[0]).NLayers
    cdef int* MLLENGTH=(HS[0]).MLLENGTH
    cdef int** MLCOMP=(HS[0]).MLCOMP
    cdef int* MLREP=(HS[0]).MLREP
    cdef CLayer* LR=(HS[0]).LR
    cdef int i,j
    cdef int Lower, Upper
    cdef int ML_is_diagonal=1
    cdef double complex PHI1[4]
    cdef double complex PSI1[4]
    cdef double complex PHI2[4]
    cdef double complex PSI2[4]
    cdef double complex vz1, vz2, vz3, vz4
    cdef double complex D34,D21,D31
    cdef double complex res, root, b, d, exzexz
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

    p[0][1]=0
    p[1][0]=0
   # print "0"
   # cdef double complex test=LR[MLCOMP[0][0]].ey

    cdef int Cap=NLAYERS-1
    cdef int previously_magnetic=0
   # print "1"
    if(NLAYERS==1):
        Upper=MLCOMP[0][0]
        if(LR[Upper].magdir==3): # magnetic

            MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz1, &vz2, &PHI1, &PHI2, vy, vyvy)

            PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, 1)


            PHI1[0]=1
            PHI1[1]=sintheta
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1
            PHI2[2]=sintheta
            PHI2[1]=0
            PHI2[0]=0
            #Vacuum:
            Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime, 1, vz1, vz2, sintheta,sintheta,LR[MLCOMP[0][0]].Roughness, k0)
        else: # Not magnetic

            NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz1, &vz2, &PHI1, &PHI2, vyvy)
            PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, 0)
        #Vacuum:
            PHI1[0]=1
            PHI1[1]=sintheta
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1
            PHI2[2]=sintheta
            PHI2[1]=0
            PHI2[0]=0
            Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime, 0, vz1, vz2, sintheta,sintheta,LR[MLCOMP[0][0]].Roughness, k0)
    #    print vz1, vz2
    else:
        Lower=MLCOMP[0][0]
        Upper=MLCOMP[1][0]
        if((LR[Lower].magdir==3)|(LR[Upper].magdir==3)):
            if(LR[Lower].magdir==3): # magnetic

                MagneticPhi_z(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, LR[Lower].eg, &vz1, &vz2, &PHI1, &PHI2, vy, vyvy)
                previously_magnetic=1
            else:

                NormalPhi(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, &vz1, &vz2, &PHI1, &PHI2, vyvy)
                previously_magnetic=0
            PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)

            if(LR[Upper].magdir==3):

                MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                previously_magnetic=1
            else:
              #  print "here are the phi"
                NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                previously_magnetic=0
           # p=tbtmatrix(exp(2j*k0*LR[1].Thickness*vz3),0,0,exp(2j*k0*LR[1].Thickness*vz4))
            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
#            print "here p is calculated for the normal structure", p[0][0], p[1][1]
#            print "Components:"
#            print k0, Upper, LR[1].Thickness, vz3, vz4
            Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, rtot, &rprime, &t, &tprime, 1 , vz1, vz2, vz3,vz4,LR[MLCOMP[0][0]].Roughness, k0)
        #    print "here p is calculated for the normal structure", p[0][0], p[1][1]
       #     print "Components:"
         #   print k0, Upper, LR[Upper].Thickness, vz3, vz4

        else:

          #  print "lowest interface is not magnetic"
            vz1=CalculateVZsigma(vyvy, LR[MLCOMP[0][0]].ex)
            vz2=CalculateVZpi(vyvy, LR[MLCOMP[0][0]].ey, LR[MLCOMP[0][0]].ez)
            vz3=CalculateVZsigma(vyvy, LR[MLCOMP[1][0]].ex)
            vz4=CalculateVZpi(vyvy, LR[MLCOMP[1][0]].ey, LR[MLCOMP[1][0]].ez)
            PHI1[0]=1.
            PHI1[1]=vz3
            PHI1[2]=0
            PHI1[3]=0
            PHI2[3]=1.
            PHI2[2]=vz4/LR[MLCOMP[1][0]].ey
            PHI2[0]=0
            PHI2[1]=0

            p[0][0]=exp(1j*k0*LR[MLCOMP[1][0]].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[MLCOMP[1][0]].Thickness*vz4)

            (rtot[0])[0][0]=((vz3-vz1)/(vz3+vz1))*exp(-2.0*quadr(k0)*vz1*vz3*quadr(LR[MLCOMP[0][0]].Roughness))
            (rtot[0])[1][0]=0
            (rtot[0])[0][1]=0
            (rtot[0])[1][1]=((vz4*LR[MLCOMP[0][0]].ey-vz2*LR[MLCOMP[1][0]].ey)/(vz4*LR[MLCOMP[0][0]].ey+vz2*LR[MLCOMP[1][0]].ey))*exp(-2.0*quadr(k0)*vz2*vz4*quadr(LR[MLCOMP[0][0]].Roughness))

    i=1
    while i<NLAYERS:
       # print "loop start", i
        if(MLLENGTH[i]==1):

            vz1=vz3
            vz2=vz4
             #Inversion of Eigenvectors
            PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)

            if(i!=Cap):
                Upper=MLCOMP[i+1][0]
               # print "i, Upper", i, Upper
                if(LR[Upper].magdir==3):
                  #  print "hallo"

                    MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                 #   print "x"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)
                 #   print "y"
                    if(previously_magnetic):
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                    previously_magnetic=0
            else:
              #  print "vacuum"
                PHI1[0]=1. #Eigenvectors
                PHI1[1]=sintheta
                PHI1[2]=0
                PHI1[3]=0
                PHI2[3]=1.
                PHI2[2]=sintheta
                PHI2[0]=0
                PHI2[1]=0
                if(previously_magnetic):
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][0]].Roughness, k0)
                else:
                  #  print "else"
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][0]].Roughness, k0)

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
#            print "rtot ", i
#            print (rtot[0])[0][0]
#            print (rtot[0])[0][1]
#            print (rtot[0])[1][0]
#            print (rtot[0])[1][1]


            if(i!=Cap):
                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
        else:
            PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
          #  print "2"
            Upper=MLCOMP[i][1]
            vz1=vz3
            vz2=vz4

            if(LR[Upper].magdir==3):
                ML_is_diagonal=0
              #  print "x"
                MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)

                previously_magnetic=1

            else:
           #     print "Hallo"
                NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                if(previously_magnetic):
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)
                else:
               #     print "else"
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r_ML_in_1, &rprime, &t_ML_in_1, &t_ML_back_1, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][0]].Roughness, k0)

                previously_magnetic=0


            Mult2x2_leftside(&t_ML_back_1, &p) # t'(AB)*p(A)


            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)

#            print "r_ML_in_1 0"
#            print r_ML_in_1[0][0]
#            print r_ML_in_1[0][1]
#            print r_ML_in_1[1][0]
#            print r_ML_in_1[1][1]



            j=1
            while j<MLLENGTH[i]:
                Upper=MLCOMP[i][(j+1)%MLLENGTH[i]]
              #  print "j and Upper", j, (j+1)%MLLENGTH[i], Upper
                vz1=vz3
                vz2=vz4
                PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "5"
              #  Upper=MLCOMP[i][1]

                if(LR[Upper].magdir==3):
                  #  print "hallo"
                    MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)
                    ML_is_diagonal=0
                #    print "y"
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j]].Roughness, k0)
                    previously_magnetic=0

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
                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
                j=j+1
         #   print vz1, vz2, vz3, vz4

            p[0][0]=exp(1j*k0*LR[MLCOMP[i][MLLENGTH[i]-1]].Thickness*vz1)
            p[1][1]=exp(1j*k0*LR[MLCOMP[i][MLLENGTH[i]-1]].Thickness*vz2)
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
                Lower=MLCOMP[i][j-1]
                Upper=MLCOMP[i][j]


                if(LR[Lower].magdir==3):
                    MagneticPhi_z(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, LR[Lower].eg, &vz1, &vz2, &PHI1, &PHI2, vy, vyvy)
                    ML_is_diagonal=0
                    previously_magnetic=1
                else:
                    NormalPhi(LR[Lower].ex, LR[Lower].ey, LR[Lower].ez, &vz1, &vz2, &PHI1, &PHI2, vyvy)
                    previously_magnetic=0


                PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)

                if(LR[Upper].magdir==3):
                  #  print "hallo"
                    MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)
                    ML_is_diagonal=0
                #    print "y"
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[Lower].Roughness, k0)
                  #  previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[Lower].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[Lower].Roughness, k0)



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



               # Mult2x2_rightside(&p, &t_ML_in_1)
                Mult2x2_rightside(&C1, &t_ML_in_1)
             #   Mult2x2_rightside(&p, &t_ML_in_1)
                Mult2x2_rightside(&t, &t_ML_in_1)
                p[0][0]=exp(1j*k0*LR[Lower].Thickness*vz1)
                p[1][1]=exp(1j*k0*LR[Lower].Thickness*vz2)

#                if(j==2):
#                    print "p B"
#                elif(j==1):
#                    print "p A"
#                print p[0][0]
#                print p[1][1]

                Mult2x2_rightside(&p, &t_ML_in_1)

                j=j-1


            Mult2x2_rightside(&p, &r_ML_back_1)
            Mult2x2_leftside(&r_ML_back_1, &p)


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


            Mult2x2_rightside(&t_ML_back_2, rtot) #Here the ML part is missing
            Mult2x2_leftside(rtot, &C0)
            Mult2x2_leftside(rtot, &t_ML_in_2) # (t'(CA) p_C t'(BC) p_B t'(AB)*p(A))^N rtot (p(A) * t(AB) p_B t(BC) p_C t(CA))^N
            (rtot[0])[0][0]+=r_ML_in_2[0][0]
            (rtot[0])[1][0]+=r_ML_in_2[1][0]
            (rtot[0])[0][1]+=r_ML_in_2[0][1]
            (rtot[0])[1][1]+=r_ML_in_2[1][1]



            Upper=MLCOMP[i][0]
            if(LR[Upper].magdir==3):
              #  print "hallo"
                MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                previously_magnetic=1
            else:
           #     print "Hallo"
                NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                previously_magnetic=0
            j=1

            while j<MLLENGTH[i]:

                vz1=vz3
                vz2=vz4
                Upper=MLCOMP[i][j]
                PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "10"
                if(LR[Upper].magdir==3):
                  #  print "hallo"
                    MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)


                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][j-1]].Roughness, k0)
                    previously_magnetic=0
           #     print "4", vz1, vz2, vz3, vz4
              #  print "11"
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

                p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
                p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)
                j=j+1
              #  print "12"

            if(i==(NLAYERS-1)):
                vz1=vz3
                vz2=vz4
                PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
              #  print "13"
              #  vz3=CalculateVZsigma(vyvy, LR[Upper].ex)
              #  vz4=CalculateVZpi(vyvy, LR[Upper].ey, LR[Upper].ez)
                PHI1[0]=1. #Eigenvectors
                PHI1[1]=sintheta
                PHI1[2]=0
                PHI1[3]=0
                PHI2[3]=1.
                PHI2[2]=sintheta
                PHI2[0]=0
                PHI2[1]=0
                if(previously_magnetic):
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                else:
                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, sintheta, sintheta, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
              #  print "14"

            else:
                vz1=vz3
                vz2=vz4
                Upper=MLCOMP[i+1][0]
                PHI_to_PSI_z(&PSI1, &PHI1, &PSI2, &PHI2, previously_magnetic)
            #    print "15"
                if(LR[Upper].magdir==3):
                  #  print "hallo"
                    MagneticPhi_z(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, LR[Upper].eg, &vz3, &vz4, &PHI1, &PHI2, vy, vyvy)

                    Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    previously_magnetic=1
                else:
               #     print "Hallo"
                    NormalPhi(LR[Upper].ex, LR[Upper].ey, LR[Upper].ez, &vz3, &vz4, &PHI1, &PHI2, vyvy)

                    if(previously_magnetic):
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 1, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    else:
                   #     print "else"
                        Calculate_rt_z(&PSI1, &PHI1, &PSI2, &PHI2, &r, &rprime, &t, &tprime, 0, vz1, vz2, vz3, vz4, LR[MLCOMP[i][MLLENGTH[i]-1]].Roughness, k0)
                    previously_magnetic=0
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


            p[0][0]=exp(1j*k0*LR[Upper].Thickness*vz3)
            p[1][1]=exp(1j*k0*LR[Upper].Thickness*vz4)

        i=i+1






cdef int FindLayerNumber(string, MaxLayer):

    if(string=="default"):
        return MaxLayer

    a=string.split(",")

    N=len(a)
    j=0
    IsInBracket=False
    for i in range(N):
        if (a[i].find('(')!=-1):
            j=j+1
            IsInBracket=True
    #        print a[i], "w", j
        elif (a[i].find(')')!=-1):
            IsInBracket=False
   #         print a[i], "x", j
        else:

            if(IsInBracket==False):
                j=j+1

    return j

cdef MakeConsistencyCheck(string, Heterostructure *STR, MaxLayer):
    #Check for brackets, Maximum layer number and substrate
    Ntot=len(string)
    IsInBracket=False
    for i in range(Ntot):
        if(IsInBracket):
            if(string[i]=='('):
                raise Exception("Brackets inside Brackets are not supported.")
            if(string[i]==')'):
                IsInBracket=False
        else:
            if(string[i]==')'):
                raise Exception("Something is wrong with the Multilayer brackets")
            if(string[i]=='('):
                IsInBracket=True
    if(IsInBracket):
        raise Exception("Something is wrong with the Multilayer brackets")


    for i in range((STR[0]).NLayers):
        for j in range((STR[0]).MLLENGTH[i]):
            if((STR[0]).MLCOMP[i][j]>=MaxLayer):
                raise Exception("Layer " + str((STR[0]).MLCOMP[i][j]) + " in Multilayer structure string is not defined")

    if((STR[0]).MLLENGTH[0]!=1):
        raise Exception("Substrate as Multilayer is ill-defined")

    return 0



cdef int FindComposition(string, int MaxLayer, int LayerNumber, Heterostructure *STR):


  #  LayerNumber=FindLayerNumber(string, MaxLayer)
    (STR[0]).NLayers=LayerNumber

  #  print LayerNumber

    if(string=="default"):
        for i in range(MaxLayer):
            (STR[0]).MLREP[i]=1
            (STR[0]).MLLENGTH[i]=1
            (STR[0]).MLCOMP[i]= <int*>malloc(sizeof(int))
            (STR[0]).MLCOMP[i][0]=i
        return 0

    SaveNumbers=[]

  #  print "not default"
    a=string.split(",")

    N=len(a)
    cdef int j=0
    IsInBracket=False
    cdef int k=0
    cdef int l


    for i in range(N):
        if (a[i].find('(')!=-1):
            b=a[i].split('*(')
            (STR[0]).MLREP[j]=int(b[0])
         #   print "j", j
            SaveNumbers=SaveNumbers+[int(b[1])]
            IsInBracket=True
            k=1
        elif (a[i].find(')')!=-1):
            b=a[i].split(')')
            SaveNumbers=SaveNumbers+[int(b[0])]
            IsInBracket=False
            k=k+1
            (STR[0]).MLLENGTH[j]=k
            (STR[0]).MLCOMP[j]= <int*>malloc(int((STR[0]).MLLENGTH[j])*sizeof(int))
         #   print "j", j
            for l in range((STR[0]).MLLENGTH[j]):
                (STR[0]).MLCOMP[j][l]=int(SaveNumbers[l])
            SaveNumbers=[]
            j=j+1
        #    print a[i], "x", j
        else:
            if(IsInBracket):
                k=k+1
                SaveNumbers=SaveNumbers+[int(a[i])]
            else:
                (STR[0]).MLREP[j]=1
                (STR[0]).MLLENGTH[j]=1
                (STR[0]).MLCOMP[j]= <int*>malloc(sizeof(int))
                (STR[0]).MLCOMP[j][0]=int(a[i])
                j=j+1




    MakeConsistencyCheck(string, &(STR[0]), MaxLayer)




def Generate_structure(int NLayers_types, MLstructure="default"):




#    cdef Heterostructure STR
    if(NLayers_types<=0):
        raise Exception("Please generate at least one layer!")

    cdef int i
    HS=[0 for i in range(NLayers_types)]
    HS[0]=Lowestlayer(0,0,MLstructure, NLayers_types)
    for i in range(1,NLayers_types):
        HS[i]=Layer(0,0)

  #  HS[0].initstructure(MLstructure, NLayers_types)

    return HS

def Reflectivity(HS, th, wavelength, MultipleScattering=1, MagneticCutoff=1.0E-6, Output="T"):

    if( HS[0].isthisthelowestlayer()!=1   ):
        raise Exception("Underlying structure not initialized. Please generate the layer list with Generate_Structure!")

    if(hasattr(th, "__len__")==False):
        th=np.array([th])
    cdef int NAngles=len(th)

    cdef long a=(HS[0].motherpointer())
    cdef Heterostructure* A=<Heterostructure*>a #Recover the C storage of the structure
    cdef CLayer* B

    cdef int NLayers=(A[0]).NLayers

    cdef int NLayers_types=(A[0]).NLayers_types
    cdef int i


    cdef int allx=0
    cdef int ally=0
    cdef int allz=0

    cdef double Cutoffquad=quadr(MagneticCutoff)

    cdef int Setting1
    cdef int Setting2
    cdef int Setting3=0 #0 Means Heterostructure has no magnetic parts

    if(Output=="s"): #return complex rss
        Setting1=1
        Setting2=0
    elif(Output=="S"): #return real Rss
        Setting1=1
        Setting2=1
    elif(Output=="p"): #return complex rpp
        Setting1=2
        Setting2=0
    elif(Output=="P"): #return real Rpp
        Setting1=2
        Setting2=1
    elif(Output=="t"): #return complex rss,rsp,rps,rpp
        Setting1=3
        Setting2=0
    elif(Output=="T"): #return real Rs, Rp, Rl, Rr
        Setting1=3
        Setting2=1
    for i in range(NLayers_types):  #Download all layers that the user has filled into the Structure
        a=(HS[i]).layercontentpointer()
        B=<CLayer*>a
       # LR[i]=B[0]
        (A[0]).LR[i]=B[0]
       # print "content", (A[0]).LR[i].ex
        if(B[0].type==3):
            Setting3=1
            if(B[0].magdir==1):
                allx=1
            elif(B[0].magdir==2):
                ally=1
            elif(B[0].magdir==3):
                allz=1
            if(Setting1!=3):
                raise Exception('Exception! Magnetic heterostructures must have "t" or "T" as an output parameter')
            if(cabsquadr(B[0].eg)<Cutoffquad):#Apply the magnetic cutoff

               # HS[i].seteps([HS[i].epsxx(), HS[i].epsyy(), HS[i].epszz()])

                (A[0]).LR[i].eg=0
                (A[0]).LR[i].type=2
                (A[0]).LR[i].magdir=0




    if((ally&allz)|(ally&allx)|(allx&allz)):
        raise Exception('Exception! Multiple magnetization directions are so far not supported!')







    cdef double complex rss, rpp
    cdef double complex rmat[2][2]

    if(Setting3==0):#Not magnetic

        if(Setting1!=2): #2 means only pi


            if(MultipleScattering):

                if(Setting2): #means absquadr
                    routs=np.zeros(NAngles)
                    for i in range(NAngles):
                        rss=LinDicParatt_Sigma_MS(A, th[i], wavelength)
                        routs[i]=cabsquadr(rss)
                else:
                    routs=np.zeros(NAngles, dtype=complex)
                    for i in range(NAngles):
                        rss=LinDicParatt_Sigma_MS(A, th[i], wavelength)
                        routs[i]=rss
            else:

                if(Setting2): #means absquadr
                    routs=np.zeros(NAngles)

                    for i in range(NAngles):
                        rss=LinDicParatt_Sigma(A, th[i], wavelength)
                        routs[i]=cabsquadr(rss)
                else:
                    routs=np.zeros(NAngles, dtype=complex)
                    for i in range(NAngles):
                        rss=LinDicParatt_Sigma(A, th[i], wavelength)
                        routs[i]=rss
        if(Setting1!=1): #1 means only sigma
            if(MultipleScattering):
                if(Setting2): #means absquadr
                    routp=np.zeros(NAngles)
                    for i in range(NAngles):
                        rpp=LinDicParatt_Pi_MS(A, th[i], wavelength)
                        routp[i]=cabsquadr(rpp)
                else:
                    routp=np.zeros(NAngles, dtype=complex)
                    for i in range(NAngles):
                        rpp=LinDicParatt_Pi_MS(A, th[i], wavelength)
                        routp[i]=rpp
            else:
                if(Setting2): #means absquadr
                    routp=np.zeros(NAngles)
                    for i in range(NAngles):
                        rpp=LinDicParatt_Pi(A, th[i], wavelength)
                        routp[i]=cabsquadr(rpp)
                else:
                    routp=np.zeros(NAngles, dtype=complex)
                    for i in range(NAngles):
                        rpp=LinDicParatt_Pi(A, th[i], wavelength)
                        routp[i]=rpp
    else:#Magnetic
        if(allx):
            if(MultipleScattering):
                if(Setting2): #means absquadr
                    routs=np.zeros(NAngles)
                    routp=np.zeros(NAngles)
                    rempty=np.zeros(NAngles)
                    for i in range(NAngles):
                        rss=LinDicParatt_Sigma_MS(A, th[i], wavelength)
                        rpp=LinDicParatt_Pi_xmag_MS(A, th[i], wavelength)
                        routs[i]=cabsquadr(rss)
                        routp[i]=cabsquadr(rpp)
                else:
                    routs=np.zeros(NAngles, dtype=complex)
                    routp=np.zeros(NAngles, dtype=complex)
                    rempty=np.zeros(NAngles)
                    for i in range(NAngles):
                        rss=LinDicParatt_Sigma_MS(A, th[i], wavelength)
                        rpp=LinDicParatt_Pi_xmag_MS(A, th[i], wavelength)
                        routs[i]=rss
                        routp[i]=rpp
            else:
                if(Setting2): #means absquadr
                    routs=np.zeros(NAngles)
                    routp=np.zeros(NAngles)
                    rempty=np.zeros(NAngles)
                    for i in range(NAngles):
                        rss=LinDicParatt_Sigma(A, th[i], wavelength)
                        rpp=LinDicParatt_Pi_xmag(A, th[i], wavelength)
                        routs[i]=cabsquadr(rss)
                        routp[i]=cabsquadr(rpp)
                else:
                    routs=np.zeros(NAngles, dtype=complex)
                    routp=np.zeros(NAngles, dtype=complex)
                    rempty=np.zeros(NAngles)
                    for i in range(NAngles):
                        rss=LinDicParatt_Sigma(A, th[i], wavelength)
                        rpp=LinDicParatt_Pi_xmag(A, th[i], wavelength)
                        routs[i]=rss
                        routp[i]=rpp
        elif(ally):
            if(MultipleScattering):
                if(Setting2): #means absquadr
                    routs=np.zeros(NAngles)
                    routp=np.zeros(NAngles)
                    routl=np.zeros(NAngles)
                    routr=np.zeros(NAngles)
                    for i in range(NAngles):
                        Paratt_magnetic_y_MS(A, th[i], wavelength, &rmat)
                        routs[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[1][0])
                        routp[i]=cabsquadr(rmat[1][0])+cabsquadr(rmat[1][1])
                        routl[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
                        routr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                else:
                    routs=np.zeros(NAngles, dtype=complex)
                    routp=np.zeros(NAngles, dtype=complex)
                    routl=np.zeros(NAngles, dtype=complex)
                    routr=np.zeros(NAngles, dtype=complex)
                    for i in range(NAngles):
                        Paratt_magnetic_y_MS(A, th[i], wavelength, &rmat)
                        routs[i]=rmat[0][0]
                        routp[i]=rmat[0][1]
                        routl[i]=rmat[1][0]
                        routr[i]=rmat[1][1]
            else:
                if(Setting2): #means absquadr
                    routs=np.zeros(NAngles)
                    routp=np.zeros(NAngles)
                    routl=np.zeros(NAngles)
                    routr=np.zeros(NAngles)
                    for i in range(NAngles):
                        Paratt_magnetic_y(A, th[i], wavelength, &rmat)
                        routs[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[1][0])
                        routp[i]=cabsquadr(rmat[1][0])+cabsquadr(rmat[1][1])
                        routl[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
                        routr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                else:
                    routs=np.zeros(NAngles, dtype=complex)
                    routp=np.zeros(NAngles, dtype=complex)
                    routl=np.zeros(NAngles, dtype=complex)
                    routr=np.zeros(NAngles, dtype=complex)
                    for i in range(NAngles):
                        Paratt_magnetic_y(A, th[i], wavelength, &rmat)
                        routs[i]=rmat[0][0]
                        routp[i]=rmat[0][1]
                        routl[i]=rmat[1][0]
                        routr[i]=rmat[1][1]
        elif(allz):
            if(MultipleScattering):
                if(Setting2): #means absquadr
                    routs=np.zeros(NAngles)
                    routp=np.zeros(NAngles)
                    routl=np.zeros(NAngles)
                    routr=np.zeros(NAngles)
                    for i in range(NAngles):
                        Paratt_magnetic_z_MS(A, th[i], wavelength, &rmat)
                        routs[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[1][0])
                        routp[i]=cabsquadr(rmat[1][0])+cabsquadr(rmat[1][1])
                        routl[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
                        routr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                else:
                    routs=np.zeros(NAngles, dtype=complex)
                    routp=np.zeros(NAngles, dtype=complex)
                    routl=np.zeros(NAngles, dtype=complex)
                    routr=np.zeros(NAngles, dtype=complex)
                    for i in range(NAngles):
                        Paratt_magnetic_z_MS(A, th[i], wavelength, &rmat)
                        routs[i]=rmat[0][0]
                        routp[i]=rmat[0][1]
                        routl[i]=rmat[1][0]
                        routr[i]=rmat[1][1]
            else:
                if(Setting2): #means absquadr
                    routs=np.zeros(NAngles)
                    routp=np.zeros(NAngles)
                    routl=np.zeros(NAngles)
                    routr=np.zeros(NAngles)
                    for i in range(NAngles):
                        Paratt_magnetic_z(A, th[i], wavelength, &rmat)
                        routs[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[1][0])
                        routp[i]=cabsquadr(rmat[1][0])+cabsquadr(rmat[1][1])
                        routl[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
                        routr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                else:
                    routs=np.zeros(NAngles, dtype=complex)
                    routp=np.zeros(NAngles, dtype=complex)
                    routl=np.zeros(NAngles, dtype=complex)
                    routr=np.zeros(NAngles, dtype=complex)
                    for i in range(NAngles):
                        Paratt_magnetic_z(A, th[i], wavelength, &rmat)
                        routs[i]=rmat[0][0]
                        routp[i]=rmat[0][1]
                        routl[i]=rmat[1][0]
                        routr[i]=rmat[1][1]


    if(allx):
        if(Setting2):
            rempty=0.5*(routp+routs)

    if(Setting3==0):
        if(Setting1==1):
            return routs
        elif(Setting1==2):
            return routp
        elif(Setting1==3):
            return [routs, routp]
    else:
        if(allx):
            if(Setting2):
                return [routs, routp, rempty, rempty]
            else:
                return [routs, rempty, rempty, routp]
        else:
            return [routs, routp, routl, routr]






