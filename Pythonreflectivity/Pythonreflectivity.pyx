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

#cython: boundscheck=False, wraparound=False, nonecheck=False
import cython
from libc.stdlib cimport malloc, free
import ctypes
from numpy cimport ndarray
from numpy import zeros, array
from libc.string cimport memcpy
from Structural cimport *
from Mathematical_Functions_Reflectivity cimport *
from Multilayer_Functions_Reflectivity cimport *
from Reflectivity_Sigma cimport *
from Reflectivity_Pi cimport *
from MOKE_transversal cimport *
from MOKE_longitudinal cimport *
from MOKE_polar cimport *
from Full_Matrix cimport *
from Kramers_Kronig cimport *
#cdef struct CLayer:
#    double Thickness, Roughness
#    double complex cx, cy, cz, cg
#    double complex cxy, cyx, cxz, czx, cyz, czy
#    int type, magdir


cdef class Layer:
    cdef CLayer Layercontent
    cdef Heterostructure Mother
    cdef int islowestlayer
    def __cinit__(self, double d, double sigma, MLstructure, int NL_types, int IsLowest):

        self.Layercontent.Thickness = d
        self.Layercontent.Roughness = sigma
        self.Layercontent.type=1
        self.Layercontent.magdir = 0
        self.Layercontent.cx=0
        self.Layercontent.cy=0
        self.Layercontent.cz=0
        self.Layercontent.cg=0
        self.islowestlayer=IsLowest
        cdef int NLayers
        cdef int *MLLENGTH
        cdef int *MLREP
        cdef int **MLCOMP
        cdef CLayer *LR
        if(self.islowestlayer):
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
        if(self.islowestlayer):
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
    def setchi(self, chi):

        if(hasattr(chi, "__len__")):

            if(len(chi) in (1,3,4,9)):

                if(len(chi)==1):
                    self.Layercontent.type=1
                    #self.Layercontent.cps[0]=chi[0]
                    self.Layercontent.cx=chi[0]
                    self.Layercontent.cy=chi[0]
                    self.Layercontent.cz=chi[0]
                    self.Layercontent.cg=0
                    self.Layercontent.magdir=0
                elif(len(chi)==3):
                    self.Layercontent.type=2
                    self.Layercontent.cx=chi[0]
                    self.Layercontent.cy=chi[1]
                    self.Layercontent.cz=chi[2]
                    self.Layercontent.cg=0
                    self.Layercontent.magdir=0


                elif(len(chi)==4):
                    self.Layercontent.type=3
                    self.Layercontent.cx=chi[0]
                    self.Layercontent.cy=chi[1]
                    self.Layercontent.cz=chi[2]
                    if(chi[3]==0):
                        self.Layercontent.magdir=0
                        self.Layercontent.type=2
                    elif(self.Layercontent.magdir==0):
                        raise Exception("Please use setmag to set a magnetization direction for this layer before you set a magnetooptic constant")
                    self.Layercontent.cg=chi[3]
                    if( self.Layercontent.magdir== 1):
                        self.Layercontent.cyz=-chi[3]
                        self.Layercontent.czy=chi[3]
                    elif( self.Layercontent.magdir== 2):
                        self.Layercontent.cxz=chi[3]
                        self.Layercontent.czx=-chi[3]
                    elif( self.Layercontent.magdir== 3):
                        self.Layercontent.cxy=chi[3]
                        self.Layercontent.cyx=-chi[3]

                elif(len(chi)==9 ):

                    self.Layercontent.type=4
                    self.Layercontent.cx=chi[0]
                    self.Layercontent.cxy=chi[1]
                    self.Layercontent.cxz=chi[2]
                    self.Layercontent.cyx=chi[3]
                    self.Layercontent.cy=chi[4]
                    self.Layercontent.cyz=chi[5]
                    self.Layercontent.czx=chi[6]
                    self.Layercontent.czy=chi[7]
                    self.Layercontent.cz=chi[8]

            else:
                raise Exception('chi must be a number or an array of length 1, 3 or 4')
        else:

            self.Layercontent.type=1
            self.Layercontent.cx=chi
            self.Layercontent.cy=chi
            self.Layercontent.cz=chi
    def chi(self, wl=-1):
        if(self.Layercontent.type==1):
            return self.Layercontent.cx
        elif(self.Layercontent.type==2):
            return array([ self.Layercontent.cx, self.Layercontent.cy, self.Layercontent.cz ])
        elif(self.Layercontent.type==3):
            return array([ self.Layercontent.cx, self.Layercontent.cy, self.Layercontent.cz , self.Layercontent.cg ])
        elif(self.Layercontent.type==4):
            return array([ self.Layercontent.cx, self.Layercontent.cxy, self.Layercontent.cxz ,  self.Layercontent.cyx, self.Layercontent.cy, self.Layercontent.cyz ,  self.Layercontent.czx, self.Layercontent.czy, self.Layercontent.cz  ])

    cpdef double complex chixx(self):
        return self.Layercontent.cx
    cpdef double complex chiyy(self):
        return self.Layercontent.cy
    cpdef double complex chizz(self):
        return self.Layercontent.cz
    cpdef double complex chig(self):
        return self.Layercontent.cg

    cpdef double complex chixy(self):
        return self.Layercontent.cxy
    cpdef double complex chiyx(self):
        return self.Layercontent.cyx

    cpdef double complex chixz(self):
        return self.Layercontent.cxz
    cpdef double complex chizx(self):
        return self.Layercontent.czx

    cpdef double complex chizy(self):
        return self.Layercontent.czy
    cpdef double complex chiyz(self):
        return self.Layercontent.cyz
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
            elif(self.Layercontent.magdir==4):
                return "The material is described by a full chi matrix"
            else:
                return "No magnetization has been set"
        else:
            return "This layer is not magnetic"
    cdef int GetType(self):
        return self.Layercontent.type
    cdef int dir(self):
        return self.Layercontent.magdir
    def setchixx(self, x):
        self.Layercontent.cx=x
    def setchiyy(self, x):
        self.Layercontent.cy=x
    def setchizz(self, x):
        self.Layercontent.cz=x

    def setchixy(self, x):
        self.Layercontent.cxy=x
        self.Layercontent.type=4
    def setchiyx(self, x):
        self.Layercontent.cyx=x
        self.Layercontent.type=4
    def setchixz(self, x):
        self.Layercontent.cxz=x
        self.Layercontent.type=4
    def setchizx(self, x):
        self.Layercontent.czx=x
        self.Layercontent.type=4
    def setchiyz(self, x):
        self.Layercontent.cyz=x
        self.Layercontent.type=4
    def setchizy(self, x):
        self.Layercontent.czy=x
        self.Layercontent.type=4

    def setchig(self,  x):
        if(x==0):
            self.Layercontent.magdir=0
            self.Layercontent.type=2
        elif(self.Layercontent.magdir==0):
            raise Exception("Please use setmag to set a magnetization direction for this layer before you set a gyrotropy")
        self.Layercontent.cg=x

        if(self.Layercontent.type==4):
            if( self.Layercontent.magdir== 1):
                self.Layercontent.cyz=-x
                self.Layercontent.czy=x
            elif( self.Layercontent.magdir== 2):
                self.Layercontent.cxz=x
                self.Layercontent.czx=-x
            elif( self.Layercontent.magdir== 3):
                self.Layercontent.cxy=x
                self.Layercontent.cyx=-x
        else:
            self.Layercontent.type=3




cdef void copyR( ndarray[double, ndim=1, mode="c"] R, double *ptr,  int N ):
    memcpy( &R[0] , ptr, N*sizeof( double )  )

cdef void copyRcomp( ndarray[double complex, ndim=1, mode="c"] R, double complex *ptr,  int N ):
    memcpy( &R[0] , ptr, N*sizeof( double complex )  )



def Generate_structure(int NLayers_types, MLstructure="default"):




#    cdef Heterostructure STR
    if(NLayers_types<=0):
        raise Exception("Please generate at least one layer!")

    cdef int i
    HS=[0 for i in range(NLayers_types)]
    HS[0]=Layer(0,0,MLstructure, NLayers_types, 1)
    for i in range(1,NLayers_types):
        HS[i]=Layer(0,0,0,0,0)

  #  HS[0].initstructure(MLstructure, NLayers_types)

    return HS




def Reflectivity(HS, th, wavelength, MultipleScattering=1, MagneticCutoff=1.0e-50, Output="T"):


    if( HS[0].isthisthelowestlayer()!=1   ):
        raise Exception("Underlying structure not initialized. Please generate the layer list with Generate_Structure!")




    if(hasattr(th, "__len__")==False):
        th=ndarray([th])
    cdef int NAngles=len(th)

    if((th[0]<=0) or (th[NAngles-1]>90)):
        raise Exception("Theta must be in the range 0<theta<=90")

    cdef long a=(HS[0].motherpointer())
    cdef Heterostructure* A=<Heterostructure*>a #Recover the C storage of the structure
    cdef CLayer* B

    cdef int NLayers=(A[0]).NLayers

    cdef int NLayers_types=(A[0]).NLayers_types
    cdef int i


    cdef int allx=0
    cdef int ally=0
    cdef int allz=0
    cdef int UseFullMatrix=0
    cdef int Number_of_matrix_type_layers=0

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
     #   print("Layer ", B[0].type, B[0].cg)
        if(B[0].type==3):
            Setting3=1
            if(B[0].magdir==1):
                allx=1
            elif(B[0].magdir==2):
                ally=1
            elif(B[0].magdir==3):
                allz=1

            if(Cmaxnorm(B[0].cg)<MagneticCutoff):#Apply the magnetic cutoff
               # print("apply cutoff to layer", i)
               # HS[i].seteps([HS[i].epsxx(), HS[i].epsyy(), HS[i].epszz()])

                ((A[0]).LR[i]).cg=0
                ((A[0]).LR[i]).type=2
                ((A[0]).LR[i]).magdir=0

        if(B[0].type==4):
            Reduce_complexity_of_chi( &((A[0]).LR[i]), MagneticCutoff , &allx, &ally, &allz)
        if( ((A[0]).LR[i]).type==4 ):
            UseFullMatrix=1
            Number_of_matrix_type_layers+=1

        if( ((A[0]).LR[i]).type<3 ):
            ((A[0]).LR[i]).magdir=0

  #  print("hallo", allx, ally, allz)



    if(Setting1!=3  and ( ally or allz or UseFullMatrix ) ):

        raise Exception('Exception! Magnetic heterostructures must have "t" or "T" as an output parameter')

    #print("debug mode")
   # UseFullMatrix=1


    if((ally&allz)|(ally&allx)|(allx&allz)):
        #raise Exception('Exception! Multiple magnetization directions are so far not supported!')
        UseFullMatrix=1
    elif(allx or ally or allz):
        Setting3=1
    cdef MatrixSafer *AllMS
    cdef int *Layer_type_to_Matrixsafe

    if(UseFullMatrix):
      #  raise Exception('Exception! Multiple magnetization directions are so far not supported!')
        AllMS=<MatrixSafer*>malloc( Number_of_matrix_type_layers*sizeof(MatrixSafer))
        Layer_type_to_Matrixsafe = <int*>malloc( NLayers_types*sizeof(int) )
        j=0
        for i in range( NLayers_types ):
            if( ((A[0]).LR[i]).type==4 ):
                Fill_Matrixsafer( &(AllMS[j]), (A[0]).LR[i] )
                Layer_type_to_Matrixsafe[i]=j
                j+=1
   # print("filled all saves", Number_of_matrix_type_layers )

#    print("Matrix-Safer")
#    print("there are", Number_of_matrix_type_layers, "matrixlayers" )
#    for i in range(NLayers_types):
#        print("Layer", i, "safed in", Layer_type_to_Matrixsafe[i] )




    cdef double complex rss, rpp
    cdef double complex rmat[2][2]

    if(UseFullMatrix):
    #    print("use full")
        if( Setting2 ):
            Result=zeros( (4, NAngles) )
        else:
            Result=zeros( (4, NAngles) , dtype='complex')
        if(Setting2): #means absquadr

            rpiarr     = <double*>malloc(NAngles*sizeof(double))
            rsigmaarr     = <double*>malloc(NAngles*sizeof(double))
            rleftarr     = <double*>malloc(NAngles*sizeof(double))
            rrightarr     = <double*>malloc(NAngles*sizeof(double))
            for i in range(NAngles):
                Full_Matrix(A, AllMS, Layer_type_to_Matrixsafe, th[i], wavelength, &rmat)
                rsigmaarr[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[0][1])
                rpiarr[i]=cabsquadr(rmat[1][1])+cabsquadr(rmat[1][0])
                rleftarr[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
                rrightarr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
            copyR( Result[ 0 ], &rsigmaarr[ 0 ], NAngles )
            copyR( Result[ 1 ], &rpiarr[ 0 ], NAngles )
            copyR( Result[ 2 ], &rleftarr[ 0 ], NAngles )
            copyR( Result[ 3 ], &rrightarr[ 0 ], NAngles )
            free(rpiarr)
            free(rsigmaarr)
            free(rleftarr)
            free(rrightarr)

        else:
            r11 = <double complex*>malloc(NAngles*sizeof(double complex))
            r22 = <double complex*>malloc(NAngles*sizeof(double complex))
            r12 = <double complex*>malloc(NAngles*sizeof(double complex))
            r21 = <double complex*>malloc(NAngles*sizeof(double complex))
            for i in range(NAngles):
                Full_Matrix(A, AllMS, Layer_type_to_Matrixsafe, th[i], wavelength, &rmat)
                r11[i]=rmat[0][0]
                r12[i]=rmat[0][1]
                r21[i]=rmat[1][0]
                r22[i]=rmat[1][1]
            copyRcomp( Result[3], &r22[ 0 ], NAngles )
            copyRcomp( Result[0], &r11[ 0 ], NAngles )
            copyRcomp( Result[1], &r12[ 0 ], NAngles )
            copyRcomp( Result[2], &r21[ 0 ], NAngles )
            free(r11)
            free(r22)
            free(r12)
            free(r21)

    elif(Setting3==0):#Not magnetic

        if( Setting2 ):
            Result=zeros( (2, NAngles) )
        else:
            Result=zeros( (2, NAngles) , dtype='complex')

        if(Setting1!=2): #2 means only pi


            if(MultipleScattering):

                if(Setting2): #means absquadr
                    routs=zeros(NAngles)
                    rsigmaarr  = <double*>malloc(NAngles*sizeof(double))
                    for i in range(NAngles):
                        #rss=LinDicParatt_Sigma_MS(A, th[i], wavelength)
                        #routs[i]=cabsquadr(rss)
                        rsigmaarr[i]=cabsquadr( LinDicParatt_Sigma_MS(A, th[i], wavelength) )
                    copyR( Result[0], &rsigmaarr[ 0 ], NAngles )
                    free(rsigmaarr)
                else:
                    #routs=zeros(NAngles, dtype=complex)
                    r11 = <double complex*>malloc(NAngles*sizeof(double complex))
                    for i in range(NAngles):
                        #rss=
                        #routs[i]=rss
                        r11[i] = LinDicParatt_Sigma_MS(A, th[i], wavelength)
                    copyRcomp( Result[0], &r11[ 0 ], NAngles )
                    free(r11)
            else:

                if(Setting2): #means absquadr
                    routs=zeros(NAngles)
                    rsigmaarr  = <double*>malloc(NAngles*sizeof(double))
                    for i in range(NAngles):
                        #rss=LinDicParatt_Sigma(A, th[i], wavelength)
                        #routs[i]=cabsquadr(rss)
                        rsigmaarr[i]=cabsquadr( LinDicParatt_Sigma(A, th[i], wavelength) )
                    copyR( Result[0], &rsigmaarr[ 0 ], NAngles )
                    free(rsigmaarr)
                else:
                    #routs=zeros(NAngles, dtype=complex)
                    r11 = <double complex*>malloc(NAngles*sizeof(double complex))
                    for i in range(NAngles):
                        #rss=
                        #routs[i]=rss
                        r11[i] = LinDicParatt_Sigma(A, th[i], wavelength)
                    copyRcomp( Result[0], &r11[ 0 ], NAngles )
                    free(r11)
        if(Setting1!=1): #1 means only sigma

            if(MultipleScattering):
                if(Setting2): #means absquadr
                    #routp=zeros(NAngles)
                    rpiarr     = <double*>malloc(NAngles*sizeof(double))
                    for i in range(NAngles):
                        #rpp=LinDicParatt_Pi_MS(A, th[i], wavelength)
                        #routp[i]=cabsquadr(rpp)
                        rpiarr[i]=cabsquadr( LinDicParatt_Pi_MS(A, th[i], wavelength) )
                    copyR( Result[ 1 ], &rpiarr[ 0 ], NAngles )
                    free(rpiarr)

                else:
                    #routp=zeros(NAngles, dtype=complex)
                    r22 = <double complex*>malloc(NAngles*sizeof(double complex))
                    for i in range(NAngles):
                        #rpp=LinDicParatt_Pi_MS(A, th[i], wavelength)
                        #routp[i]=rpp
                        r22[i] = LinDicParatt_Pi_MS(A, th[i], wavelength)
                    copyRcomp( Result[1], &r22[ 0 ], NAngles )
                    free(r22)
            else:
                if(Setting2): #means absquadr
                    #routp=zeros(NAngles)
                    rpiarr     = <double*>malloc(NAngles*sizeof(double))
                    for i in range(NAngles):
                       # rpp=LinDicParatt_Pi(A, th[i], wavelength)
                       # routp[i]=cabsquadr(rpp)
                        rpiarr[i]=cabsquadr( LinDicParatt_Pi(A, th[i], wavelength) )
                    copyR( Result[ 1 ], &rpiarr[ 0 ], NAngles )
                    free(rpiarr)
                else:
                    #routp=zeros(NAngles, dtype=complex)
                    r22 = <double complex*>malloc(NAngles*sizeof(double complex))
                    for i in range(NAngles):
                        #rpp=LinDicParatt_Pi(A, th[i], wavelength)
                        #routp[i]=rpp
                        r22[i] = LinDicParatt_Pi(A, th[i], wavelength)
                    copyRcomp( Result[1], &r22[ 0 ], NAngles )
                    free(r22)
    else:#Magnetic
    #    print("use mag")
        if(allx):
            if( Setting2 ):
                Result=zeros( (2, NAngles) )
            else:
                Result=zeros( (2, NAngles) , dtype='complex')
        else:
            if( Setting2 ):
                Result=zeros( (4, NAngles) )
            else:
                Result=zeros( (4, NAngles) , dtype='complex')

        if(allx):
            if(Setting1!=2): #2 means only pi


                if(MultipleScattering):

                    if(Setting2): #means absquadr
                        routs=zeros(NAngles)
                        rsigmaarr  = <double*>malloc(NAngles*sizeof(double))
                        for i in range(NAngles):
                            #rss=LinDicParatt_Sigma_MS(A, th[i], wavelength)
                            #routs[i]=cabsquadr(rss)
                            rsigmaarr[i]=cabsquadr( LinDicParatt_Sigma_MS(A, th[i], wavelength) )
                        copyR( Result[0], &rsigmaarr[ 0 ], NAngles )
                        free(rsigmaarr)
                    else:
                        #routs=zeros(NAngles, dtype=complex)
                        r11 = <double complex*>malloc(NAngles*sizeof(double complex))
                        for i in range(NAngles):
                            #rss=
                            #routs[i]=rss
                            r11[i] = LinDicParatt_Sigma_MS(A, th[i], wavelength)
                        copyRcomp( Result[0], &r11[ 0 ], NAngles )
                        free(r11)
                else:

                    if(Setting2): #means absquadr
                        routs=zeros(NAngles)
                        rsigmaarr  = <double*>malloc(NAngles*sizeof(double))
                        for i in range(NAngles):
                            #rss=LinDicParatt_Sigma(A, th[i], wavelength)
                            #routs[i]=cabsquadr(rss)
                            rsigmaarr[i]=cabsquadr( LinDicParatt_Sigma(A, th[i], wavelength) )
                        copyR( Result[0], &rsigmaarr[ 0 ], NAngles )
                        free(rsigmaarr)
                    else:
                        #routs=zeros(NAngles, dtype=complex)
                        r11 = <double complex*>malloc(NAngles*sizeof(double complex))
                        for i in range(NAngles):
                            #rss=
                            #routs[i]=rss
                            r11[i] = LinDicParatt_Sigma(A, th[i], wavelength)
                        copyRcomp( Result[0], &r11[ 0 ], NAngles )
                        free(r11)
            if(Setting1!=1): #1 means only sigma

                if(MultipleScattering):

                    if(Setting2): #means absquadr
                        #routp=zeros(NAngles)
                        rpiarr     = <double*>malloc(NAngles*sizeof(double))
                        for i in range(NAngles):
                            #rpp=LinDicParatt_Pi_MS(A, th[i], wavelength)
                            #routp[i]=cabsquadr(rpp)
                            rpiarr[i]=cabsquadr( LinDicParatt_Pi_xmag_MS(A, th[i], wavelength) )
                        copyR( Result[ 1 ], &rpiarr[ 0 ], NAngles )
                        free(rpiarr)

                    else:
                        #routp=zeros(NAngles, dtype=complex)
                        r22 = <double complex*>malloc(NAngles*sizeof(double complex))
                        for i in range(NAngles):
                            #rpp=LinDicParatt_Pi_MS(A, th[i], wavelength)
                            #routp[i]=rpp
                            r22[i] = LinDicParatt_Pi_xmag_MS(A, th[i], wavelength)
                        copyRcomp( Result[1], &r22[ 0 ], NAngles )
                        free(r22)
                else:
                    if(Setting2): #means absquadr
                        #routp=zeros(NAngles)
                        rpiarr     = <double*>malloc(NAngles*sizeof(double))
                        for i in range(NAngles):
                           # rpp=LinDicParatt_Pi(A, th[i], wavelength)
                           # routp[i]=cabsquadr(rpp)
                            rpiarr[i]=cabsquadr( LinDicParatt_Pi_xmag(A, th[i], wavelength) )
                        copyR( Result[ 1 ], &rpiarr[ 0 ], NAngles )
                        free(rpiarr)
                    else:
                        #routp=zeros(NAngles, dtype=complex)
                        r22 = <double complex*>malloc(NAngles*sizeof(double complex))
                        for i in range(NAngles):
                            #rpp=LinDicParatt_Pi(A, th[i], wavelength)
                            #routp[i]=rpp
                            r22[i] = LinDicParatt_Pi_xmag(A, th[i], wavelength)
                        copyRcomp( Result[1], &r22[ 0 ], NAngles )
                        free(r22)


            ##################################################################################


        elif(ally):
          #  print("ally")

            if(MultipleScattering):
                if(Setting2): #means absquadr


#                    routs=zeros(NAngles)
#                    routp=zeros(NAngles)
#                    routl=zeros(NAngles)
#                    routr=zeros(NAngles)
#                        routs[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[0][1])
#                        routp[i]=cabsquadr(rmat[1][0])+cabsquadr(rmat[1][1])
#                        routl[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
#                        routr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                    rpiarr     = <double*>malloc(NAngles*sizeof(double))
                    rsigmaarr     = <double*>malloc(NAngles*sizeof(double))
                    rleftarr     = <double*>malloc(NAngles*sizeof(double))
                    rrightarr     = <double*>malloc(NAngles*sizeof(double))
                    for i in range(NAngles):
                        Paratt_magnetic_y_MS(A, th[i], wavelength, &rmat)
                        rsigmaarr[i]=cabsquadr(rmat[0][0])++cabsquadr(rmat[0][1])
                        rpiarr[i]=cabsquadr(rmat[1][1])+cabsquadr(rmat[1][0])
                        rleftarr[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
                        rrightarr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                    copyR( Result[ 0 ], &rsigmaarr[ 0 ], NAngles )
                    copyR( Result[ 1 ], &rpiarr[ 0 ], NAngles )
                    copyR( Result[ 2 ], &rleftarr[ 0 ], NAngles )
                    copyR( Result[ 3 ], &rrightarr[ 0 ], NAngles )
                    free(rpiarr)
                    free(rsigmaarr)
                    free(rleftarr)
                    free(rrightarr)

                else:
#                    routs=zeros(NAngles, dtype=complex)
#                    routp=zeros(NAngles, dtype=complex)
#                    routl=zeros(NAngles, dtype=complex)
#                    routr=zeros(NAngles, dtype=complex)
#                    for i in range(NAngles):
#                        Paratt_magnetic_y_MS(A, th[i], wavelength, &rmat)
#                        routs[i]=rmat[0][0]
#                        routp[i]=rmat[0][1]
#                        routl[i]=rmat[1][0]
#                        routr[i]=rmat[1][1]
                    r11 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r22 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r12 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r21 = <double complex*>malloc(NAngles*sizeof(double complex))
                    for i in range(NAngles):
                        Paratt_magnetic_y_MS(A, th[i], wavelength, &rmat)
                        r11[i]=rmat[0][0]
                        r12[i]=rmat[0][1]
                        r21[i]=rmat[1][0]
                        r22[i]=rmat[1][1]
                    copyRcomp( Result[3], &r22[ 0 ], NAngles )
                    copyRcomp( Result[0], &r11[ 0 ], NAngles )
                    copyRcomp( Result[1], &r12[ 0 ], NAngles )
                    copyRcomp( Result[2], &r21[ 0 ], NAngles )
                    free(r11)
                    free(r22)
                    free(r12)
                    free(r21)
            else:

                if(Setting2): #means absquadr
#                    routs=zeros(NAngles)
#                    routp=zeros(NAngles)
#                    routl=zeros(NAngles)
#                    routr=zeros(NAngles)
#                    for i in range(NAngles):
#                        Paratt_magnetic_y(A, th[i], wavelength, &rmat)
#                        routs[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[1][0])
#                        routp[i]=cabsquadr(rmat[1][0])+cabsquadr(rmat[1][1])
#                        routl[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
#                        routr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                    rpiarr     = <double*>malloc(NAngles*sizeof(double))
                    rsigmaarr     = <double*>malloc(NAngles*sizeof(double))
                    rleftarr     = <double*>malloc(NAngles*sizeof(double))
                    rrightarr     = <double*>malloc(NAngles*sizeof(double))
                    for i in range(NAngles):
                        Paratt_magnetic_y(A, th[i], wavelength, &rmat)
                        rsigmaarr[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[0][1])
                        rpiarr[i]=cabsquadr(rmat[1][0])+cabsquadr(rmat[1][1])
                        rleftarr[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
                        rrightarr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                    copyR( Result[ 0 ], &rsigmaarr[ 0 ], NAngles )
                    copyR( Result[ 1 ], &rpiarr[ 0 ], NAngles )
                    copyR( Result[ 2 ], &rleftarr[ 0 ], NAngles )
                    copyR( Result[ 3 ], &rrightarr[ 0 ], NAngles )
                    free(rpiarr)
                    free(rsigmaarr)
                    free(rleftarr)
                    free(rrightarr)

                else:
#                    routs=zeros(NAngles, dtype=complex)
#                    routp=zeros(NAngles, dtype=complex)
#                    routl=zeros(NAngles, dtype=complex)
#                    routr=zeros(NAngles, dtype=complex)
#                    for i in range(NAngles):
#                        Paratt_magnetic_y(A, th[i], wavelength, &rmat)
#                        routs[i]=rmat[0][0]
#                        routp[i]=rmat[0][1]
#                        routl[i]=rmat[1][0]
#                        routr[i]=rmat[1][1]
                    r11 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r22 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r12 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r21 = <double complex*>malloc(NAngles*sizeof(double complex))
                    for i in range(NAngles):
                        Paratt_magnetic_y(A, th[i], wavelength, &rmat)
                        r11[i]=rmat[0][0]
                        r12[i]=rmat[0][1]
                        r21[i]=rmat[1][0]
                        r22[i]=rmat[1][1]
                    copyRcomp( Result[3], &r22[ 0 ], NAngles )
                    copyRcomp( Result[0], &r11[ 0 ], NAngles )
                    copyRcomp( Result[1], &r12[ 0 ], NAngles )
                    copyRcomp( Result[2], &r21[ 0 ], NAngles )
                    free(r11)
                    free(r22)
                    free(r12)
                    free(r21)
        elif(allz):
            if(MultipleScattering):
                if(Setting2): #means absquadr
#                    routs=zeros(NAngles)
#                    routp=zeros(NAngles)
#                    routl=zeros(NAngles)
#                    routr=zeros(NAngles)
#                    for i in range(NAngles):
#                        Paratt_magnetic_z_MS(A, th[i], wavelength, &rmat)
#                        routs[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[0][1])
#                        routp[i]=cabsquadr(rmat[1][0])+cabsquadr(rmat[1][1])
#                        routl[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
#                        routr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                    rpiarr     = <double*>malloc(NAngles*sizeof(double))
                    rsigmaarr     = <double*>malloc(NAngles*sizeof(double))
                    rleftarr     = <double*>malloc(NAngles*sizeof(double))
                    rrightarr     = <double*>malloc(NAngles*sizeof(double))
                    for i in range(NAngles):
                        Paratt_magnetic_z_MS(A, th[i], wavelength, &rmat)
                        rsigmaarr[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[0][1])
                        rpiarr[i]=cabsquadr(rmat[1][0])+cabsquadr(rmat[1][1])
                        rleftarr[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
                        rrightarr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                    copyR( Result[ 0 ], &rsigmaarr[ 0 ], NAngles )
                    copyR( Result[ 1 ], &rpiarr[ 0 ], NAngles )
                    copyR( Result[ 2 ], &rleftarr[ 0 ], NAngles )
                    copyR( Result[ 3 ], &rrightarr[ 0 ], NAngles )
                    free(rpiarr)
                    free(rsigmaarr)
                    free(rleftarr)
                    free(rrightarr)
                else:
#                    routs=zeros(NAngles, dtype=complex)
#                    routp=zeros(NAngles, dtype=complex)
#                    routl=zeros(NAngles, dtype=complex)
#                    routr=zeros(NAngles, dtype=complex)
#                    for i in range(NAngles):
#                        Paratt_magnetic_z_MS(A, th[i], wavelength, &rmat)
#                        routs[i]=rmat[0][0]
#                        routp[i]=rmat[0][1]
#                        routl[i]=rmat[1][0]
#                        routr[i]=rmat[1][1]
                    r11 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r22 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r12 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r21 = <double complex*>malloc(NAngles*sizeof(double complex))
                    for i in range(NAngles):
                        Paratt_magnetic_z_MS(A, th[i], wavelength, &rmat)
                        r11[i]=rmat[0][0]
                        r12[i]=rmat[0][1]
                        r21[i]=rmat[1][0]
                        r22[i]=rmat[1][1]
                    copyRcomp( Result[3], &r22[ 0 ], NAngles )
                    copyRcomp( Result[0], &r11[ 0 ], NAngles )
                    copyRcomp( Result[1], &r12[ 0 ], NAngles )
                    copyRcomp( Result[2], &r21[ 0 ], NAngles )
                    free(r11)
                    free(r22)
                    free(r12)
                    free(r21)
            else:
                if(Setting2): #means absquadr
#                    routs=zeros(NAngles)
#                    routp=zeros(NAngles)
#                    routl=zeros(NAngles)
#                    routr=zeros(NAngles)
#                    for i in range(NAngles):
#                        Paratt_magnetic_z(A, th[i], wavelength, &rmat)
#                        routs[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[0][1])
#                        routp[i]=cabsquadr(rmat[1][0])+cabsquadr(rmat[1][1])
#                        routl[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
#                        routr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                    rpiarr     = <double*>malloc(NAngles*sizeof(double))
                    rsigmaarr     = <double*>malloc(NAngles*sizeof(double))
                    rleftarr     = <double*>malloc(NAngles*sizeof(double))
                    rrightarr     = <double*>malloc(NAngles*sizeof(double))
                    for i in range(NAngles):
                        Paratt_magnetic_z(A, th[i], wavelength, &rmat)
                        rsigmaarr[i]=cabsquadr(rmat[0][0])+cabsquadr(rmat[0][1])
                        rpiarr[i]=cabsquadr(rmat[1][1])+cabsquadr(rmat[1][0])
                        rleftarr[i]=0.5*(cabsquadr(rmat[0][0]-1j*rmat[0][1])+cabsquadr(rmat[1][0]-1j*rmat[1][1]) )
                        rrightarr[i]=0.5*(cabsquadr(rmat[0][0]+1j*rmat[0][1])+cabsquadr(rmat[1][0]+1j*rmat[1][1]) )
                    copyR( Result[ 0 ], &rsigmaarr[ 0 ], NAngles )
                    copyR( Result[ 1 ], &rpiarr[ 0 ], NAngles )
                    copyR( Result[ 2 ], &rleftarr[ 0 ], NAngles )
                    copyR( Result[ 3 ], &rrightarr[ 0 ], NAngles )
                    free(rpiarr)
                    free(rsigmaarr)
                    free(rleftarr)
                    free(rrightarr)

                else:
#                    routs=zeros(NAngles, dtype=complex)
#                    routp=zeros(NAngles, dtype=complex)
#                    routl=zeros(NAngles, dtype=complex)
#                    routr=zeros(NAngles, dtype=complex)
#                    for i in range(NAngles):
#                        Paratt_magnetic_z(A, th[i], wavelength, &rmat)
#                        routs[i]=rmat[0][0]
#                        routp[i]=rmat[0][1]
#                        routl[i]=rmat[1][0]
#                        routr[i]=rmat[1][1]
                    r11 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r22 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r12 = <double complex*>malloc(NAngles*sizeof(double complex))
                    r21 = <double complex*>malloc(NAngles*sizeof(double complex))
                    for i in range(NAngles):
                        Paratt_magnetic_z(A, th[i], wavelength, &rmat)
                        r11[i]=rmat[0][0]
                        r12[i]=rmat[0][1]
                        r21[i]=rmat[1][0]
                        r22[i]=rmat[1][1]
                    copyRcomp( Result[3], &r22[ 0 ], NAngles )
                    copyRcomp( Result[0], &r11[ 0 ], NAngles )
                    copyRcomp( Result[1], &r12[ 0 ], NAngles )
                    copyRcomp( Result[2], &r21[ 0 ], NAngles )
                    free(r11)
                    free(r22)
                    free(r12)
                    free(r21)




    if(UseFullMatrix):
        free(AllMS)
        free( Layer_type_to_Matrixsafe )
        return Result
   # if(allx):
    #    if(Setting2):
    #        rempty=0.5*(routp+routs)
    if(Setting3==0 or allx):
        if(Setting1==1):
            return Result[0] #routs
        elif(Setting1==2):
            return Result[1] #routp
        elif(Setting1==3):
            return Result
            #return [routs, routp]

    else:
        return Result

def KramersKronig(e, im):


    L1=len(e)
    L2=len(im)
    if(L1 != L2):
        raise Exception("The Kramers-Kronig transformation failed because the input arrays have different lengths")

    re=zeros(L1)
    Kramers_Kronig_Transformation_internal(e, im, re, L1)
    return re

