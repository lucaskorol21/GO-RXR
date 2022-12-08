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


from Mathematical_Functions_Reflectivity cimport *

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

cdef void Calculate_Multilayer(double complex *t_comp1_up, double complex *t_comp2_up, double complex *t_comp1_do, double complex *t_comp2_do, double complex *r_ML_in1, double complex *r_ML_in2, double complex *r_ML_ba1, double complex *r_ML_ba2, int N):

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



cdef void Calculate_ANXBN(double complex (*A)[2][2], double complex (*B)[2][2], double complex (*X)[2][2], int N):
    cdef int expite
    cdef int i,j
    cdef double complex  resA[2][2]
    cdef double complex  resB[2][2]

    expite=N


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

        expite = expite//2
        Mult2x2_leftside(A, A);
        Mult2x2_leftside(B, B);




    for i in range(2):
        for j in range(2):
            (A[0])[i][j]=resA[i][j];
            (B[0])[i][j]=resB[i][j];


    Mult2x2_leftside(X, B)
    Mult2x2_rightside(A, X)


cdef void Matrixexp( double complex (*A)[4][4], int N):
    cdef int expite
    cdef int i,j
    cdef double complex  res[4][4]
   # cdef double norm, d

    expite=N;
    for i in range(4):
        for j in range(4):
            if(i==j):
                res[i][j]=1;
            else:
                res[i][j]=0
    expite=N
    while(expite):
        if(expite%2):
            Mult4x4_leftside( &res, A )
        expite=expite//2
        Mult4x4_leftside( A, A )
#        norm=0
#        for i in range(4):
#            for j in range(4):
#                d=Cmaxnorm( res[i][j])
#                if( d>norm ):
#                    norm=d
#        for i in range(4):
#            for j in range(4):
#                res[i][j]/=d



    for i in range(4):
        for j in range(4):
            (A[0])[i][j]=res[i][j]




cdef void Calculate_Multilayer_equation(double complex  (*A)[2][2], double complex  (*B)[2][2], double complex (*X)[2][2], double complex  (*result)[2][2], int N):
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


cdef void Calculate_Multilayer_with_Matrices(double complex (*t_comp1_up)[2][2], double complex (*t_comp2_up)[2][2], double complex (*t_comp1_do)[2][2], double complex (*t_comp2_do)[2][2], double complex (*r_ML_in1)[2][2], double complex (*r_ML_in2)[2][2], double complex (*r_ML_ba1)[2][2], double complex (*r_ML_ba2)[2][2], int N):

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
