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

cdef int MakeConsistencyCheck(string, Heterostructure *STR, int MaxLayer):
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




