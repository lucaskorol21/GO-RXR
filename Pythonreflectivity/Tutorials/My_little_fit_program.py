import numpy as np
import Pythonreflectivity as pr
from scipy import linalg
from scipy.interpolate import interp1d
from joblib import Parallel, delayed
import copy
import matplotlib.pyplot as plt


hbar_times_c=1973.16 #ev times Angstrom
re= 2.8179e-5 # Thomson scattering length

contribution_of_xmcd_to_fit=20
numerical_derivative_factor=1.0e-9


def Lorentzian(a,s,t,x):
    return a*s/(s*s+(t-x)**2 )


def Fitparameters_to_Heterostructure( HS, fitparameters, energy, fitparameters_to_scriptparameters ):
    ##This function takes the fitparameters and fills the heterostructure HS with thicknesses, roughnesses and optical constants chi

    number_of_distinct_layers=len(HS)
    number_of_fitparameters=len(fitparameters)
    script=np.zeros(number_of_fitparameters+100)

    for i in range(number_of_fitparameters):
        script[ fitparameters_to_scriptparameters[i] ]=fitparameters[i]


    en=np.linspace(600, 750, 301)
    peak_height=script[0]
    peak_width=script[1]
    peak_position=script[2]
    peak_step=script[3]

    ##Model for the absorption peak
    #This is possibly a rather inefficient way of doing it because the Kramers-Kronig transformation should
    #not really by done on every energy point. But this is merely an example.
    
    ffim=Lorentzian( peak_height, peak_width, peak_position, en)+ peak_step*(np.arctan( (en-peak_position)/peak_width )+np.pi/2)
    ##Kramers-Kronig transformation
    ffre=pr.KramersKronig(en, ffim)

    ffre_offset=script[7]
    ffre=ffre+ffre_offset
    
    ff=ffre+1j*ffim

    mag_height=script[4]
    mag_width=script[5]
    mag_position=script[6]
    ffmag_im=Lorentzian( mag_height, mag_width, mag_position, en)
    ffmag_re=pr.KramersKronig(en, ffmag_im)
    ffmag=-1j*( ffmag_re+1j*ffmag_im )

    ##find the form factor  the input energy
    ff_func=interp1d(en, ff)
    ffmag_func=interp1d(en, ffmag)
    
    localff=ff_func(energy)
    localffmag=ffmag_func(energy)


    k=hbar_times_c/energy #vacuum k
    
    atom_density=0.2 # moles per cubic centimetre

    
    
    chi_diag = localff*atom_density*re*4*np.pi/(k*k)
    chi_mag  = localffmag*atom_density*re*4*np.pi/(k*k)
    
    chi1=np.array([chi_diag, chi_diag, chi_diag, chi_mag])
    HS[1].setmag("y")
    HS[1].setchi(chi1)
    
    chi_substrate=-0.0005+0.0005j
    HS[0].setchi(chi_substrate)
    
    chi_cap=(-0.0002+0.0002j)*script[8]

    HS[2].setchi(chi_cap)

    Layer_thickness=script[9]
    HS[1].setd(Layer_thickness)
    Cap_thickness=script[10]
    
    HS[2].setd(Cap_thickness)

    Substrate_roughness=script[11]
    Layer_roughness=script[12]
    Cap_roughness=script[13]

    
    
    HS[0].setsigma(Substrate_roughness)
    HS[1].setsigma(Layer_roughness)
    HS[2].setsigma(Cap_roughness)
    






def Read_fitparameters_input_file(filename):
    #This function reads the fitparameters for the start of the iteration and also the upper and lower limits in the fit
    print(("Reading Fitparameters from file " + filename ))
    #Find number of Fitparameters
    with open(filename) as f:
        j=0
        for line in f:
            j+=1

    number_of_fitparameters=j

    #Make arrays to store the input:

    upper_limits=np.zeros(number_of_fitparameters)
    lower_limits=np.zeros(number_of_fitparameters)
    fitparameters=np.zeros(number_of_fitparameters)
    fitparameters_to_scriptparameters=np.zeros(number_of_fitparameters, dtype=int)

    with open(filename) as f:
        j=0
        for line in f:
            #Read and evaluate the input file
            line=line.split(" ")
            if(line[0]=="script"):
                fitparameters_to_scriptparameters[j]=line[1]
                fitparameters[j]=float(line[2])
                lower_limits[j]=float(line[3])
                upper_limits[j]=float(line[4])
                if(upper_limits[j]<=lower_limits[j]):
                    raise Exception("Error! Upper limits<lower limits is not reasonable input!")
            else:
                raise Exception("Error while reading input file for fitparameters. script expected, found " + str(line[0]) + "instead" )
            j+=1
    #Return everything:
    return fitparameters, upper_limits, lower_limits,fitparameters_to_scriptparameters


def Read_data():
    #Reads the experimental data for the X-ray reflectivity fit
    print("Reading experimental data")
    
    files=([\
        'Data630.0.dat', \
        'Data632.0.dat', \
        'Data634.0.dat', \
        'Data636.0.dat', \
        'Data638.0.dat', \
        'Data640.0.dat', \
        'Data642.0.dat', \
        'Data644.0.dat', \
        'Data646.0.dat', \
        'Data648.0.dat', \
        'Data650.0.dat', \
        'Data652.0.dat', \
        'Data654.0.dat', \
        'Data656.0.dat', \
        'Data658.0.dat', \
        'Data660.0.dat', \
        'Data662.0.dat', \
        'Data664.0.dat', \
        'Data666.0.dat', \
        'Data668.0.dat', \
        'Data670.0.dat', \
        ])
    en=np.linspace(630,670,21)

    #store the data here
    Nfiles=len(files)
    allth=[0]*Nfiles
    allrleft=[0]*Nfiles
    allrright=[0]*Nfiles
    allxmcd=[0]*Nfiles
    
    total_number_of_data_points=0
    Nangles=np.zeros(Nfiles, dtype='int')
    for i in range(Nfiles):
        #Find the file length
        
        with open( files[i] ) as f:
            j=0
            for line in f:
                j+=1
        Nangles[i]=j
        
        total_number_of_data_points+=3*Nangles[i]
        
        ## store individual file data here
        allth[i]=np.zeros(Nangles[i])
        allrleft[i]=np.zeros(Nangles[i])
        allrright[i]=np.zeros(Nangles[i])
        allxmcd[i]=np.zeros(Nangles[i])
        
        with open( files[i] ) as f:
            j=0
            for line in f:
                #Read the data files
                line=line.split(" ")
                allth[i][j]=float( line[0] )
                allrleft[i][j]=float( line[1] )
                allrright[i][j]=float( line[2] )
                allxmcd[i][j]=float( line[3] )
                j+=1
    ##store all in one large array. This keeps the fit procedure simple later
    Rall=np.zeros(total_number_of_data_points)
    m=0
    for i in range(Nfiles):
        Rall[m:(m+Nangles[i])]=np.log( allrleft[i] )
        m+=Nangles[i]
        Rall[m:(m+Nangles[i])]=np.log( allrright[i] )
        m+=Nangles[i]
        Rall[m:(m+Nangles[i])]=contribution_of_xmcd_to_fit*allxmcd[i]
        m+=Nangles[i]

    return en, total_number_of_data_points, allth, allrleft, allrright, allxmcd, Rall

def Fitparameters_to_Reflectivity(en, allth,   fitparameters, fitparameters_to_scriptparameters, total_number_of_data_points, outputtype):
    ##This function takes the fitparameters and calculates the reflectivity for them
    number_of_fitparameters=len(fitparameters)
    number_of_energies=len(en)

    ##Find multiplier and background fit parameters
    script=np.zeros(number_of_fitparameters+100)

    for i in range(number_of_fitparameters):

        script[ fitparameters_to_scriptparameters[i] ]=fitparameters[i]
    
    background=script[20]
    multiplier=script[21]

    ##Fill this heterostructure and reflectivity array in each step
    
    HS=pr.Generate_structure(3)
    Rout=np.zeros(total_number_of_data_points)
    allR=[0]*number_of_energies
    m=0
    for i in range(number_of_energies):
        Fitparameters_to_Heterostructure( HS, fitparameters, en[i], fitparameters_to_scriptparameters )
        wavelength=2*np.pi*hbar_times_c/en[i]

        allR[i]=pr.Reflectivity( HS, allth[i], wavelength )
        allR[i][2]=multiplier*(allR[i][2]+background)
        allR[i][3]=multiplier*(allR[i][3]+background)


        Nangles=len(allth[i])
        Rout[m:(m+Nangles)]=np.log( allR[i][2] )
        m+=Nangles
        Rout[m:(m+Nangles)]=np.log( allR[i][3] )
        m+=Nangles
        Rout[m:(m+Nangles)]=contribution_of_xmcd_to_fit*(allR[i][2]-allR[i][3])/(allR[i][2]+allR[i][3])
        m+=Nangles
    if(outputtype==1):
        return Rout
    else:
        return allR

def Levenberg_Marquardt_Fitter(startfitparameters, en, allth, expdata, fitparameters_to_scriptparameters, upper_limits, lower_limits, total_number_of_data_points):

    number_of_fitparameters=len(startfitparameters)
    ite=0
    aite=startfitparameters


    #this Matrix stores all the reflectivity derivatives
    DT=np.zeros((number_of_fitparameters, total_number_of_data_points))

    ##this Matrix stores fitparameters for each point that is calculated in parallel
    AllFitparameters=[copy.copy(aite) for i in range(number_of_fitparameters+1)]

    Ncores=20
    Loop=20 #This should be something like the number of threads that can run in parallel/number of cores
    #The algorithm will first find a direction for a good descent and then check this number of points on the line.
    #The best one will yield the new fit parameter set
    
    while(True):

        if(ite>0):
            #Go here once you calculated the first step

            #Calculate the fit parameters on the line of descent
            AllFitparameters=[copy.copy(aite) for i in range(Loop)]
            for i in range(Loop):
                scale=0.5**i
                AllFitparameters[i]=AllFitparameters[i]-scale*DTr - scale*(1-scale)*DTr2

            #[aite-(0.5**i)*DTr for i in range(Loop)]

            for j in range(Loop):
                for i in range(number_of_fitparameters):
                    if(AllFitparameters[j][i]<lower_limits[i]):
                        AllFitparameters[j][i]=lower_limits[i]
                    elif(AllFitparameters[j][i]>upper_limits[i]):
                        AllFitparameters[j][i]=upper_limits[i]

          ##  Calculate the reflectivity on the line of descent
            out=Parallel(n_jobs=Ncores)(delayed(Fitparameters_to_Reflectivity)( en, allth, AllFitparameters[i], fitparameters_to_scriptparameters, total_number_of_data_points, 1)  for i in range(Loop) )




            chisquared=np.zeros(Loop)
            for i in range(Loop):
                chisquared[i] = sum( (out[i] - expdata)**2 )
            #Find the lowest chisquare
            minchi=chisquared[0]
            min_i=0
            for i in range(1, Loop):
                if(chisquared[i]<minchi):
                    minchi=chisquared[i]
                    min_i=i
            aite=copy.copy(AllFitparameters[min_i])
            
            Fiterror2=chisquared[min_i]

            if( abs( (Fiterror1-Fiterror2)/(Fiterror1+Fiterror2) ) < 1.0e-7 ):
                print(( "Converged at " + str(Fiterror2) ))
                print( "Fit parameters in convergence" )
                for i in range(number_of_fitparameters):
                    print(( "script " + str(fitparameters_to_scriptparameters[i]) +  " " + str(aite[i]) + " " + str(lower_limits[i]) + " " + str(upper_limits[i]) ))
                res=Fitparameters_to_Reflectivity( en, allth, AllFitparameters[i], fitparameters_to_scriptparameters, total_number_of_data_points, 2 )
                return aite, res
            

        #Remember how many iterations have been done
        if(ite>=1):
            print(( "Iteration ", ite, " X^2 old ", Fiterror1, " X^2 new ", Fiterror2 ))
        
        ite+=1
        
        
        #Make Fit parameters for the calculation of the derivative
        AllFitparameters=[copy.copy(aite) for i in range(number_of_fitparameters+1)]
        delta=np.zeros(number_of_fitparameters)
        for i in range(number_of_fitparameters):
            if(AllFitparameters[i][i]==0 or (lower_limits[i]<0<upper_limits[i]) ):
                delta[i]=numerical_derivative_factor*max( abs(upper_limits[i]), abs(lower_limits[i] ) )
            else:
                delta[i]=numerical_derivative_factor*AllFitparameters[number_of_fitparameters][i]
            AllFitparameters[i][i]+=delta[i]


        ##Calculate the reflectivity at each delta-step, in parallel

        out=Parallel(n_jobs=Ncores)(delayed(Fitparameters_to_Reflectivity)( en, allth, AllFitparameters[i], fitparameters_to_scriptparameters, total_number_of_data_points, 1 ) for i in range(number_of_fitparameters+1) )



        #Calculate the first fit error
        if(ite==1):
            Fiterror1=sum( (out[ number_of_fitparameters ] - expdata)**2 )
        else:
            Fiterror1=Fiterror2
        #Calculate the derivative of the reflectivity
        for i in range(number_of_fitparameters):
            DT[i] = ( out[i]-out[number_of_fitparameters] )/delta[i]

        ##Calculate the gradient
        
        out[number_of_fitparameters]-=expdata
        A=np.dot(DT,DT.T)
        b=np.dot(DT,(out[number_of_fitparameters]).T )

        #If one of the derivatives is entirely zero, the fit parameter is essentially meaningless. That may happen for a number of reasons. However, Gauss-Newton fails for this case.
        for i in range(number_of_fitparameters):
            if(b[i]==0):
                print(("WARNING! No gradient component", i, "! Singular matrix!"))
                print(b)
                return

        #Solve this system of equations to calculate the descent vector that is used for large steps
        DTr=linalg.solve(A,b,sym_pos=True)
        #This is another good descent vector that is used for small steps
        DTr2=np.zeros(number_of_fitparameters)
        for i in range(number_of_fitparameters):
            DTr2[i]=b[i]/A[i][i]


##This one here is important if you want to parallelize python. 
if __name__ == '__main__':
    fitparameters, upper_limits, lower_limits,fitparameters_to_scriptparameters=Read_fitparameters_input_file("Fitparameters.txt")
    en, total_number_of_data_points, allth, allrleft, allrright, allxmcd, Rall=Read_data()



    aite, res=Levenberg_Marquardt_Fitter(fitparameters, en, allth, Rall, fitparameters_to_scriptparameters, upper_limits, lower_limits, total_number_of_data_points)
    

    ###The "experimental data" have been calculated with this fit parameter set:
    realdataparams=np.array([261, 7, 651, 0.8, 32, 4, 649, -6, 1.2, 103, 15,2.5, 3.5, 4.1, 3.0e-7, 0.34])


    for i in range(len(en) ):
        plt.plot( allth[i], np.log( res[i][2]) )
        plt.plot( allth[i], np.log( allrleft[i] )   )
        plt.show()
        plt.plot( allth[i], (res[i][2]-res[i][3])/(res[i][3]+res[i][2]) )
        plt.plot( allth[i], allxmcd[i]  )
        plt.show()






##      
##    fitparameters=np.array([261, 7, 651, 0.8, 32, 4, 649, -6, 1.2, 103, 15,2.5, 3.5, 4.1, 3.0e-7, 0.34])
##    HS=pr.Generate_structure(3)
##
##    Nfiles=21
##    Nangles=100
##    th=np.linspace( 10, 70, Nangles)
##    Rleft=np.zeros( (Nfiles,Nangles) )
##    Rright=np.zeros( (Nfiles,Nangles) )
##    XMCD=np.zeros( (Nfiles,Nangles) )
##    endata=np.linspace(630, 670, Nfiles)
##
##
##    import matplotlib.pyplot as plt
##    i=0
##    backgr=3.0e-7
##    mult=0.34
##
##    for en in endata:
##        Fitparameters_to_Heterostructure( HS, fitparameters, en, fitparameters_to_scriptparameters )
##        wavelength=2*np.pi*hbar_times_c/en
##
##        R=pr.Reflectivity( HS, th, wavelength)
##        
##        Rleft[i]=mult*(R[2]+backgr)
##        noise1=Rleft[i]*0.01*(np.random.rand(Nangles) -0.5)
##        Rleft[i]+=noise1
##        
##        Rright[i]=mult*(R[3]+backgr)
##        noise2=Rright[i]*0.01*(np.random.rand(Nangles) -0.5)
##        Rright[i]+=noise2
##        XMCD[i]=(Rleft[i]-Rright[i])/(Rleft[i]+Rright[i])
##        filename="Data" + str(en) + ".dat"
##        f1=open(filename, "w")
##        for j in range(Nangles):
##            f1.write(str(th[j]))
##            f1.write(" ")
##            f1.write(str(Rleft[i][j]))
##            f1.write(" ")
##            f1.write(str(Rright[i][j]))
##            f1.write(" ")
##            f1.write(str(XMCD[i][j]))
##            f1.write("\n")
##        f1.close()        
##            
        i+=1
