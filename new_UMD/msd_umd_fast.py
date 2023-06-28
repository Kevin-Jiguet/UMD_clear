
### Author : Kevin Jiguet-Covex


import sys,getopt,os.path
import umd_processes_fast as umdpf
import time
import numpy as np
import concurrent.futures
from functools import partial
import ctypes
from os.path import join
import math


current_path=os.path.abspath(__file__)#For all this to work, the file c_msd.so must be in the same directory than this script
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u
msd_lib = ctypes.cdll.LoadLibrary(join(path_new, 'c_msd.so'))
msd_lib.compute_msd_tilted.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_double)]
msd_lib.compute_msd_tilted.restype = ctypes.POINTER(ctypes.c_double)

msd_lib.free_memory.argtypes = [ctypes.POINTER(ctypes.c_double)]
msd_lib.free_memory.restype = None

msd_lib.compute_msd.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
msd_lib.compute_msd.restype = ctypes.POINTER(ctypes.c_double)



def msdAtom_C(n,PosAr,numsteps,hh,vv,ballistic):
    t1=time.time()
    nitmax=(numsteps//2-ballistic)
    Posp = (PosAr.flatten()).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    msd=[]
    t2=time.time()
    msdP=msd_lib.compute_msd(Posp,hh,vv,ballistic,nitmax)
    for i in range(nitmax):
        msd.append(msdP[i])
    t3=time.time()
    msd_lib.free_memory(msdP)
    print("step "+str(n)+" : conversions "+str(t2-t1)+" s ; calculations "+str(t3-t2)+" s" )
    return msd

def msdAtom_CAxes(n,PosAr,numsteps,hh,vv,ballistic,Axes):
    t1=time.time()
    nitmax=(numsteps//2-ballistic)
    Axesp = (9 * ctypes.c_double)(*Axes)
    Posp = (PosAr.flatten()).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    msd=msd_lib.compute_msd_tilted(Posp,hh,vv,ballistic,nitmax,Axesp)
    t2=time.time()    
    msdA1,msdA2,msdA3 = [], [], []
    for i in range(int(nitmax/vv)):
        msdA1.append(msd[i])
        msdA2.append(msd[i+int(nitmax/vv)])
        msdA3.append(msd[i+2*int(nitmax/vv)])
    msd_lib.free_memory(msd)
    t3=time.time()
    print("step "+str(n)+" : calc "+str(t2-t1)+" s ; copy "+str(t3-t2)+" s" )
    return msdA1,msdA2,msdA3


def main(argv):
    hh = 1
    vv = 1
    ballistic = 0
    TimeStep = 1
    umdpf.headerumd()
    start=time.time()
    umdfile=''
    Auto = False
#    Axes=[1,0,0,0,1,0,0,0,1]
    Axes = None
    try:
        opts, arg = getopt.getopt(argv,"hf:z:v:b:m:a:",["fumdfile","zHorizontalJump","vVerticalJump","bBallistic","mMode","aAxes"])
    except getopt.GetoptError:
        print ('msd_umd.py -f <XYZ_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic> -m <Mode> -a <Axes>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Program to compute the Mean Square Displacement. Usage: ')
            print ('msd_umd_fast.py -f <UMD_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic> -x <Atoms>')
            print ('UMD_filename = input file with the trajectory, in UMD format.')
            print (' As the MSD is measured with a sliding window of various size up to half the trajectory\'s length, ')
            print (' we can accelerate the calculation using a selected reduced sampling ')
            print ('HorizontalJump = discretization for the start of the sampling window.')
            print ('VerticalJump = discretization for the length of the sampling window.')
            print ('Ballistic = estimation of the ballistic part of the trajectory. Default is 0. Typical values of 100 are sufficient.')
            print ('Mode = atoms or elements ; parameter to either print the msd of each individual atom (+ the mean msd of each element) or only the msd of the elements.')
            print ('Axes = a list of float that describes the axes along which the msd should be calculated. If this argument is absent, the default will be the traditional cartesian coordinates, and the execution will be slightly faster.\n Use < -a Auto > to use the axes as defined by the UMD lattice.')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            print ('I will use the ',umdfile,' for input','\n')
        elif opt in ("-z", "--zHorizontalJump"):
            hh = int(arg)
            if hh<=0 :
                print('ERROR : horizontal jump -z has to be strictly positive.')
                sys.exit()
        elif opt in ("-v", "--vVerticalJump"):
            vv = int(arg)
            if vv<=0 :
                print('ERROR : vertical jump -v has to be strictly positive.')                
                sys.exit()
        elif opt in ("-b", "--bBallistic"):
            ballistic = int(arg)
        elif opt in ("-m","--mMode"):
            x=str(arg)
        elif opt in ("-a","--aAxes"):
            if arg == "Auto":
                Auto =True
            else :
                Axes = eval(arg)
            
                l1 = math.sqrt(Axes[0]**2+Axes[1]**2+Axes[2]**2)
                l2 = math.sqrt(Axes[3]**2+Axes[4]**2+Axes[5]**2)
                l3 = math.sqrt(Axes[6]**2+Axes[7]**2+Axes[8]**2)
                #Normalization        
                Axes[0]/=l1 
                Axes[1]/=l1 
                Axes[2]/=l1 
                Axes[3]/=l2 
                Axes[4]/=l2 
                Axes[5]/=l2 
                Axes[6]/=l3 
                Axes[7]/=l3
                Axes[8]/=l3 

            
    if (os.path.isfile(umdfile)):
        
        MyCrystal,SnapshotsValues,TimeStep,numsteps = umdpf.read_values(umdfile,"absxcart","chunks",1,"all")
        print(len(SnapshotsValues))

        if Auto :
            Axes = MyCrystal.rprim[0] + MyCrystal.rprim[1] + MyCrystal.rprim[2]    

        
        ArrayAbsXcart = np.array(SnapshotsValues).flatten()

        ListAbsXcart = [ArrayAbsXcart[at::MyCrystal.natom] for at in range(MyCrystal.natom)]

        

        print ('Number of atoms of each type is ',MyCrystal.types)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            if Axes !=None :
                msdAtomRed=partial(msdAtom_CAxes,numsteps=numsteps,hh=hh,vv=vv,ballistic=ballistic,Axes=Axes)
            else :
                msdAtomRed=partial(msdAtom_C,numsteps=numsteps,hh=hh,vv=vv,ballistic=ballistic)

            msdfile = umdfile[:-8] + '.msd_fasterA'
                              
            msdArray=list(executor.map(msdAtomRed,[iatom for iatom in range(MyCrystal.natom)],[ListAbsXcart[iatom] for iatom in range(MyCrystal.natom)]))
            headerstring='MSD : -f '+umdfile+' -m '+ x +'-a '+str(Axes)+'\n'
            if(x=="elements"):
                
                if Axes != None :
                    msdfile+='.axes'
                    MSD=np.array([[np.zeros(len(msdArray[0][0])),np.zeros(len(msdArray[0][0])),np.zeros(len(msdArray[0][0])),np.zeros(len(msdArray[0][0]))] for _ in range(MyCrystal.ntypat)])
                    for iatom in range(len(msdArray)) :
                        for kk in range(3):
                            MSD[MyCrystal.typat[iatom]][1+kk]+=msdArray[iatom][kk]
                            MSD[MyCrystal.typat[iatom]][0]+=msdArray[iatom][kk]
                else :
                    MSD=np.array([[np.zeros(len(msdArray[0]))] for _ in range(MyCrystal.ntypat)])
                    for iatom in range(len(msdArray)) :
                        MSD[MyCrystal.typat[iatom]][0]+=msdArray[iatom]
                        
                msdfile+='.elements.dat'
                f = open(msdfile,'w')
                f.write(headerstring)
                string='time_(fs)\t'
                for itypat in range(MyCrystal.ntypat):
                    string +=  MyCrystal.elements[itypat] + '\t' + MyCrystal.elements[itypat] + ' (1)\t' + MyCrystal.elements[itypat] + '(2)\t'+ MyCrystal.elements[itypat] + '(3)\t'
                string = string + '\n'
                f.write(string)
 #                   print(niter, hh, ballistic)
                weight=int((int(numsteps/2)-ballistic)/hh)-1
 #                   print("weight=",weight*MyCrystal.natom)                    
                for ii in range(len(MSD[0][0])):     
                    instant = (float(ii)*TimeStep*vv)+ballistic
                    string = str(instant)
                    for jj in range(MyCrystal.ntypat):
                        string+='\t' + str(MSD[jj][0][ii]/(float(MyCrystal.types[jj])*float(weight)))
                        if Axes != None :
                            for kk in range(1,4):
                                string += '\t' + str(MSD[jj][kk][ii]/(float(MyCrystal.types[jj])*float(weight)))
                    string = string + '\n'
                    f.write(string)
                print ('MSDs printed in file ',msdfile)
                
            elif(x=="atoms"):
                
                if Axes != None :
                    msdfile+=".axes"
                    MSD=np.array([[np.zeros(len(msdArray[0][0])),np.zeros(len(msdArray[0][0])),np.zeros(len(msdArray[0][0])),np.zeros(len(msdArray[0][0]))] for _ in range(MyCrystal.natom)])
                    for iatom in range(len(msdArray)) :
                        for kk in range(3):
                            MSD[iatom][1+kk]=msdArray[iatom][kk]
                            MSD[iatom][0]+=msdArray[iatom][kk]
                else :
                    MSD=np.array([[np.zeros(len(msdArray[0]))] for _ in range(MyCrystal.natom)])
                    for iatom in range(len(msdArray)) :
                        MSD[iatom][0] = msdArray[iatom]
                msdfile+='.atoms.dat'
                f = open(msdfile,'w')
                f.write(headerstring)
                string='time_(fs)\t'
                lastEl=''
                for iatom in range(MyCrystal.natom):
                    if lastEl!=MyCrystal.elements[MyCrystal.typat[iatom]]:
                        indexat = 0
                        lastEl=MyCrystal.elements[MyCrystal.typat[iatom]]
                    label = lastEl+str(indexat)
                    if Axes != None :
                        string +=  label + '\t' + label + '(1)\t' + label + '(2)\t'+ label + '(3)\t'
                    else :
                        string += label +'\t'
                    
                    indexat+=1
                    
                string = string + '\n'
                f.write(string)

                weight=int((int(numsteps/2)-ballistic)/hh)-1

                for ii in range(len(MSD[0][0])):     
                    instant = (float(ii)*TimeStep*vv)+ballistic
                    string = str(instant)
                    for jj in range(MyCrystal.natom):
                        string+='\t' + str(MSD[jj][0][ii]/(float(weight)))
                        if Axes != None :
                            for kk in range(1,4):
                                string += '\t' + str(MSD[jj][kk][ii]/(float(weight)))
                    string +='\n'
                    f.write(string)


                print ('MSDs printed in file ',msdfile)
   
         
    else:
        print ('umd file ',umdfile,'does not exist')
        sys.exit()
    
    end=time.time()
    
    print("runtime :",end-start)
 

if __name__ == "__main__":
    main(sys.argv[1:])
