
### Author : Kevin Jiguet-Covex


import sys,getopt,os.path
import crystallography as cr
import umd_process as umdp
import time
import numpy as np
import concurrent.futures
from functools import partial
import ctypes
from os.path import join


current_path=os.path.abspath(__file__)#For all this to work, the file c_msd.so must be in the same directory than this script
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u
msd_lib = ctypes.cdll.LoadLibrary(join(path_new, 'c_msd.so'))
msd_lib.compute_msd.argtypes = [ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double), ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
msd_lib.compute_msd.restype = None



def msdAtom_C(n,pos,hh,vv,ballistic):
    PosAr=np.array(pos)
    nitmax=(len(pos)//2-ballistic)
    msd=np.array([0.0 for _ in range(nitmax//vv)])    
    Posp = PosAr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    MSDp = msd.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    msd_lib.compute_msd(Posp,MSDp,hh,vv,ballistic,nitmax)
    return msd



def main(argv):
    hh = 1
    vv = 1
    ballistic = 0
    TimeStep = 1
    umdp.headerumd()
    start=time.time()
    umdfile=''
    try:
        opts, arg = getopt.getopt(argv,"hf:z:v:b:x:",["fumdfile","zHorizontalJump","vVerticalJump","bBallistic","xAtoms"])
    except getopt.GetoptError:
        print ('msd_umd.py -f <XYZ_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic> -a <Atoms>')
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
            print ('Atoms = atoms or elements ; parameter to either print the msd of each individual atom (+ the mean msd of each element) or only the msd of the elements.')
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
        elif opt in ("-x","--xAtoms"):
            x=str(arg)
            
    if (os.path.isfile(umdfile)):
        
        MyCrystal=cr.Lattice()
        ff=open(umdfile,'r')
        while True :
            line=ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'natom':
                    MyCrystal.natom = int(entry[1])
                    MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
                    dicoAtoms={at:[] for at in range(MyCrystal.natom)}

                if entry[0] == 'ntypat':
                    MyCrystal.ntypat = int(entry[1])
                    MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                    MyCrystal.masses = [1.0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.zelec = [1.0 for _ in range(MyCrystal.ntypat)]
                if entry[0] == 'types':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.types[ii]=int(entry[ii+1])
                if entry[0] == 'elements':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.elements[ii]=entry[ii+1]
                if entry[0] == 'masses':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.masses[ii]=float(entry[ii+1])
                if entry[0] == 'Zelectrons':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.zelec[ii]=float(entry[ii+1])
                if entry[0] == 'typat':
                    for ii in range(MyCrystal.natom):
                        MyCrystal.typat[ii] = int(entry[ii+1])
                if entry[0] == 'timestep':
                    TimeStep = float(entry[1])
                if entry[0] == 'atoms:':
                    for iatom in range(MyCrystal.natom):
                        line=ff.readline()
                        line=line.strip()
                        entry=line.split()
                        dicoAtoms[iatom].append([float(entry[6]),float(entry[7]),float(entry[8])])
                        
        print ('Number of atoms of each type is ',MyCrystal.types)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            msdAtomRed=partial(msdAtom_C,hh=hh,vv=vv,ballistic=ballistic)
            msdfile = umdfile[:-8] + '.msd_par'
                
               
            msdArray=list(executor.map(msdAtomRed,[iatom for iatom in range(MyCrystal.natom)],[dicoAtoms[iatom] for iatom in range(MyCrystal.natom)]))
                
            if(x=="elements"):
                
                MSD=np.array([[0.0 for _ in range(len(msdArray[0]))] for _ in range(MyCrystal.ntypat)])
                for iatom in range(len(msdArray)) :
                    MSD[MyCrystal.typat[iatom]]+=msdArray[iatom]
                msdfile+='.elements.dat'
                f = open(msdfile,'w')
                string='time_(fs)\t'
                for itypat in range(MyCrystal.ntypat):
                    string=string + MyCrystal.elements[itypat] + '\t'
                string = string + '\n'
                f.write(string)
                niter=len(dicoAtoms[0])
 #                   print(niter, hh, ballistic)
                weight=int((int(niter/2)-ballistic)/hh)-1
 #                   print("weight=",weight*MyCrystal.natom)                    
                for ii in range(len(MSD[0])):     
                    instant = (float(ii)*TimeStep*vv)+ballistic
                    string = str(instant)
                    for jj in range(MyCrystal.ntypat):
                        string = string + '\t' + str(MSD[jj][ii]/(float(MyCrystal.types[jj])*float(weight)))
                    string = string + '\n'
                    f.write(string)
                print ('MSDs printed in file ',msdfile)
   
            elif(x=="atoms"):
                
                msdfile+='.atoms.dat'
                f = open(msdfile,'w')
                string='time_(fs)\t'
                niter=len(dicoAtoms[0])
                weight=int((int(niter/2)-ballistic)/hh)-1
                string='time (fs)'
                for iatom in range(MyCrystal.natom):
                    string+='\t'+MyCrystal.elements[MyCrystal.typat[iatom]]+str(iatom)
                f.write(string+'\n')
                for ii in range(len(msdArray[0])):     
                    instant = (float(ii)*TimeStep*vv)+ballistic
                    string = str(instant)
                    for iatom in range(MyCrystal.natom):
                        string = string + '\t' + str(msdArray[iatom][ii]/(float(weight)))
                    string = string + '\n'
                    f.write(string)
                print ('MSDs printed in file ',msdfile)
         
    else:
        print ('umd file ',umdfile,'does not exist')
        sys.exit()
    
    end=time.time()
    
    print("runtime :",end-start)
 

if __name__ == "__main__":
    main(sys.argv[1:])
