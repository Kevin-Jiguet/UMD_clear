#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:55:47 2023

@author: Kevin Jiguet-Covex
"""

import sys,getopt,os.path
import math
import time
from functools import partial
import concurrent.futures
import numpy as np
import ctypes
from os.path import join
import umd_processes_fast as umdpf
import platform

OS = platform.system()
LibraryName =""

if OS == "Linux":
    LibraryName = 'c_bonds_full.so'
elif OS == "Windows":
    LibraryName = 'c_bonds_full.dll'
elif OS == "Darwin":
    LibraryName = 'c_bonds_full.dylib'
    
current_path=os.path.abspath(__file__)#For all this to work, the file c_bonds_fullD.so must be in the same directory than gofr_umd
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u
#c_msdbtest = 

fullbond_lib = ctypes.cdll.LoadLibrary(join(path_new, LibraryName))

fullbond_lib.compute_fBonds.argtypes = [ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int]
fullbond_lib.compute_fBonds.restype = ctypes.POINTER(ctypes.c_int)

fullbond_lib.compute_fBonds_specific.argtypes = [ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.POINTER(ctypes.c_double)]
fullbond_lib.compute_fBonds_specific.restype = ctypes.POINTER(ctypes.c_int)

fullbond_lib.free_memory.argtypes = [ctypes.POINTER(ctypes.c_int)]
fullbond_lib.free_memory.restypes = None



def read_inputfile(InputFile,MyCrystal):#Creates a matrix from the .dat file containing the bond length for each pair of atom types
    BondTable = [[0.0 for _ in range(MyCrystal.ntypat)] for _ in range(MyCrystal.ntypat)]
    with open(InputFile) as ff:
        for line in ff:
            line=line.strip()
            entry=line.split()
            if (len(entry)==3):
                for ii in range(MyCrystal.ntypat):
                    if MyCrystal.elements[ii]==entry[0]:
                        for jj in range(MyCrystal.ntypat):
                            if MyCrystal.elements[jj]==entry[1]:
                                BondTable[ii][jj]=float(entry[2])*float(entry[2])
                                BondTable[jj][ii]=float(entry[2])*float(entry[2])
    return BondTable

def WriteBonding_C_full(MySnapshotL,step,MyCrystal,BondTable,timestep,natom,numCells,acell,specifics):
    
    print("calculating bonds in snapshot n. "+str(step))
    #Preparing data to be used by the C function
    MSp = np.array(list(MySnapshotL)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))    
    BTp = np.array(BondTable).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    CrystalTypes = []
    for ty in range(MyCrystal.ntypat):
        CrystalTypes+=[ty for _ in range(MyCrystal.types[ty])]
        
    CTp = (ctypes.c_int * natom)(*CrystalTypes)

    if specifics ==None :
        
        Bonds = fullbond_lib.compute_fBonds(MSp,BTp,CTp,natom,MyCrystal.ntypat,acell[0],acell[1],acell[2],numCells)
    
    else :
        SPp = (ctypes.c_double * 6)(*specifics)
        Bonds = fullbond_lib.compute_fBonds_specific(MSp,BTp,CTp,natom,MyCrystal.ntypat,acell[0],acell[1],acell[2],numCells,SPp)
    
    LBonds = Bonds[0]

    #Converting bond profile from C to python list
    BondsList=[[at] for at in range(natom)]
    atom=0
    for i in range(1,LBonds):
        if Bonds[i]!=-1:
            BondsList[Bonds[i]].append(atom)
        else : 
            atom+=1
                
    fullbond_lib.free_memory(Bonds)
#    fullbond_lib.free_memory_int(CTp)
#    fullbond_lib.free_memory_double(BTp)
#    fullbond_lib.free_memory_double(MSp)

    return BondsList        

def main(argv):
    umdpf.headerumd()
    UMDname='output.umd.dat'
    Nsteps = 1
    InputFile = ''
    header = ''
    numCells=None
    maxlength=None
    specifics = None
    try:
        opts, arg = getopt.getopt(argv,"hf:s:l:i:n:p:",["fUMDfile","sSampling_Frequency","lMaxLength","iInputFile","nNumCells","pSpecifics"])
    except getopt.GetoptError:
        print ('Bond_fast_specific.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -i <InputFile> -n <NumCells> -p <Specifics>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Computation of the bonding map for each snapshot of a umd file')
            print ('bond.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -i <InputFile> -n <NumCells>')
            print ('default values: -f output.umd.dat -s 1 -l None')
            print ('The input file contains the bond lengths for the different atom pairs. \n The option -l overwrites all the values of this file.')
            print ('-l : unique bonding length. Use in the case of not having an input file -i. ')
            print ('-p : restricts the computation of the bonds in a specific place of the simulation. Takes a list of floats under the form -p [xmin,ymin,zmin,xmax,ymax,zmax], with the list contains the maximal and minimal coordinates values for the atoms to be considered. Default None (no restriction).')
            print ('-n : states the number of sub-cells the script will work with. The default value is automatically computed for optimal performances.')
            sys.exit()
        elif opt in ("-l","--lmaxlength"):
            if(arg!=""):
                maxlength = float(arg)    
                header = header + "-l=" + str(maxlength)
        elif opt in ("-f", "--fUMDfile"):
            UMDname = str(arg)
            header = header + 'FILE: -f=' + UMDname
        elif opt in ("-s","--sNsteps"):
            Nsteps = int(arg)
            header = header + ' -s=' + arg
            print('Will sample the trajectories every ',Nsteps,' steps')
        elif opt in ("-i", "--iInputFile"):
            InputFile = str(arg)
            header = header + ' -i=' + InputFile
            print ('Bonding cutoffs to be read from file ',InputFile)
        elif opt in ("-n", "--nNumCells"):
            numCells = int(arg)
            header = header + ' -n=' + arg
        elif opt in ("-p","--pSpecifics"):
            specifics = eval(arg)

    if not (os.path.isfile(UMDname)):
        print ('the UMD files ',UMDname,' does not exist')            
        sys.exit()

    start=time.time()

            
#    (MyCrystal,AllSnapshotsLA,acell,TimeStep)=read_xcart_only(UMDname,Nsteps)#reads the cartesian coordinates and put then in AllSnapshotsL

    MyCrystal,AllSnapshotsL,TimeStep,lensnap = umdpf.read_values(UMDname,"xcart","line",Nsteps)
    acell=MyCrystal.acell

    if maxlength==None and len(InputFile)>0 :
        BondTable = read_inputfile(InputFile,MyCrystal)#Converts the input file into a matrix
    elif maxlength!=None :
        BondTable = [[maxlength**2 for _ in range(MyCrystal.ntypat)] for _ in range(MyCrystal.ntypat)]#Creates a matrix filled with -l value


    M=math.sqrt(max([max(Bondlengths) for Bondlengths in BondTable]))#maximal length of any bond

    if numCells == None :
        numCells = int(min(acell)/M)
        print("Number of cells automatically fixed to "+str(numCells))
        
    if acell[0]/numCells<M or acell[1]/numCells<M or acell[2]/numCells<M:
        print('WARNING : one or more dimension(s) of the sub-cells smaller than the greatest bond length.')            

    natom=MyCrystal.natom
    if UMDname[-8:] == ".umd.dat":
        FileAll=UMDname[:-8]+'.bonding.dat'
    else :
        FileAll=UMDname+'.bonding.dat'

    print ('Bonds will be written in <',FileAll,'> file')
    
    umdpf.copyhead(FileAll,UMDname,Nsteps)     
    
    WriteBondingRed=partial(WriteBonding_C_full,MyCrystal=MyCrystal,BondTable=BondTable,timestep=TimeStep,natom=natom,numCells=numCells,acell=acell,specifics = specifics)
        
    
    with concurrent.futures.ProcessPoolExecutor() as executor :
        Lines=list(executor.map(WriteBondingRed,AllSnapshotsL,[step*Nsteps for step in range(len(AllSnapshotsL))]))#Calculates the bond profile for each snapshot


    step=0
    print("Writing...")
    fa=open(FileAll,'a')
    
    if specifics != None :
        line = "Space interval ["+str(specifics[0])+","+str(specifics[3])+"] "+"["+str(specifics[1])+","+str(specifics[4])+"] "+"["+str(specifics[2])+","+str(specifics[5])+"]\n"
        fa.write(line)
    for lines in Lines :                
        header = 'time '+str(step*TimeStep*Nsteps)+' fs\nstep '+str(step)+'\n'
        fa.write(header)
        for line in lines :
            l=''
            for atom in line :
                l+=str(atom)+'\t'
            l+='\n'
            fa.write(l)
                
        fa.write('end\n')
        step+=1
                            
    fa.close()
    
    
    print ('Bonds written in <',FileAll,'> file')

    end = time.time()
    print("runtime:",end-start," s")    

if __name__ == "__main__":
    main(sys.argv[1:])