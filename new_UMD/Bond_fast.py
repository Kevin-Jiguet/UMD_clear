#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:55:47 2023

@author: Kevin Jiguet-Covex
"""

import sys,getopt,os.path,itertools
import math
import time
from functools import partial
import concurrent.futures
import numpy as np
import ctypes
from os.path import join
import umd_processes_fast as umdpf

current_path=os.path.abspath(__file__)#For all this to work, the file c_bonds_fullD.so must be in the same directory than gofr_umd
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u
#c_msdbtest = 

fullbond_lib = ctypes.cdll.LoadLibrary(join(path_new, 'c_bonds_fullD.so'))
fullbond_lib.compute_fBonds.argtypes = [ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int]
fullbond_lib.compute_fBonds.restype = ctypes.POINTER(ctypes.c_int)

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

def WriteBonding_C_full(MySnapshotL,step,MyCrystal,BondTable,timestep,natom,numCells,acell):
    
    print("calculating bonds in snapshot n. "+str(step))
    #Preparing data to be used by the C function
    MSp = np.array(list(MySnapshotL)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))    
    BTp = np.array(BondTable).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    CrystalTypes = []
    for ty in range(MyCrystal.ntypat):
        CrystalTypes+=[ty for _ in range(MyCrystal.types[ty])]
        
    CTp = (ctypes.c_int * natom)(*CrystalTypes)

    Bonds = fullbond_lib.compute_fBonds(MSp,BTp,CTp,natom,MyCrystal.ntypat,acell[0],acell[1],acell[2],numCells)
    LBonds = Bonds[0]

    #Converting bond profile from C to python list
    BondsList=[[at] for at in range(natom)]
    atom=0
    for i in range(1,LBonds):
        if Bonds[i]!=-1:
            BondsList[Bonds[i]].append(atom)
        else : 
            atom+=1
                
    
    return BondsList        

def main(argv):
    umdpf.headerumd()
    UMDname='output.umd.dat'
    Nsteps = 1
    InputFile = ''
    header = ''
    numCells=5
    maxlength=None
    try:
        opts, arg = getopt.getopt(argv,"hf:s:l:i:n:",["fUMDfile","sSampling_Frequency","lMaxLength","iInputFile","nNumCells"])
    except getopt.GetoptError:
        print ('Bond_par.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -i <InputFile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Computation of the bonding map for each snapshot of a umd file')
            print ('bond.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -i <InputFile> -n <NumCells>')
            print ('  default values: -f output.umd.dat -s 1 -l None -n 5')
            print ('The input file contains the bond lengths for the different atom pairs. \n The option -l overwrites all the values of this file.')
            print ('The option -n states the number of sub-cells the script will work with. \n For an optimal functionment, -n should be as high as possible \n as long as the dimension of the sub-cells is higher than the biggest bond length.')
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
            print('Will sample the MD trajectory every ',Nsteps,' steps')
        elif opt in ("-i", "--iInputFile"):
            InputFile = str(arg)
            header = header + ' -i=' + InputFile
            print ('Bonding cutoffs to be read from file ',InputFile)
        elif opt in ("-n", "--nNumCells"):
            numCells = int(arg)
            header = header + ' -n=' + arg

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

    if acell[0]/numCells<M or acell[1]/numCells<M or acell[2]/numCells<M:
        print('WARNING : one or more dimension(s) of the sub-cells smaller than the greatest bond length.')
            

    natom=MyCrystal.natom
    FileAll=UMDname[:-8]+'.bondingB.dat'
    print ('Bonds will be written in <',FileAll,'> file')
    fa=open(FileAll,'w')
    ff=open(UMDname,'r')
    
    for i in range(20):    
        line=ff.readline()
        fa.write(line)
        
    ff.close()
    Lines=[]
    WriteBondingRed=partial(WriteBonding_C_full,MyCrystal=MyCrystal,BondTable=BondTable,timestep=TimeStep,natom=natom,numCells=numCells,acell=acell)
        
    with concurrent.futures.ProcessPoolExecutor() as executor :
        Lines=list(executor.map(WriteBondingRed,AllSnapshotsL,[step*Nsteps for step in range(len(AllSnapshotsL))]))#Calculates the bond profile for each snapshot


    step=0
    print("Writing...")
    for lines in Lines :                
        header = 'time '+str(step*TimeStep)+' fs\nstep '+str(step)+'\n'
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

    
    end=time.time()

    print("runtime:",end-start," s")    

if __name__ == "__main__":
    main(sys.argv[1:])
