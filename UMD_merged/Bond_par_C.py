#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:55:47 2023

@author: Kevin Jiguet-Covex
"""

import sys,getopt,os.path,itertools
import crystallography as cr
import umd_process as umdp
import math
import time
from functools import partial
import concurrent.futures
import numpy as np
import ctypes
from os.path import join

current_path=os.path.abspath(__file__)#For all this to work, the file c_gofr.so must be in the same directory than gofr_umd
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u
#c_msdbtest = 

fullbond_lib = ctypes.cdll.LoadLibrary(join(path_new, 'c_bonds_full.so'))
fullbond_lib.compute_fBonds.argtypes = [ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int]
fullbond_lib.compute_fBonds.restype = ctypes.POINTER(ctypes.c_int)

def read_inputfile(InputFile,MyCrystal):
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
    #print(BondTable)
    return BondTable

def WriteBonding_C_full(MySnapshotL,step,MyCrystal,BondTable,timestep,natom,numCells,acell):
    
    print(step)
#    print(MySnapshotL[:10])
    MSp = np.array(MySnapshotL).ctypes.data_as(ctypes.POINTER(ctypes.c_double))    
    BTp = np.array(BondTable).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    CrystalTypes = []
    for ty in range(MyCrystal.ntypat):
        CrystalTypes+=[ty for _ in range(MyCrystal.types[ty])]
        
    CTp = (ctypes.c_int * natom)(*CrystalTypes)

    Bonds = fullbond_lib.compute_fBonds(MSp,BTp,CTp,natom,MyCrystal.ntypat,acell[0],acell[1],acell[2],numCells)
    M = Bonds[0]

    index = 1
    
    BondsList=[[at] for at in range(natom)]
    atom=0
    
    while index<1+M*27*natom : 
        if Bonds[index]!=-1:
            BondsList[Bonds[index]].append(atom)
            index+=1
        else : 
            atom+=1
            index = M*27*atom+1
    
    return BondsList


def read_xcart_only(umdfile):
    niter = 0
    MyCrystal = cr.Lattice()
    AllSnapshotsList = []
    with open(umdfile,'r') as ff:
        while True:
            line = ff.readline()
            if not line: break
            #print(line,len(line))
            if len(line) > 1:
                line=line.strip()
                entry=line.split()
                if entry[0] == 'natom':
                    MyCrystal.natom = int(entry[1])
                    MyCrystal.typat = [0 for _ in range(MyCrystal.natom)]
                elif entry[0] == 'ntypat':
                    MyCrystal.ntypat = int(entry[1])
                    MyCrystal.types = [0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.elements = ['X' for _ in range(MyCrystal.ntypat)]
                    MyCrystal.masses = [1.0 for _ in range(MyCrystal.ntypat)]
                    MyCrystal.zelec = [1.0 for _ in range(MyCrystal.ntypat)]
                elif entry[0] == 'types':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.types[ii]=int(entry[ii+1])
                elif entry[0] == 'elements':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.elements[ii]=entry[ii+1]
                elif entry[0] == 'masses':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.masses[ii]=float(entry[ii+1])
                elif entry[0] == 'Zelectrons':
                    for ii in range(MyCrystal.ntypat):
                        MyCrystal.zelec[ii]=float(entry[ii+1])
                elif entry[0] == 'typat':
                    for ii in range(MyCrystal.natom):
                        MyCrystal.typat[ii] = int(entry[ii+1])
                elif entry[0] == 'timestep':
                    TimeStep = float(entry[1])
                if entry[0] == 'acell':
                    acell = [0,0,0]
                    for ii in range(3):
                        acell[ii] = float(entry[ii+1])
                

                if entry[0] == 'atoms:':
                    print('reading file : current iteration no.',niter)
                    if niter==5000:
                        break
                    MySnapshot = np.zeros((MyCrystal.natom,3),dtype=float)
                    for iatom in range(MyCrystal.natom):
                        line = ff.readline()
                        line=line.strip()
                        entry=line.split()
                        coord=np.array([0.0,0.0,0.0])
                        for jj in range(3):
                            coord[jj]=float(entry[jj+3])
                        MySnapshot[iatom]=coord
                    AllSnapshotsList.append(MySnapshot)
                    niter += 1
    #!!! remove first element from AllSnapshots
    print('len of allsnapshots is',len(AllSnapshotsList))
#    print([AllSnapshotsList[i][:10] for i in range(len(AllSnapshotsList))])
    return(MyCrystal,AllSnapshotsList,acell,TimeStep)
        

def main(argv):
    umdp.headerumd()
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
            print ('computation of the bonding map')
            print ('bond.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -c <Cations> -a <Anions> -m <MinLife>  -i <InputFile> -r <Rings>')
            print ('  default values: -f output.umd.dat -s 1 -s 1 -l None -n 5')
            print (' the input file contains the bond lengths for the different atom pairs. \n the values overwrite the option -l')
            print (' rings = 1 default, polymerization, all anions and cations bond to each other; rings = 0 only individual cation-anion groups')
            sys.exit()
        elif opt in ("-l","--lmaxlength"):
            if(arg!=""):
                maxlength = float(arg)    
                header = header + "-l=" + str(maxlength)
        elif opt in ("-f", "--fUMDfile"):
            UMDname = str(arg)
            header = header + 'FILE: -f=' + UMDname
#            maxlength=None
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
        elif opt == "-v":
            v=str(arg)

    if not (os.path.isfile(UMDname)):
        print ('the UMD files ',UMDname,' does not exist')            
        sys.exit()

    start=time.time()

        

    (MyCrystal,AllSnapshotsL,acell,TimeStep)=read_xcart_only(UMDname)
    if maxlength==None and len(InputFile)>0 :
        BondTable = read_inputfile(InputFile,MyCrystal)
    elif maxlength!=None :
        BondTable = [[maxlength**2 for _ in range(MyCrystal.ntypat)] for _ in range(MyCrystal.ntypat)]


    M=math.sqrt(max([max(Bondlengths) for Bondlengths in BondTable]))

    if acell[0]/numCells<M or acell[1]/numCells<M or acell[2]/numCells<M:
        print('WARNING : dimension of the cell smaller than the greatest bondlength.')
        print(acell[0]/numCells)
    
        

    natom=MyCrystal.natom
    FileAll=UMDname[:-8]+'.bondingfile.dat'
    print ('Bondings will be written in <',FileAll,'> file')
    fa=open(FileAll,'w')
    ff=open(UMDname,'r')
    
    for i in range(20):    
        line=ff.readline()
        fa.write(line)
        
    ff.close()
    Lines=[]
    WriteBondingRed=partial(WriteBonding_C_full,MyCrystal=MyCrystal,BondTable=BondTable,timestep=TimeStep,natom=natom,numCells=numCells,acell=acell)
        
    with concurrent.futures.ProcessPoolExecutor() as executor :
        Lines=list(executor.map(WriteBondingRed,AllSnapshotsL,[step for step in range(0,len(AllSnapshotsL),Nsteps)]))

    step=0
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
    
    end=time.time()

    print("runtime:",end-start)    

if __name__ == "__main__":
    main(sys.argv[1:])
