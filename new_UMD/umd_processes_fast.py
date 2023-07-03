#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:22:44 2023

@author: k
"""
###
##AUTHORS: KEVIN JIGUET
###


import sys, os.path
import crystallography as cr
from functools import partial
import concurrent.futures
from os.path import join
import ctypes
import numpy as np
import time

current_path=os.path.abspath(__file__)#For all this to work, the file < c_UMDprocess.so > must be in the same directory than this script
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u

read_lib = ctypes.cdll.LoadLibrary(join(path_new, 'c_UMDprocess.so'))
read_lib.defineBonds.argtypes = [ctypes.POINTER(ctypes.c_char),ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
read_lib.defineBonds.restype = ctypes.POINTER(ctypes.c_int)
read_lib.read_umd_values.argtypes = [ctypes.POINTER(ctypes.c_char),ctypes.c_int,ctypes.c_int]
read_lib.read_umd_values.restype = ctypes.POINTER(ctypes.c_double)
read_lib.free_memory.argtypes = [ctypes.POINTER(ctypes.c_double)]
read_lib.free_memory.restype = None


def read_snapshot_values_C(octettop,octetbot,step,natom,File,X,mode):#Extracts the coordinates from a snapshot 
                                                                #whose data is contained between the octets "octettop" and "octetbot" in the umd file (thus the prep_read_coord function)
                                                                #The parameter X is the first index of the relevant 3 coordinate on each line in the umd file
    ff=open(File,"r")
    ff.seek(octettop,0)
    snapshot=ff.read(octetbot-octettop)
    ff.close()
    snapPointer = ctypes.c_char_p(snapshot.encode('utf-8'))
    
    CList=read_lib.read_umd_values(snapPointer,natom,X)
    SnapshotValues=[]
    if X==-1:
        k=12
    else :
        k=3
    #Converts the C data into a python list
    if mode=="lists":#Returns a list of 3-elements list (x,y,z) for each atom
        for i in range(natom):
            List=[]
            for j in range(k):
                List.append(round(CList[k*i+j],5))
            SnapshotValues.append(List)
    elif mode=="line":#Returns a unique list of values for all the atoms, atom after atom
        for i in range(k*natom):
            SnapshotValues.append(round(CList[i],5))
    elif mode=="chunks":#Returns a unique list of values for all the atoms, axis after axis
        for j in range(k):
            for i in range(natom):
                SnapshotValues.append(round(CList[k*i+j],5))
    
    read_lib.free_memory(CList)
    
    return SnapshotValues

def read_snapshotbonds_C(octettop,octetbot,CentMin,CentMax,OutMin,OutMax,File):#Converts the bonding information of a given snapshot into a list of atoms 
                                                                               #the data is contained between the octets "octettop" and "octetbot" in the bonding file
    
    ff=open(File,"r")
    ff.seek(octettop,0)
    snapshot=ff.read(octetbot-octettop)
    ff.close()

    snapPointer = ctypes.c_char_p(snapshot.encode('utf-8'))
    
    BList = read_lib.defineBonds(snapPointer,len(snapshot),CentMin,CentMax,OutMin,OutMax)
    
    SnapshotBondingTable = []
    SnapshotBondIndexes = []
   
    #Converts the C data into two python lists
    for i in range(2,BList[0]+1):
        SnapshotBondingTable.append(BList[i])#This one contains the atoms
    for i in range(BList[0]+2,BList[0]+BList[1]+2):
        SnapshotBondIndexes.append(BList[i])#This one contains the information about which is bound to which in the other list
        
    return SnapshotBondingTable,SnapshotBondIndexes

def Octets_from_File(File,key,nlines):#Returns the indexes of octets pertaining to each snapshot in the file.
    ff=open(File,"r")                 #"key" is the pattern that separates the successive snapshots in the file
    OctIndexes = []                   #"nlines" is the number of lines pertaining to atoms in each snapshot

    print("Preparing the file ",File," to be read...")
    while True :
        line = ff.readline()
        if not line :
            break
        line=line.strip().split()
        if len(line)>0 and line[0]==key:
            OctIndexes.append(ff.tell())
            for _ in range(nlines):
                ff.readline()
    OctIndexes.append(ff.tell())
    ff.close()
    return OctIndexes

def Crystallization(File):#builds the Crystal
    MyCrystal = cr.Lattice()
    ff=open(File,"r")
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
            elif entry[0] == 'acell':
                MyCrystal.acell=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprim_a':
                MyCrystal.rprim[0]=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprim_b':
                MyCrystal.rprim[1]=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprim_c':
                MyCrystal.rprim[2]=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprimd_a':
                MyCrystal.rprimd[0]=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprimd_b':
                MyCrystal.rprimd[1]=[float(entry[1]),float(entry[2]),float(entry[3])]
            elif entry[0] == 'rprimd_c':
                MyCrystal.rprimd[2]=[float(entry[1]),float(entry[2]),float(entry[3])]
            if entry[0] == 'atoms:':
                break
    ff.close()
    return MyCrystal, TimeStep


def read_values(UMDfile,key,mode="line",Nsteps=1,cutoff="all"):
    
    t=time.time()
    MyCrystal,TimeStep = Crystallization(UMDfile)
    OctIndexesRaw = Octets_from_File(UMDfile, 'atoms:', MyCrystal.natom)
    
    if cutoff=="all":
        cutoff = len(OctIndexesRaw)-1
        
    OctIndexes=OctIndexesRaw[:cutoff+1]
    OctTop = [OctIndexes[i*Nsteps] for i in range(int((len(OctIndexes)-1)/Nsteps))]#Indicates the beginning of each relevant snapshot
    OctBot = [OctIndexes[i*Nsteps+1] for i in range(int((len(OctIndexes)-1)/Nsteps))]#Indicates it's end


    if key == "xred":
        X=0
        print('Extracting the reduced coordinates from the file ',UMDfile)
 
    elif key == "xcart":
        X=3
        print('Extracting the cartesian coordinates from the file ',UMDfile)

    elif key == "absxcart":
        X=6
        print('Extracting the absolute coordinates from the file ',UMDfile)

    elif key == "velocity":
        X=9
        print('Extracting the velocities from the file ',UMDfile)
    elif key == "everything":
        X=-1

    else :
        print("Parameter < ",key," > not recognized")
        sys.exit()

    if mode!="line" and mode !="lists" and mode!="chunks":
        print("Parameter < ",mode," > not recognized")
        sys.exit()
    
        
    
    readvalues = partial(read_snapshot_values_C,natom=MyCrystal.natom,File=UMDfile,X=X,mode=mode)
    with concurrent.futures.ProcessPoolExecutor() as executor :
        SnapshotsValuesList = list(executor.map(readvalues,OctTop,OctBot,[step for step in range(len(OctIndexes)-1)]))        

    print("Values extracted in ",(time.time()-t)," s")
    return MyCrystal, SnapshotsValuesList, TimeStep, len(OctIndexes)-1


def data_type(SnapshotsValuesList,natom,mode="line",datatype="list"):

    if mode=="atoms":
        SnapshotsValuesX = [[[x[3*i],x[3*i+1],x[3*i+2]]  for i in range(natom)] for x in SnapshotsValuesList]
    else :
        SnapshotsValuesX = SnapshotsValuesList

    if datatype=="array":
            SnapshotsValues = [np.array(x) for x in SnapshotsValuesX]
    else :
            SnapshotsValues = SnapshotsValuesX
    
    return SnapshotsValues

def read_bonds(BondFile,centEl,adjEl,Nsteps=1):
    
    ff=open(BondFile,"r")
    natom = int(ff.readline().strip().split()[1])
    
    while True :
        
        line = ff.readline().strip().split()
        if not line : break
        if len(line)>0:
            if line[0] == "types":
                types = [0]+[int(x) for x in line[1:]]
            elif line[0] == "elements":    
                elements = [""]+line[1:]
            elif line[0] == "atoms":
                break
            
    ff.close()
    print(elements[1:])
    
    try:
        indCent = elements.index(centEl)
    except :
        print("ERROR : element ",centEl," not present in simulation")
        sys.exit()        
    try:
        indAdj = elements.index(adjEl)
    except :
        print("ERROR : element ",adjEl," not present in simulation")
        sys.exit()

    
    CentMin = sum(types[:indCent])
    CentMax = CentMin + types[indCent]-1
    AdjMin = sum(types[:indAdj])
    AdjMax = AdjMin + types[indAdj]-1
        
    MyCrystal,TimeStep = Crystallization(BondFile)
    OctIndexes = Octets_from_File(BondFile,"step",natom)

    OctTop = [OctIndexes[i*Nsteps] for i in range(int((len(OctIndexes)-1)/Nsteps))]#Indicates the beginning of each relevant snapshot
    OctBot = [OctIndexes[i*Nsteps+1] for i in range(int((len(OctIndexes)-1)/Nsteps))]#Indicates it's end

        
    print("Extracting bonds from file ",BondFile)
    
    bondsRed = partial(read_snapshotbonds_C,CentMin=CentMin, CentMax=CentMax, OutMin=AdjMin, OutMax=AdjMax, File=BondFile)
    
    with concurrent.futures.ProcessPoolExecutor() as executor :
        Data = list(executor.map(bondsRed,OctTop,OctBot))
    Bonds = [D[0] for D in Data]
    BondsIndexes = [D[1] for D in Data]
    
    return CentMin, CentMax, AdjMin, AdjMax, MyCrystal, Bonds, BondsIndexes, TimeStep


def read_stresses_4visc(UMDfile):
    ff=open(UMDfile,"r")
    AllSnapshots=[]
    
    line = ff.readline().strip().split()
    natom=int(line[1])
    TimeStep=0
    
    while True : 
        SnapshotCrystal = cr.Lattice()
        flagheader=1
        while flagheader :
            l = ff.readline()
            if not l :
                l=ff.readline()
                if not l :#If we encounter two empty lines in a row, we're at the end of the file
                    print('len of allsnapshots is',len(AllSnapshots))
                    return AllSnapshots,TimeStep
            line = l.strip().split()
            if len(line)>0 :
                if line[0]=='timestep':
                    TimeStep = float(line[1])
                if line[0]=='atoms:':
                    flagheader=0
                elif line[0]=='rprimd_a':
                    SnapshotCrystal.rprimd[0]=[float(line[ii]) for ii in range(1,4)]
                elif line[0]=='rprimd_b':
                    SnapshotCrystal.rprimd[1]=[float(line[ii]) for ii in range(1,4)]
                elif line[0]=='rprimd_c':
                    SnapshotCrystal.rprimd[2]=[float(line[ii]) for ii in range(1,4)]
                elif line[0]=='StressTensor':
                    SnapshotCrystal.stress=[float(line[i]) for i in range(1,7)]
                elif line[0]=='Temperature':
                    SnapshotCrystal.temperature=float(line[1])
 
        SnapshotCrystal.cellvolume = SnapshotCrystal.makevolume()
        AllSnapshots.append(SnapshotCrystal)
        
        for _ in range(natom):
            next(ff)

def headerumd():
    print ('\n * * * * * ')
    print ('    UMD package for analyzing Molecular Dynamics simulations.')
    print ('distributed under GNU-GNU General Public License 3')
    print ('please cite as: ')
    print ('    Razvan Caracas, Jean-Alexis Hernandez, Anais Kobsch, Zhi Li, Natalia Solomatova, and Francois Soubiran' )
    print ('    UMD: An open source package for analyzing Molecular Dynamics simulations' )
    print ('    Journal of Visualized Experiments, in press/media prep (2020)' )
    print (' ')
    
def print_header(FileName,MyCrystal):
    newfile = FileName + '.umd.dat'
    nf = open(newfile,'a')
    string = 'natom ' + str(MyCrystal.natom) + '\n'
    nf.write(string)
    string = 'ntypat ' + str(MyCrystal.ntypat) + '\n'
    nf.write(string)
    string = 'no.of.electrons ' + str(MyCrystal.noelectrons) + '\n'
    nf.write(string)
    string = 'types '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.types[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'elements '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.elements[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'masses '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.masses[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'Zelectrons '
    for ii in range(MyCrystal.ntypat):
        string = string + str(MyCrystal.zelec[ii]) + ' '
    string = string + '\n'
    nf.write(string)
    string = 'typat '
    for ii in range(MyCrystal.natom):
        string = string + str(MyCrystal.typat[ii]) + ' '
    string = string + '\n\n'
    nf.write(string)
    string = 'lambda_ThermoInt' + ' ' +  str(MyCrystal.lambda_ThermoInt)
    string = string + '\n\n'
    nf.write(string)
    nf.close()

def print_snapshots(FileName,MyCrystal,TimeStep,CurrentTime,diffcoords):
    newfile = FileName + '.umd.dat'
    nf = open(newfile,'a')
    string = 'timestep ' + str(TimeStep) + ' fs\n' + 'time ' + str(CurrentTime) + ' fs\n'
    nf.write(string)
    string = 'InternalEnergy ' + str(round(MyCrystal.internalenergy,6)) + ' eV\n'
    nf.write(string)
    string = 'ElectronicEntropy ' + str(round(MyCrystal.electronicentropy,6)) + ' eV\n'
    nf.write(string)
    string = 'KineticEnergyIons ' + str(round(MyCrystal.kineticenergy,6)) + ' eV\n'
    nf.write(string)
    string = 'EnergyWithDrift ' + str(round(MyCrystal.energywithdrift,6)) + ' eV\n'
    nf.write(string)
    string = 'Enthalpy ' + str(round(MyCrystal.enthalpy,6)) + ' eV\n'
    nf.write(string)
    string = 'Magnetization ' + str(round(MyCrystal.magnetization,6)) + ' Magneton-Bohr\n'
    nf.write(string)
    string = 'Temperature ' + str(round(MyCrystal.temperature,1)) + ' K\n'
    nf.write(string)
    string = 'Pressure ' + str(round(MyCrystal.pressure,4)) + ' GPa\n'
    nf.write(string)
    string = 'Density ' + str(round(MyCrystal.density,3)) + ' g.cm-3\n'
    nf.write(string)
    string = 'StressTensor '
    for ii in range(6):
        string = string + str(round(MyCrystal.stress[ii],4)) + ' '
    string = string + ' GPa\n'
    nf.write(string)
#    string = 'acell ' + str(round(MyCrystal.acell[0],3)) + ' ' + str(round(MyCrystal.acell[1],3) + ' ' + str(round(MyCrystal.acell[2],3)) + ' A\n'
    string = 'acell ' + str(MyCrystal.acell[0]) + ' ' + str(MyCrystal.acell[1]) + ' ' + str(MyCrystal.acell[2]) + ' A\n'
    nf.write(string)
    string = 'rprim_a ' + str(MyCrystal.rprimd[0][0]/MyCrystal.acell[0]) + '  ' +str(MyCrystal.rprimd[0][1]/MyCrystal.acell[0]) + '  ' +str(MyCrystal.rprimd[0][2]/MyCrystal.acell[0]) + '\n'
    nf.write(string)
    string = 'rprim_b ' + str(MyCrystal.rprimd[1][0]/MyCrystal.acell[1]) + '  ' +str(MyCrystal.rprimd[1][1]/MyCrystal.acell[1]) + '  ' +str(MyCrystal.rprimd[1][2]/MyCrystal.acell[1]) + '\n'
    nf.write(string)
    string = 'rprim_c ' + str(MyCrystal.rprimd[2][0]/MyCrystal.acell[2]) + '  ' +str(MyCrystal.rprimd[2][1]/MyCrystal.acell[2]) + '  ' +str(MyCrystal.rprimd[2][2]/MyCrystal.acell[2]) + '\n'
    nf.write(string)
    string = 'rprimd_a ' + str(MyCrystal.rprimd[0][0]) + '  ' +str(MyCrystal.rprimd[0][1]) + '  ' +str(MyCrystal.rprimd[0][2]) + ' A\n'
    nf.write(string)
    string = 'rprimd_b ' + str(MyCrystal.rprimd[1][0]) + '  ' +str(MyCrystal.rprimd[1][1]) + '  ' +str(MyCrystal.rprimd[1][2]) + ' A\n'
    nf.write(string)
    string = 'rprimd_c ' + str(MyCrystal.rprimd[2][0]) + '  ' +str(MyCrystal.rprimd[2][1]) + '  ' +str(MyCrystal.rprimd[2][2]) + ' A\n'
    nf.write(string)
    string = 'atoms: reduced*3 cartesian*3(A) abs.diff.*3(A) velocity*3(A/fs) force*3(eV/A) charge(no.elec) magnetization(magneton-Bohr) \n'
    nf.write(string)
    for iatom in range(MyCrystal.natom):
        string = str(round(MyCrystal.atoms[iatom].xred[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xred[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xred[2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].xcart[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xcart[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].xcart[2],5)) + ' '
        string = string + str(round(diffcoords[iatom][0],5)) + ' ' + str(round(diffcoords[iatom][1],5)) + ' ' + str(round(diffcoords[iatom][2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].vels[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].vels[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].vels[2],5)) + ' '
        string = string + str(round(MyCrystal.atoms[iatom].forces[0],5)) + ' ' + str(round(MyCrystal.atoms[iatom].forces[1],5)) + ' ' + str(round(MyCrystal.atoms[iatom].forces[2],5)) + ' '
        string = string + str(MyCrystal.atoms[iatom].charge) + ' ' + str(MyCrystal.atoms[iatom].magnet) + '\n'
        nf.write(string)
    string='\n'
    nf.write(string)
    nf.close()
    return(CurrentTime,TimeStep)

def read_bigheader_umd(umdfile, short=0):
    MyCrystal = cr.Lattice()
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
                    break
    return(MyCrystal,TimeStep)

def copyhead(Target,File):
    
    ff = open(File,"r")
    fa = open(Target,"w")
    
    while True : 
        line = ff.readline()
        if not line :
            ff.close()
            fa.close()
            break
        l = line.strip().split()
        
        if len(l)>0 :
            if l[0] == "atoms:":
                ff.close()
                fa.close()
                break                
        
        fa.write(line)
