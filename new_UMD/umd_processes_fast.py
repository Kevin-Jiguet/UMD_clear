#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:22:44 2023

@author: k
"""
###
##AUTHORS: KEVIN JIGUET
###


import sys,getopt,os.path,itertools
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
read_lib.read_umd_values.argtypes = [ctypes.POINTER(ctypes.c_char),ctypes.c_int,ctypes.c_int,ctypes.c_int]
read_lib.read_umd_values.restype = ctypes.POINTER(ctypes.c_double)

def read_snapshot_values_C(octettop,octetbot,step,natom,File,X):#Extracts the coordinates from a snapshot 
                                                                #the data is contained between the octets "octettop" and "octetbot" in the umd file (thus the prep_read_coord function)
                                                                #The parameter X is the first index of the relevant coordinate on a line in the umd file
    ff=open(File,"r")
    ff.seek(octettop,0)
    snapshot=ff.read(octetbot-octettop)
    ff.close()
    snapPointer = ctypes.c_char_p(snapshot.encode('utf-8'))
    
    CList = read_lib.read_umd_values(snapPointer,natom,len(snapshot),X)

    SnapshotValues=[]
    #Converts the C data into a python list


    for i in range(3*natom):
        SnapshotValues.append(round(CList[i],5))
    
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
    ff=open(File,"r")                 #"key" is the pattern that separates snapshots in the file
    OctIndexes = []                   #"nlines"is the number of lines (atoms) in each snapshot

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
    for _ in range(20):
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
            if entry[0] == 'acell':
                MyCrystal.acell=[float(entry[1]),float(entry[2]),float(entry[3])]
    ff.close()
    return MyCrystal, TimeStep



def read_values(UMDfile,key,Nsteps=1):
    MyCrystal,TimeStep = Crystallization(UMDfile)
    OctIndexes = Octets_from_File(UMDfile, 'atoms:', MyCrystal.natom)
    
    OctTop = [OctIndexes[i*Nsteps] for i in range(int((len(OctIndexes)-1)/Nsteps))]#Indicates the beginning of each relevant snapshot
    OctBot = [OctIndexes[i*Nsteps+1] for i in range(int((len(OctIndexes)-1)/Nsteps))]#Indicates it's end

    if key == "xcart":
        X=3
        print('Extracting the cartesian coordinates from the file ',UMDfile)

    elif key == "absxcart":
        X=6
        print('Extracting the absolute coordinates from the file ',UMDfile)

    elif key == "velocity":
        X=9
        print('Extracting the velocities from the file ',UMDfile)

    else :
        print("Parameter < ",key," > not recognized")
        sys.exit()

    readvalues = partial(read_snapshot_values_C,natom=MyCrystal.natom,File=UMDfile,X=X)
    with concurrent.futures.ProcessPoolExecutor() as executor :
        SnapshotsValuesList = list(executor.map(readvalues,OctTop,OctBot,[step for step in range(len(OctIndexes)-1)]))        

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

def read_bonds(BondFile,centEl,adjEl):
    
    ff=open(BondFile,"r")
    natom = int(ff.readline().strip().split()[1])
    ff.readline()
    typesStr = ff.readline().strip().split()[1:]
    types = [0]+[int(x) for x in typesStr]
    elements = [""]+ff.readline().strip().split()[1:]
    ff.close()
    print(elements)
    
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

        
    print("Extracting bonds from file ",BondFile)
    
    bondsRed = partial(read_snapshotbonds_C,CentMin=CentMin, CentMax=CentMax, OutMin=AdjMin, OutMax=AdjMax, File=BondFile)
    
    with concurrent.futures.ProcessPoolExecutor() as executor :
        Data = list(executor.map(bondsRed,OctIndexes[:-1],OctIndexes[1:]))
    Bonds = [D[0] for D in Data]
    BondsIndexes = [D[1] for D in Data]
    
    return CentMin, CentMax, AdjMin, AdjMax, MyCrystal, Bonds, BondsIndexes, TimeStep
