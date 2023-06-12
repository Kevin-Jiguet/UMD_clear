#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:44:43 2023
"""
#!/usr/bin/env python3
###
##AUTHORS: KEVIN JIGUET
###


import sys,getopt,os.path,itertools
import crystallography as cr
import umd_process as umdp
import math
import time
from functools import partial
import concurrent.futures
from os.path import join
import ctypes


current_path=os.path.abspath(__file__)#For all this to work, the files < c_clusters_B.so > and < c_define_bonds.so > must be in the same directory than gofr_umd
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u

clust_lib = ctypes.cdll.LoadLibrary(join(path_new, 'c_clusters_B.so'))
clust_lib.fullclustering.argtypes = [ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
clust_lib.fullclustering.restype = ctypes.POINTER(ctypes.c_int)

bond_lib = ctypes.cdll.LoadLibrary(join(path_new, 'c_define_bonds.so'))
bond_lib.defineBonds.argtypes = [ctypes.POINTER(ctypes.c_char),ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
bond_lib.defineBonds.restype = ctypes.POINTER(ctypes.c_int)

def analysis_subtab(Clusters,Step):#Creates a dictionnary whose keys are the individual clusters, and the values are a list of list of their *consecutive* steps of existence. 
        population_parallel={}
        for Snapshot in Clusters :
            for cluster in Snapshot :
                if str(cluster) in population_parallel.keys():
                    if population_parallel[str(cluster)][-1][-1]==Step-1:
                        population_parallel[str(cluster)][-1].append(Step)
                    else :
                        population_parallel[str(cluster)].append([Step])
                else :
                    population_parallel[str(cluster)]=[[Step]]
            Step+=1
        return population_parallel
    

def clustering(SnapshotBonds,SnapshotBondIndexes,step,CentMin,CentMax,OutMin,OutMax,r):
    
    print("calculating species from snapshot n. "+str(step))
    M=SnapshotBondIndexes[-1]#maximum number of bound atoms for any atoms
    #Preparing data for the C script        
    SBp = (ctypes.c_int * len(SnapshotBonds))(*SnapshotBonds)
    BIp = (ctypes.c_int * (len(SnapshotBondIndexes)-1))(*SnapshotBondIndexes[:-1])

    Np = clust_lib.fullclustering(SBp,BIp,CentMin,CentMax,OutMin,OutMax,M,r)#Actually computes the clusters
    
    Clusters = [[]]
    length = Np[0]
    
    #Converting C data into a python list of clusters
    for i in range(1,length):
        atom=Np[i]
        if atom ==-1 :
            Clusters[-1].sort()
            Clusters.append([])
        else :   
            Clusters[-1].append(atom)

    if Clusters ==[[]]:
        return []
        
    return Clusters
    
def prep_read(BondFile):#Returns the indexes of octets pertaining to each snapshot in the bonding file ; builds the Crystal
    MyCrystal = cr.Lattice()
    ff=open(BondFile,"r")
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
    
    OctIndexes = []
    while True :
        line = ff.readline()
        if not line :
            break
        line=line.strip().split()
        if line[0]=="step":
            print("preparing snapshot n. "+line[1]+" to be read")
            OctIndexes.append(ff.tell())
            for _ in range(MyCrystal.natom):
                ff.readline()
    OctIndexes.append(ff.tell())
    ff.close()
    return MyCrystal, OctIndexes, TimeStep

def read_snapshotbonds_C(octettop,octetbot,step,CentMin,CentMax,OutMin,OutMax,File):#Converts the bonding information of a given snapshot into a list of atoms 
    print("extracting bonds from snapshot n. "+str(step))                           #the data is contained between the octets "octettop" and "octetbot" in the bonding file
    
    ff=open(File,"r")
    ff.seek(octettop,0)
    snapshot=ff.read(octetbot-octettop)
    ff.close()

    snapPointer = ctypes.c_char_p(snapshot.encode('utf-8'))
    
    BList = bond_lib.defineBonds(snapPointer,len(snapshot),CentMin,CentMax,OutMin,OutMax)
    
    SnapshotBondingTable = []
    SnapshotBondIndexes = []
   
    #Converts the C data into two python lists
    for i in range(2,BList[0]+1):
        SnapshotBondingTable.append(BList[i])#This one contains the atoms
    for i in range(BList[0]+2,BList[0]+BList[1]+2):
        SnapshotBondIndexes.append(BList[i])#This one contains the information about which is bound to which in the other list
        
    return SnapshotBondingTable,SnapshotBondIndexes

def main(argv):
    umdp.headerumd()
    BondFile='bonding.umd.dat'
    Nsteps = 1
    ClusterAtoms = []
    Cations = []
    Anions = []
    minlife = 5
    rings = 1
    header = ''
    nChunks=8
    start=time.time()
    try:
        opts, arg = getopt.getopt(argv,"hf:s:c:a:m:r:n:",["fBondFile","sSampling_Frequency", "cCations","aAnions","mMinlife","rRings","nnChunks"])
    except getopt.GetoptError:
        print ('speciation.py -f <bond_filename> -s <Sampling_Frequency> -c <Cations> -a <Anions> -m <MinLife> -r <Rings> -n <nChunks>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('speciation_fast.py program to identify speciation')
            print ('speciation_fast.py -f <Bonding_filename> -s <Sampling_Frequency> -l <MaxLength> -c <Cations> -a <Anions> -m <MinLife> -r <Rings> -n <nChunks>')
            print ('  default values: -f bonding.umd.dat -s 1 -m 0 -r 1 -n 4')
            print (' the bond file contains the bonds relations for each snapshot. Computed with Bond_par_C.py.')
            print (' -r : rings = 1 default, all anions bind to cations ; rings = 0, polymerization, all anions bind to cations AND anions ; rings = x>0, all anions bind to cations then to other anions to form a xth-coordination polyhedra')
            print (' -m : minimal duration of existence for a chemical species to be taken into account (fs) ; default 5')
            print (' -n : number of chunks for the praallelization of the clusters data ; default 8')
            sys.exit()
        elif opt in ("-f", "--fBondFile"):
            BondFile = str(arg)
            header = header + 'FILE: -f=' + BondFile
        elif opt in ("-s","--sNsteps"):
            Nsteps = int(arg)
            header = header + ' -s=' + arg
            print('Will sample the MD trajectory every ',Nsteps,' steps')
        elif opt in ("-m","--mMinlife"):
            minlife = float(arg)
        elif opt in ("-n","--nnChunks"):
            nChunks = int(arg)
        elif opt in ("-c","--Cations"):
            header = header + ' -c=' + arg
            Cations = arg.split(",")
            #print ('Cation list is: ',Cations)
        elif opt in ("-a","--Anions"):
            header = header + ' -a=' + arg
            Anions= arg.split(",")
            #print ('Anion list is: ',Anions)
        elif opt in ("-r","--rRings"):
            rings = int(arg)
            if rings == 0:
                print ('Calculation of polymerized coordination polyhedra')
            elif rings > 0:
                print ('Calculation of order '+str(arg)+" polymerized coordination")
            else :
                print ('Undefined calculation')
                print ('-r should be a positive integer')
                sys.exit()

    header = header + ' -r=' + str(rings)

    if not (os.path.isfile(BondFile)):
        print ('the UMD files ',BondFile,' does not exist')            
        sys.exit()

    for ii in range(len(Cations)):
        ClusterAtoms.append(Cations[ii])
    for ii in range(len(Anions)):
        if Anions[ii] not in ClusterAtoms:
            ClusterAtoms.append(Anions[ii])
                  
    (MyCrystal,OctIndexes,TimeStep)=prep_read(BondFile)    

    centralatoms = []
    outeratoms = []
    ligands = []
    for iatom in range(MyCrystal.natom):
        for jatom in range(len(ClusterAtoms)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==ClusterAtoms[jatom]:
                ligands.append(iatom)              #contains the list with the index of the ligand atoms from the 0 ... natom
    for iatom in range(MyCrystal.natom):
        for jatom in range(len(Cations)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==Cations[jatom]:
                centralatoms.append(iatom)         #contains the list with the index of the central atoms from the 0 ... natom
        for jatom in range(len(Anions)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==Anions[jatom]:
                outeratoms.append(iatom)           #contains the list with the index of the coordinating atoms from the 0 ... natom

    print('All ligands are: ',ligands)
    print('Central atoms are :',centralatoms)
    print('Coordinating atoms are :',outeratoms)
    
    if centralatoms == []:
        print("ERROR : element ",Cations," not present in simulation")
        sys.exit()
    if outeratoms == []:
        print("ERROR : element ",Anions," not present in simulation")
        
    readRed=partial(read_snapshotbonds_C,CentMin=centralatoms[0],CentMax=centralatoms[-1],OutMin=outeratoms[0],OutMax=outeratoms[-1],File=BondFile)
    

    with concurrent.futures.ProcessPoolExecutor() as executor :
        Data = list(executor.map(readRed,OctIndexes[:-1],OctIndexes[1:],[step for step in range(len(OctIndexes)-1)])) #Converting bond file into python list

    Bonds = [Data[i][0] for i in range(len(Data))]
    BondsIndexes = [Data[i][1] for i in range(len(Data))]#BondIndexes[x][i] indicates how many atoms, for the snapshot n° x, are listed in Bonds before the ones bound to the atom n° i. 
    
    clusteringRed=partial(clustering, CentMin=centralatoms[0],CentMax=centralatoms[-1],OutMin=outeratoms[0],OutMax=outeratoms[-1],r=rings)
       
    with concurrent.futures.ProcessPoolExecutor() as executor :
        clusters=list(executor.map(clusteringRed,Bonds,BondsIndexes, [step for step in range(len(Bonds))])) #Computes the clusters of atoms for each snapshot separately
    
    
    with concurrent.futures.ProcessPoolExecutor() as executor :#divides the succession of snapshots in chunks then establish the life duration of each cluster of atoms within these chunks
        L=len(clusters)
        Clusters=[clusters[i*int(L/nChunks):(i+1)*int(L/nChunks)] for i in range(nChunks)]
        if int(L/nChunks)!=L/nChunks:
            Clusters+=[clusters[(nChunks)*int(L/nChunks):]]
        populations=list(executor.map(analysis_subtab,Clusters,[Step*int(L/nChunks) for Step in range(nChunks)]))

        
    #Fuses the chunks together and updates the life duration information of each cluster
    population=populations[0]
    for pop in populations[1:]: #for each chunk
        for key in pop:
            if key in population.keys():#if the cluster has been previously seen
                if population[key][-1][-1]==pop[key][0][0]-1:#if the cluster was alive at the end of the previous chung and at the beginning of this one
                    population[key][-1]+=pop[key][0]
                    population[key]+=pop[key][1:]
                else :
                    population[key]+=pop[key] 
            else : 
                population[key]=pop[key]

                
        #Creating the output files        
    FileAll = BondFile[:-4] +'.r=' + str(rings) + '.popul.dat'
    print ('Population will be written in ',FileAll,' file')
    FileStat = BondFile[:-4] + '.r=' + str(rings) + '.stat.dat'
    print ('Statistics will be written in ',FileStat,' file')
    header+="\n"            
                
    print("Writing...")
    dicoNames={}
    dicoTimes={}
    dicoStats={}
    total=0
    for key in population : 
        index=[0 for _ in range(MyCrystal.ntypat)]
        atomsCluster=key.strip('[]').split(', ')
        name=''
        for at in atomsCluster:
            index[MyCrystal.typat[int(at)]]+=1

        for elem in range(len(index)) :
            if index[elem]!=0:
                name+=MyCrystal.elements[elem]+'_'+str(index[elem])

        for life in population[key]:
#            print(life,"tot+=",Nsteps*TimeStep*(life[-1]-life[0]+1),"from ", life[-1]-life[0]+1)
            if name in dicoNames: 
                dicoNames[name].append([key,life[0],life[-1],(life[-1]-life[0]+1)*TimeStep])
                dicoStats[name]['lifetime']+=Nsteps*TimeStep*(life[-1]-life[0]+1)
                total+=Nsteps*TimeStep*(life[-1]-life[0]+1)
            else :
                dicoNames[name]=[[key,life[0],life[-1],(life[-1]-life[0]+1)*TimeStep]]
                dicoStats[name]={'lifetime':Nsteps*TimeStep*(life[-1]-life[0]+1), '#atoms': len(atomsCluster)}
                total+=Nsteps*TimeStep*(life[-1]-life[0]+1)

            if life[0] in dicoTimes :
                dicoTimes[life[0]].append([key,name,life[-1],(life[-1]-life[0]+1)*Nsteps*TimeStep])
            else:
                dicoTimes[life[0]]=[[key,name,life[-1],(life[-1]-life[0]+1)*Nsteps*TimeStep]]
        
    newstring="Formula\tBegin (step)\tEnd(step)\tlifetime (fs)\tcomposition\n"

 #   print(total)

    fa = open(FileAll,'w')
    fa.write(header)
    fa.write(newstring)
    for ii in range(len(Bonds)):
        if ii in dicoTimes:
                #print(dicoTimes[ii])
            for clust in dicoTimes[ii]:
                if clust[3]>minlife :
                    newstring=clust[1]+"\t"+str(ii)+"\t"+str(clust[2])+"\t"+str(clust[3])+"\t"+clust[0]+"\n"
                    fa.write(newstring)
    fa.close()
        
    newstring="Cluster\tTime(fs)\tPercent\tNumber of atoms\n"
        
    fa=open(FileStat,'w')
    fa.write(header)
    fa.write(newstring)
    for clust in dicoStats : 
        newstring=clust+"  \t"+str(dicoStats[clust]['lifetime'])+"\t"+str(float(dicoStats[clust]['lifetime'])/float(total))+"\t"+str(dicoStats[clust]['#atoms'])+"\n"
        fa.write(newstring)
    fa.close()
    
    print ('Population written in file :  ',FileAll)
    print ('Statistics written in file :  ',FileStat)
    
    
    end=time.time()
    print("total duration : ",end-start, " s")

if __name__ == "__main__":
    main(sys.argv[1:])
