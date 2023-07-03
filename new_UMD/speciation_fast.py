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
import time
from functools import partial
import concurrent.futures
from os.path import join
import ctypes
import umd_processes_fast as umdpf


current_path=os.path.abspath(__file__)#For all this to work, the files < c_clusters_B.so > and < c_define_bonds.so > must be in the same directory than gofr_umd
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u

clust_lib = ctypes.cdll.LoadLibrary(join(path_new, 'c_clusters_B.so'))
clust_lib.fullclustering.argtypes = [ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
clust_lib.fullclustering.restype = ctypes.POINTER(ctypes.c_int)

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
    

def clustering(SnapshotBonds,SnapshotBondIndexes,step,CentMin,CentMax,OutMin,OutMax,r,Nsteps):
    
    print("calculating species from snapshot n. "+str(step*Nsteps))
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

def main(argv):
    umdpf.headerumd()
    BondFile='bonding.umd.dat'
    Nsteps = 1
    ClusterAtoms = []
    minlife = 5
    Central=""
    Adjacent=""
    rings = 1
    header = ''
    nChunks = 8
    start=time.time()
    try:
        opts, arg = getopt.getopt(argv,"hf:s:c:a:m:r:n:",["fBondFile","sSampling_Frequency", "cCentral","aAdjacent","mMinlife","rRings","nnChunks"])
    except getopt.GetoptError:
        print ('speciation.py -f <bond_filename> -s <Sampling_Frequency> -c <Central atoms> -a <Adjacent atoms> -m <MinLife> -r <Rings> -n <nChunks>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('speciation_fast.py program to identify speciation')
            print ('speciation_fast.py -f <Bonding_filename> -s <Sampling_Frequency> -l <MaxLength> -c <Central atoms> -a <Adjacent atoms> -m <MinLife> -r <Rings> -n <nChunks>')
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
        elif opt in ("-c","--Central"):
            header = header + ' -c=' + arg
            Central = str(arg)
        elif opt in ("-a","--Adjacent"):
            header = header + ' -a=' + arg
            Adjacent = str(arg)
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
        print ('the file ',BondFile,' does not exist')            
        sys.exit()


    CentMin,CentMax,OutMin,OutMax,MyCrystal,Bonds,BondsIndexes,TimeStep = umdpf.read_bonds(BondFile,Central,Adjacent,Nsteps=Nsteps)

    ClusterAtoms=[Central,Adjacent]
    print(ClusterAtoms)

    centralatoms = []
    outeratoms = []
    ligands = []
    for iatom in range(MyCrystal.natom):
        for jatom in range(len(ClusterAtoms)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==ClusterAtoms[jatom]:
                ligands.append(iatom)              #contains the list with the index of the ligand atoms from the 0 ... natom
    for iatom in range(MyCrystal.natom):
        for jatom in range(len(Central)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==Central:
                centralatoms.append(iatom)         #contains the list with the index of the central atoms from the 0 ... natom
        for jatom in range(len(Adjacent)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==Adjacent:
                outeratoms.append(iatom)           #contains the list with the index of the coordinating atoms from the 0 ... natom

    print('All ligands are: ',ligands)
    print('Central atoms are :',centralatoms)
    print('Coordinating atoms are :',outeratoms)
    
    if centralatoms == []:
        print("ERROR : element ",Central," not present in simulation")
        sys.exit()
    if outeratoms == []:
        print("ERROR : element ",Adjacent," not present in simulation")
        
    
    clusteringRed=partial(clustering, CentMin=CentMin,CentMax=CentMax,OutMin=OutMin,OutMax=OutMax,r=rings,Nsteps=Nsteps)
       
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
                if population[key][-1][-1]==pop[key][0][0]-1:#if the cluster was alive at the end of the previous chunk and at the beginning of this one
                    population[key][-1]+=pop[key][0]
                    population[key]+=pop[key][1:]
                else :
                    population[key]+=pop[key] 
            else : 
                population[key]=pop[key]

                
        #Creating the output files        
    FileAll = BondFile[:-4] +'.r=' + str(rings) + '.populB.dat'
#    print ('Population will be written in ',FileAll,' file')
    FileStat = BondFile[:-4] + '.r=' + str(rings) + '.statB.dat'
#    print ('Statistics will be written in ',FileStat,' file')
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
                    newstring=clust[1]+"\t"+str(ii*Nsteps)+"\t"+str(clust[2]*Nsteps)+"\t"+str(clust[3])+"\t"+clust[0]+"\n"
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
