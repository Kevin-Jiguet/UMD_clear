#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS
###


import sys,getopt,os.path,itertools
import crystallography as cr
import umd_process as umdp
import math
import time
from functools import partial
import concurrent.futures

global population_parallel
global t

population_parallel=[]
t=[0]


def analysis_subtab(Clusters,Step):
        population_parallel={}
        for Snapshot in Clusters :
        #    print("analysis of snapshot no "+str(Step))
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
    




def neighboring(SnapshotBonds,ligands,iatom,newcluster):
    newcluster.append(iatom)
    iatombonds=SnapshotBonds.get(iatom)
    for jatom in ligands:
        if jatom in iatombonds:
            #print(SnapshotBonds)
            #print('iatom=',iatom)
            #print('jatom=',jatom)
            SnapshotBonds[iatom].remove(jatom)
            SnapshotBonds[jatom].remove(iatom)
            newcluster=neighboring(SnapshotBonds,ligands,jatom,newcluster)
    return newcluster

def clustering(SnapshotBonds,i=-1,ligands=0):
    if (i!=-1):
        print("iteration no "+str(i))
    #print('     clustering: start')
    #print ('ligands:',ligands)
    neighbors = []
    #print("TOTAL #bonds",totNbonds/2)
#    sys.setrecursionlimit(int(totNbonds/2)) #change the recursion limit ; bugs with totbondsN/2 in some cases
    for iatom in ligands :
        if(len(SnapshotBonds.get(iatom)))==0:
            neighbors.append([iatom])

    for iatom in ligands :
        
        newcluster = []
        newcluster1 = neighboring(SnapshotBonds,ligands,iatom,newcluster)
        newcluster1.sort()
        newcluster = [i[0] for i in itertools.groupby(newcluster1)]
        if len(newcluster)>1:
            neighbors.append(newcluster)
#    print("neighbors=",neighbors)
    return neighbors

def clusteringnorings(SnapshotBonds,centralatoms,outeratoms): #compute the list of clusters centered on the atoms from the centeratoms list.
    neighbors = [] #will contain the list of clusters, central atom per central atom.
    for iatom in centralatoms:
        newcluster = [iatom]
        iatombonds=SnapshotBonds.get(iatom)
        for atom in iatombonds :
            if atom in outeratoms :
                newcluster.append(atom) 
        if len(newcluster)>1:
            neighbors.append(newcluster) #newcluster contains the list of atoms of a cluster, the first being the central atom.
    
    return neighbors


def read_bonds(BondFile,s,AllAtoms):
    Bonds=[]
    ff=open(BondFile,'r')
    natom=0
    ntypat=0
    types=[]
    elements=[]
    typat=[]
    
    for _ in range(20) :
        line=ff.readline()
        linestr=line.strip()
        entry=linestr.split()
        if(entry==[]):
            line=ff.readline()
            linestr=line.strip()
            entry=linestr.split()
        if(entry[0]=='natom'):
            natom=int(entry[1])
            typat=[0 for _ in range(natom)]
        if(entry[0]=='ntypat'):
            ntypat=int(entry[1])
            types=[0 for _ in range(ntypat)]
            elements=[0 for _ in range(ntypat)]            
        if(entry[0]=='types'):
            for ii in range(ntypat):
                types[ii]=int(entry[ii+1])
        if(entry[0]=='elements'):
            for ii in range(ntypat):
                elements[ii]=entry[ii+1]
        if(entry[0]=='typat'):
            for iatom in range(natom):                
                typat[iatom]=int(entry[iatom+1])
           
           
    SnapshotBondingTable={}
    niter=0
    while True :
        niter+=1
        line=ff.readline()
        if not line : break
        lstrip=line.strip()
        atomslist=lstrip.split()
        if atomslist!=[] and atomslist[0].isdigit() and int(atomslist[0]) in AllAtoms:
            bondslist=[]
            for atom in atomslist[1:]:
                if int(atom) in AllAtoms :
                    bondslist.append(int(atom))
            SnapshotBondingTable[int(atomslist[0])]=bondslist
        if(atomslist!=[] and atomslist[0]=='end'):
            Bonds.append(SnapshotBondingTable)
            SnapshotBondingTable={}
            
    return Bonds

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
    nChunks=4
    start=time.time()
    try:
        opts, arg = getopt.getopt(argv,"hf:s:c:a:m:r:n:",["fBondFile","sSampling_Frequency", "cCations","aAnions","mMinlife","rRings","nnChunks"])
    except getopt.GetoptError:
        print ('speciation.py -f <bond_filename> -s <Sampling_Frequency> -c <Cations> -a <Anions> -m <MinLife> -r <Rings> -n <nChunks>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('speciation.py program to compute bonding maps and identify speciation')
            print ('speciation.py -f <Bonding_filename> -s <Sampling_Frequency> -l <MaxLength> -c <Cations> -a <Anions> -m <MinLife> -r <Rings> -n <nChunks>')
            print ('  default values: -f bonding.umd.dat -s 1 -m 0 -r 1 -n 4')
            print (' the bond file contains the bonds relations for each snapshot.')
            print (' rings = 1 default, polymerization, all anions and cations bond to each other; rings = 0 only individual cation-anion groups')
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
            header = header + ' -r=' + arg
            if rings == 0:
                print ('Calculation of non-polymerized coordination polyhedra')
            elif rings == 1:
                print ('Calculation of polymerized coordination')
            else :
                print ('Undefined calculation')


    if not (os.path.isfile(BondFile)):
        print ('the UMD files ',BondFile,' does not exist')            
        sys.exit()

    for ii in range(len(Cations)):
        ClusterAtoms.append(Cations[ii])
    for ii in range(len(Anions)):
        if Anions[ii] not in ClusterAtoms:
            ClusterAtoms.append(Anions[ii])
    if rings == 0:
        print('searching for cations',Cations)
        print('surrounded by anions ',Anions)
    else:
        print('all atoms bonding:',ClusterAtoms)
        
                  


                  
#reading the xc art coordinates of the atoms from the UMD file. it uses the read_xcart (i.e. only xcart) function from the umd_process library
    MyCrystal = cr.Lattice()
    AllSnapshots = [cr.Lattice]
    (MyCrystal,AllSnapshots,TimeStep)=umdp.read_xcart(BondFile)
    #in order to allow computation of large polymers
    #we need to set by hand the maximum recursion limit
    #to avoid python erro crash
    #12 is a factor f a maximum reasonable limit of a coordination polyhedron
    sys.setrecursionlimit(12*MyCrystal.natom)
    

    centralatoms = []
    outeratoms = []
    ligands = []
    for iatom in range(MyCrystal.natom):
        for jatom in range(len(ClusterAtoms)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==ClusterAtoms[jatom]:
                ligands.append(iatom)           #contains the list with the index of the ligand atoms from the 0 ... natom
    for iatom in range(MyCrystal.natom):
        for jatom in range(len(Cations)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==Cations[jatom]:
                centralatoms.append(iatom)           #contains the list with the index of the central atoms from the 0 ... natom
        for jatom in range(len(Anions)):
            if MyCrystal.elements[MyCrystal.typat[iatom]]==Anions[jatom]:
                outeratoms.append(iatom)           #contains the list with the index of the coordinating atoms from the 0 ... natom
    print('All ligands are: ',ligands)
    print('Central atoms are :',centralatoms)
    print('Coordinating atoms are :',outeratoms)
    
    Bonds=read_bonds(BondFile,Nsteps,ligands)    
        
    clusters = []

        
    FileAll = BondFile +".alt" +'.r' + str(rings) + '.popul.dat'
    print ('Population will be written in ',FileAll,' file')
    FileStat = BondFile + ".alt"+ '.r' + str(rings) + '.stat.dat'
    print ('Statistics will be written in ',FileStat,' file')
    header+="\n"
    
    if rings == 1:
        with concurrent.futures.ProcessPoolExecutor() as executor :
            clusteringRed=partial(clustering,ligands=ligands)
            clusters=list(executor.map(clusteringRed,Bonds,[i for i in range(len(Bonds))]))
    elif rings == 0:
        with concurrent.futures.ProcessPoolExecutor() as executor :
            clusteringnoringsRed=partial(clusteringnorings,centralatoms=centralatoms,outeratoms=outeratoms)
            clusters=list(executor.map(clusteringnoringsRed,Bonds))       
            
    with concurrent.futures.ProcessPoolExecutor() as executor :
        L=len(clusters)
        Clusters=[clusters[i*int(L/nChunks):(i+1)*int(L/nChunks)] for i in range(nChunks)]
        populations=list(executor.map(analysis_subtab,Clusters,[Step*int(L/nChunks) for Step in range(nChunks)]))
    
    population=populations[0]
    for pop in populations[1:]:
        for key in pop:
            if key in population.keys():
                if population[key][-1][-1]==pop[key][0][0]-1:
                    population[key][-1]+=pop[key][0]
                    population[key]+=pop[key][1:]
                else :
                    population[key]+=pop[key]
            else : 
                population[key]=pop[key]
                    
        
        
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
        newstring=clust+"\t"+str(dicoStats[clust]['lifetime'])+"\t"+str(float(dicoStats[clust]['lifetime'])/float(total))+"\t"+str(dicoStats[clust]['#atoms'])+"\n"
        fa.write(newstring)
    fa.close()
    
    
    end=time.time()
    print("total duration : ",end-start)

if __name__ == "__main__":
    main(sys.argv[1:])


