#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 09:44:43 2023
"""
#!/usr/bin/env python3
###
##AUTHORS: KEVIN JIGUET
###

import platform
import sys,getopt,os.path
import umd_processes_fast as umdpf
import time
from functools import partial
import concurrent.futures
from os.path import join
import ctypes
import numpy as np


current_path=os.path.abspath(__file__)#For all this to work, the files < c_clusters_B.so > and < c_define_bonds.so > must be in the same directory than gofr_umd
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u

OS = platform.system()
LibraryName =""


if OS == "Linux":
    LibraryName = 'c_clusters_full.so'

elif OS == "Windows":
    LibraryName = 'c_clusters_full.dll'

elif OS == "Darwin":
    LibraryName = 'c_clusters_full.dylib'
    


clust_lib = ctypes.cdll.LoadLibrary(join(path_new, LibraryName))
clust_lib.fullclustering.argtypes = [ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
clust_lib.fullclustering.restype = ctypes.POINTER(ctypes.c_int)

clust_lib.angles.argtypes = [ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double]
clust_lib.angles.restype = ctypes.POINTER(ctypes.c_double)

clust_lib.free_memory_int.argtypes = [ctypes.POINTER(ctypes.c_int)]
clust_lib.free_memory_int.restype = None

clust_lib.free_memory_double.argtypes = [ctypes.POINTER(ctypes.c_double)]
clust_lib.free_memory_double.restype = None


def analysis_subtab(Clusters,Step):#Creates a dictionnary whose keys are the individual clusters, and the values are a list of list of their *consecutive* steps of existence. 
        population={}
        for Snapshot in Clusters :
            for cluster in Snapshot :
                if str(cluster) in population.keys():
                    if population[str(cluster)][-1][-1]==Step-1:
                        population[str(cluster)][-1].append(Step)
                    else :
                        population[str(cluster)].append([Step])
                else :
                    population[str(cluster)]=[[Step]]
            Step+=1
        return population
    

def clustering(SnapshotBonds,SnapshotBondIndexes,SnapshotXCart,step,natom,nAts,CentIndexes,OutIndexes,acell,r,AngleCalc):
    
#    print(len(SnapshotBonds),len(SnapshotBondIndexes),len(SnapshotXCart))
    print("calculating species from snapshot n. "+str(step))
    M=SnapshotBondIndexes[-1]#maximum number of bound atoms for any atoms
    #Preparing data for the C script        

    SBp = (ctypes.c_int * len(SnapshotBonds))(*SnapshotBonds)
    BIp = (ctypes.c_int * (len(SnapshotBondIndexes)-1))(*SnapshotBondIndexes[:-1])
    CIp = (ctypes.c_int * len(CentIndexes))(*CentIndexes)
    OIp = (ctypes.c_int * len(OutIndexes))(*OutIndexes)
    
    Np = clust_lib.fullclustering(SBp,BIp,natom,nAts,CIp,OIp,int(len(CentIndexes)/2),int(len(OutIndexes)/2),M,r)#Computes the clusters
    Clusters = [[]]
    Angles = {}
    length = Np[0]
    #Converting C data into a python list of clusters
    for i in range(1,length+1):
        atom=Np[i]
#        print(atom)
        if atom ==-1 :
#            Clusters[-1].sort()
            Clusters.append([])
        else :   
            Clusters[-1].append(atom)

    
    if r==1 and AngleCalc:
        SXp = (ctypes.c_double * len(SnapshotXCart))(*SnapshotXCart)
        Ap = clust_lib.angles(Np,len(Clusters),M,SXp,CIp,int(len(CentIndexes)/2),OIp,int(len(OutIndexes)/2),acell[0],acell[1],acell[2])
        index=0
        flagnew=1
#        sys.exit()
#        print(Ap[0],M,len(Clusters))
        for i in range(1,int(Ap[0])):
  #          print(Ap[i],Np[i])
            if(Ap[i]==-1):
                index=index+1
                flagnew=1
            elif flagnew==1:
                flagnew=0
                Angles[str(Clusters[index])]=[Ap[i]]
            else :
                Angles[str(Clusters[index])].append(Ap[i])

        clust_lib.free_memory_double(Ap)    

    clust_lib.free_memory_int(Np)

    if Clusters ==[[]]:
        return [],{}
    
    return Clusters,Angles
    

def is_in(at,Indexes):
    for i in range(int(len(Indexes)/2)):
        if (at<=Indexes[2*i+1] and at>=Indexes[2*i]):
            return True
    return False

def main(argv):
    umdpf.headerumd()
    BondFile='bonding.umd.dat'
    UMDFile=''
    Nsteps = 1
    Centrals = []
    Adjacents = []
    minlife = 5
    rings = 1
    header = ''
    start=time.time()
    t = False
    try:
        opts, arg = getopt.getopt(argv,"hf:u:s:c:a:m:r:t:",["fBondFile","uUMDFile","sSampling_Frequency", "cCentral","aAdjacent","mMinlife","rRings","tAngles"])
    except getopt.GetoptError:
        print ('speciation_and_angles.py -f <bond_filename> -s <Sampling_Frequency> -c <Cations> -a <Anions> -m <MinLife> -r <Rings> -n <nChunks>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('speciation_and_angles.py program to compute bonding maps and identify speciation')
            print ('speciation_and_angles.py -f <Bonding_filename> -s <Sampling_Frequency> -c <Central> -a <Adjacent> -m <MinLife> -r <Rings> -t <Angles>')
            print ('default values: -f bonding.umd.dat -s 1 -m 5 -r 1')
            print ('the bond file contains the bonds relations for each snapshot. Computed with Bond_fast_specific.py.')
            print ("-c and -a : central and adjacent elements respectively. If one is 'all', every atom will be taken in account for this role.")
            print ('-r : rings = 1 default, all adjacent atoms bind to central atoms ; rings = 0, polymerization, all adjacent atoms bind to central AND other adjacent atoms ; rings = x>0, all adjacent atoms bind to central then to other adjacent atoms to form a xth-coordination polyhedra')
            print ('-m : minimal duration of existence for a chemical species to be taken into account (fs) ; default 5')
            print ('-t : nothing (default) or True. In the latter case, and only if -r = 1,the angles of the molecules will be computed, their summit being the central atom. If this option is activated, the user needs to provide the corresponding umd file.')
            print ('-u : umd file used to calculate the angles.')
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
        elif opt in ("-u","uUMDFile"):
            UMDFile = str(arg)
        elif opt in ("-c","--Central"):
            header = header + ' -c=' + arg
            Centrals = arg.split(",")
            #print ('Cation list is: ',Cations)
        elif opt in ("-a","--Adjacent"):
            header = header + ' -a=' + arg
            Adjacents = arg.split(",")
            #print ('Anion list is: ',Anions)
        elif opt in("-t","--Angles"):
            t = str(arg)
            if t == "True":
                t=True
        elif opt in ("-r","--rRings"):
            rings = int(arg)
            if rings == 0:
                print ('Calculation of polymerized coordination polyhedra')
            elif rings > 0:
                print ('Calculation of order '+str(arg)+" coordination sphere")
            else :
                print ('Undefined calculation')
                print ('-r should be a positive integer')
                sys.exit()

    header = header + ' -r=' + str(rings)

    if not (os.path.isfile(BondFile)):
        print ('the bond file ',BondFile,' does not exist')            
        sys.exit()

    if t and (not (os.path.isfile(UMDFile))):
        print ('the UMD file ',UMDFile,' does not exist')            
        sys.exit()


    CentIndexes,AdjIndexes,MyCrystal,[Bonds,BondsIndexes],TimeStep = umdpf.read_bonds(BondFile,Centrals,Adjacents,Nsteps)
    
    if CentIndexes == -1 :
        return AdjIndexes
    
    AllElements = []
    nAts = 0
    if Centrals != ["all"]:
        for el in Centrals : 
            if el not in AllElements : 
                AllElements.append(el)
                nAts += MyCrystal.types[MyCrystal.elements.index(el)]
    else : 
        AllElements = MyCrystal.elements
        nAts = MyCrystal.natom
    
    if Adjacents != ["all"]:
        for el in Adjacents :
            if el not in AllElements : 
                AllElements.append(el)
                nAts += MyCrystal.types[MyCrystal.elements.index(el)]

    else :
        AllElements = MyCrystal.elements
        nAts = MyCrystal.natom




    if t and rings==1:
        MyCrystalUMD,TimeStepUMD = umdpf.Crystallization(UMDFile)
        TimeRatio = int(TimeStep/TimeStepUMD)
        MyCrystalUMD,SnapshotsXCart,TimeStepUMD,length = umdpf.read_values(UMDFile,"xcart","line",Nsteps*TimeRatio)
    else :
        SnapshotsXCart = [[] for _ in range(len(Bonds))]#Blank list which will not be used but is needed as an argument
    
    
#    CentMax,AdjMax,CentMin,AdjMin = 0,0,0,0
    clusteringRed=partial(clustering,natom=MyCrystal.natom,nAts=nAts,CentIndexes=CentIndexes,OutIndexes=AdjIndexes,r=rings,acell=MyCrystal.acell,AngleCalc=t)


#    DataII = []
#    for i in range(len(Bonds)):
#        DataII.append(clusteringRed(Bonds[i],BondsIndexes[i],SnapshotsXCart[i],i))
#    print(DataII[0])
#    sys.exit()
    with concurrent.futures.ProcessPoolExecutor() as executor :
        DataII=list(executor.map(clusteringRed,Bonds,BondsIndexes,SnapshotsXCart, [step for step in range(len(Bonds))])) #Computes the clusters of atoms for each snapshot separately

    clusters,Angles = map(list,zip(*DataII))
    population=analysis_subtab(clusters,0)
        #Creating the output files        
    FileAll = BondFile[:-4] +'.r=' + str(rings) + '.popul.dat'
    #print ('Population will be written in ',FileAll,' file')
    FileStat = BondFile[:-4] + '.r=' + str(rings) + '.stat.dat'
    #print ('Statistics will be written in ',FileStat,' file')
    FileStep = BondFile[:-4] + '.r=' + str(rings) + '.step.dat'
    header+="\n"            
    
    print("Writing...")
    
    fs = open(FileStep,'w')
    fs.write(header)        
    headstring = "Clusters\tNumber of atoms\tComposition\n"
    for step in range(len(clusters)) :
        st = "step "+str(step*Nsteps)+"\ntime "+str(step*Nsteps*TimeStep)+"\n"
        fs.write(st)
        fs.write(headstring)
        Clusts = clusters[step]
        vapor = [at for at in range(MyCrystal.natom)]

        for clust in Clusts :
            index=[0 for _ in range(MyCrystal.ntypat)]
            name=''
            for at in clust:
                index[MyCrystal.typat[at]]+=1

            for elem in range(len(index)) :
                if index[elem]!=0:
                    name+=MyCrystal.elements[elem]+'_'+str(index[elem])
            newstring = name +'\t'+ str(len(clust)) +'\t'+ str(clust)+'\n'
            fs.write(newstring)
            for at in clust :
                vapor[at]=-1
        for at in vapor :
            if at !=-1 and (is_in(at,CentIndexes) or is_in(at,AdjIndexes)) :
                fs.write(MyCrystal.elements[MyCrystal.typat[at]]+'\t1\t'+str(at)+'\n')
        
    fs.close()
            
                    
   
    dicoNames={}
    dicoTimes={}
    dicoStats={}
    dicoMeanAngle={}
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
            if (Nsteps*TimeStep*(life[-1]-life[0]+1))>minlife :
                
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
    
    if rings==1 and t:
        newstring="Formula\tBegin (step)\tEnd(step)\tlifetime (fs)\tcomposition\tAngles\n"
    else :
        newstring="Formula\tBegin (step)\tEnd(step)\tlifetime (fs)\tcomposition\n"
     

    
    
    
    
    fa = open(FileAll,'w')
    fa.write(header)
    fa.write(newstring)
    for ii in range(len(Bonds)):
        if ii in dicoTimes:
                #print(dicoTimes[ii])
            for data in dicoTimes[ii]:
                if data[3]>minlife :
                    newstring=data[1]+"\t"+str(ii*Nsteps)+"\t"+str(data[2]*Nsteps)+"\t"+str(data[3])+"\t"+data[0]
                    
                    if rings==1 and t :
                        if data[0] in Angles[ii].keys():
                            angles=np.array([np.array(Angles[jj][data[0]]) for jj in range(ii,data[2]+1)])
                            Meanangles = np.sum(angles,0)/len(angles)
                            totmean = np.sum(Meanangles,0)/len(Meanangles)

                            if data[1] in dicoMeanAngle:
                                dicoMeanAngle[data[1]][0]+=totmean
                                dicoMeanAngle[data[1]][1]+=1
                            else:
                                dicoMeanAngle[data[1]]=[totmean,1]
                        else:
                            Meanangles=np.array([])
                    
                        for alpha in Meanangles:
                            newstring+="\t"+str(alpha)
                    
                    newstring+="\n"
                    fa.write(newstring)
    fa.close()

        
    if(rings==1 and t):
        newstring="Cluster\tTime(fs)\tPercent\tNumber of atoms\tMean Angles\n"
    else:
        newstring="Cluster\tTime(fs)\tPercent\tNumber of atoms\n"    
    
    fa=open(FileStat,'w')
    fa.write(header)
    fa.write(newstring)
    for clust in dicoStats : 
        newstring=clust+"  \t"+str(dicoStats[clust]['lifetime'])+"\t"+str(float(dicoStats[clust]['lifetime'])/float(total))+"\t"+str(dicoStats[clust]['#atoms'])
        if rings == 1 and clust in dicoMeanAngle:
            newstring+='\t'+str(dicoMeanAngle[clust][0]/dicoMeanAngle[clust][1])
        newstring+='\n'
        fa.write(newstring)
    fa.close()
    
    print ('Population written in file :  ',FileAll)
    print ('Statistics written in file :  ',FileStat)
    print ('Species for each snapshot written in file :',FileStep)
    
    end=time.time()
    print("total duration : ",end-start, " s")

if __name__ == "__main__":
    main(sys.argv[1:])