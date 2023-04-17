#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS
###
#This script calculates the polymers or coordination polyhedra of the atoms at each time. The Bondfile containing the bonds data is produced by the script Bond.py.

import sys,getopt,os.path,itertools
import crystallography as cr
import umd_process as umdp

def analysis_clusters(clusters,MyCrystal,ligands,minlife,Nsteps,FileName,rings,TimeStep):
    #this functions builds the dictionary with all the clusters
    population = {}
    for ii in range(len(clusters)):
        print ('snapshot no. ',ii,"on",len(clusters))
        #labels.append([[' ' for _ in range(2)] for _ in range(len(clusters[ii]))])
        for jj in range(len(clusters[ii])):
            word = '' #Will contain which atom of the crystal is in a given cluster
            molecule = '' #Will contain the name of the cluster (molecule style)
            indexes = [0 for _ in range(len(MyCrystal.elements))] #Will contain how much of each elements are presents in a given cluster
            for kk in range(len(clusters[ii][jj])):
                word = word + str(clusters[ii][jj][kk]) + '-'
                indexes[MyCrystal.typat[clusters[ii][jj][kk]]] +=1
            for kk in range(len(MyCrystal.elements)):
                if indexes[kk]>0:
                    molecule = molecule + str(MyCrystal.elements[kk]) + '_' + str(indexes[kk])
            clustername = molecule + '_' + word
            clusterindex = clustername + '_' + str(ii) # name + the index of the snapshot
            #print ('     ===>  current cluster to parse is ',clustername,' with index',clusterindex)
            flagalive = 0
                #if len(population)>1:
            for ll in population.keys():
                #print ('comparing currect cluster ',population[ll]['clustername']) #,' with ',population[kk][clustername])
                if clustername == population[ll]['clustername']:
                    #print ('current step',ii,'end of population',population[ll]['end'])
                    if ii - population[ll]['end'] == 1: #If the cluster already existed in the previous snapshot, we increase its length of life by 1
                        population[ll]['end'] += 1
                        population[ll]['lifetime'] += 1
                        
                        flagalive = 1
            if flagalive == 0:
                population[clusterindex]={'clustername':clustername, 'formula':molecule, 'word':word, 'begin':ii, 'end':ii, 'lifetime':1,'composition':clusters[ii][jj]}
    FileAll = FileName + '.r' + str(rings) + '.popul.dat'
    #print ('Population will be written in ',FileAll,' file')
    fa = open(FileAll,'a')
    newstring = 'Formula\t\tBegin(step)\tEnd(step)\tLifetime(fs)\t[Composition]\n'
    fa.write(newstring)
    for kk in population.keys():
        if population[kk]['lifetime'] > minlife/float(Nsteps):
            newstring = population[kk]['formula'] + '\t' + str(population[kk]['begin']*Nsteps) + '\t' + str(population[kk]['end']*Nsteps) +'\t' + str(population[kk]['lifetime']*Nsteps*TimeStep) + '\t' + str(population[kk]['composition']) + '\n'
            fa.write(newstring)
    
    statclusters = [['',0,0]]
    #print('length of statclusters is',len(statclusters))
    for kk in population.keys():
        #print('treating cluster',population[kk]['formula'])
        ll = 0
        flagnewclust = 0
        while flagnewclust == 0:
            #print('comparing to cluster',statclusters[ll][0])
            if population[kk]['formula'] == statclusters[ll][0]:
                statclusters[ll][1] += population[kk]['lifetime']*Nsteps*TimeStep
                #print('current population is',statclusters[ll][1])
                flagnewclust = 1
            else:
                if ll == len(statclusters)-1:
                    #print('adding new cluster',population[kk]['formula'])
                    statclusters.append([population[kk]['formula'],population[kk]['lifetime']*Nsteps*TimeStep,len(population[kk]['composition'])])
                    flagnewclust = 1
                else:
                    ll +=1

    totalpop = 0
    for ll in range(1,len(statclusters)):
        totalpop += statclusters[ll][1]
    FileStat = FileName + '.r' + str(rings) + '.stat.dat'
#print ('Statistics will be written in ',FileStat,' file')
    fs = open(FileStat,'a')
    newstring = 'Cluster\tTime(fs)\tPercent\n'
    fs.write(newstring)
    for ll in range(1,len(statclusters)):
        newstring = statclusters[ll][0] + '\t' + str(statclusters[ll][1]) + '\t' + str(float(statclusters[ll][1])/float(totalpop)) + '\t' + str(statclusters[ll][2]) + '\n'
        fs.write(newstring)



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

def clustering(SnapshotBonds,ligands):
    #print('     clustering: start')
    #print ('ligands:',ligands)
    neighbors = []
    #print("TOTAL #bonds",totNbonds/2)
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
    
#    print(natom,ntypat,types,elements,typat)
       
           
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
#                print(bondslist)
        if(atomslist!=[] and atomslist[0]=='end'):
            Bonds.append(SnapshotBondingTable)
            SnapshotBondingTable={}
            
    return Bonds

def main(argv):
    umdp.headerumd()
    BondFile='output.umd.dat'
    Nsteps = 1
    ClusterAtoms = []
    Cations = []
    Anions = []
    minlife = 5
    rings = 1
    header = ''
    try:
        opts, arg = getopt.getopt(argv,"hf:s:c:a:m:r:",["fBondFile","sSampling_Frequency","cCations","aAnions","mMinlife","rRings"])
    except getopt.GetoptError:
        print ('speciation.py -f <UMD_filename> -s <Sampling_Frequency> -c <Cations> -a <Anions> -m <MinLife> -i <InputFile> -r <Rings>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('speciation_alt.py program to compute bonding maps and identify speciation')
            print ('speciation_alt.py -f <UMD_filename> -s <Sampling_Frequency> -c <Cations> -a <Anions> -m <MinLife>  -i <InputFile> -r <Rings>')
            print ('default values: -f output.umd.dat -s 1 -m 0 -r 1')
            print ('the Bond file contains the bond data for each atom pair at each time.')
            print ('rings = 1 default, polymerization, all anions and cations bond to each other; rings = 0 only individual cation-anion groups')
            sys.exit()
        elif opt in ("-f", "--fBondFile"):
            BondFile = str(arg)
            print(BondFile)
            header = header + 'FILE: -b=' + BondFile
        elif opt in ("-s","--sNsteps"):
            Nsteps = int(arg)
            header = header + ' -s=' + arg
            print('Will sample the MD trajectory every ',Nsteps,' steps')
        elif opt in ("-m","--mMinlife"):
            minlife = float(arg)
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
        
                  
#writes the header of the files containing eventually the clusters
    header = header + '\n'
    FileAll = BondFile + '.r' + str(rings) + '.popul.dat'
    print ('Population will be written in ',FileAll,' file')
    fa = open(FileAll,'w')
    fa.write(header)
    fa.close()
    FileStat = BondFile + '.r' + str(rings) + '.stat.dat'
    print ('Statistics will be written in ',FileStat,' file')
    fa = open(FileStat,'w')
    fa.write(header)
    fa.close()

                  
#reading the xc art coordinates of the atoms from the UMD file. it uses the read_xcart (i.e. only xcart) function from the umd_process library
    MyCrystal = cr.Lattice()
    AllSnapshots = [cr.Lattice]
    (MyCrystal,AllSnapshots,TimeStep)=umdp.read_xcart(BondFile)
    #print('checks after reading the umd file')
    #print('no of atoms = ',MyCrystal.natom)

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
#    print("Bonds :",Bonds)
#span the entire trajectory and analyze only the Nsteps snapshots
    clusters = []
    for istep in range(0,len(Bonds),Nsteps):
        #print('analyzing new simulation snapshot, step no. ',istep)
#first build the bonding maps
        print("step ",istep," on ",len(Bonds))
        if rings == 1:
            clusters.append(clustering(Bonds[istep],ligands))
        elif rings == 0:
            clusters.append(clusteringnorings(Bonds[istep],centralatoms,outeratoms)) #Contain the list of all the clusters for each snapshot
        #print('number of identified clusters ',len(clusters))
        #print(' clusters at this point:\n',clusters)
        else:
            print ('value of rings = ',rings,' is not allowed. Only 0 (polymers) or 1 (coordinating polyhedra) are allowed',)
            sys.exit()
    analysis_clusters(clusters,MyCrystal,ligands,minlife,Nsteps,BondFile,rings,TimeStep)
#    print(clusters)



if __name__ == "__main__":
    main(sys.argv[1:])


