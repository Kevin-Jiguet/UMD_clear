#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Created on Mon Apr 24 15:20:53 2023

#author: kevin Jiguet

import sys,getopt,os.path,itertools
import crystallography as cr
import umd_process as umdp
import math
import time
from functools import partial
import concurrent.futures

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
    return BondTable



    
    
def  WriteBonding(MySnapshot,step,MyCrystal,BondTable,timestep,natom,numCells):
    print("Step :",step )
    Mesh=[[[[] for _ in range(numCells)]for _ in range(numCells)] for _ in range(numCells)]

    Places={}
    
    Lx=MyCrystal.acell[0]/numCells
    Ly=MyCrystal.acell[1]/numCells
    Lz=MyCrystal.acell[2]/numCells
            
    for iatom in range(natom):
        x=math.fmod(MySnapshot.atoms[iatom].xcart[0],MyCrystal.acell[0])
        y=math.fmod(MySnapshot.atoms[iatom].xcart[1],MyCrystal.acell[1])
        z=math.fmod(MySnapshot.atoms[iatom].xcart[2],MyCrystal.acell[2])
        xCell=int(x/Lx)
        yCell=int(y/Ly)
        zCell=int(z/Lz)
                
        Mesh[xCell][yCell][zCell].append(iatom)
        Places[iatom]=[xCell,yCell,zCell]

    lines=[[at] for at in range(natom)]
       
    
    
    
    for iatom in range(natom):
#        print("looking for atom",iatom)
        if Places[iatom]!=None:
        
            [xCell,yCell,zCell]=Places[iatom]
            Atoms=Mesh[xCell][yCell][zCell]
            PotentialNeighbors=[]
            
            l=[-1,0,1]
        
            
            for i in l:
                for j in l:
                    for k in l:                    
                        PotentialNeighbors+=list.copy(Mesh[(xCell+i)%numCells][(yCell+j)%numCells][(zCell+k)%numCells])

            PotentialNeighbors.sort()
            PotentialNeighbors=[i[0] for i in itertools.groupby(PotentialNeighbors)]
            
            for katom in Atoms :#Looking at every atom in the cell
                Places[katom]=None#removing each atom from the to-do list
                for jatom in PotentialNeighbors:
                    if(katom<jatom):
                        dx = MySnapshot.atoms[jatom].xcart[0] - MySnapshot.atoms[katom].xcart[0]
                        dy = MySnapshot.atoms[jatom].xcart[1] - MySnapshot.atoms[katom].xcart[1]
                        dz = MySnapshot.atoms[jatom].xcart[2] - MySnapshot.atoms[katom].xcart[2]
                        
                        valx = min(dx**2, (MyCrystal.acell[0]-dx)**2, (MyCrystal.acell[0]+dx)**2)
                        valy = min(dy**2, (MyCrystal.acell[1]-dy)**2, (MyCrystal.acell[1]+dy)**2)
                        valz = min(dz**2, (MyCrystal.acell[2]-dz)**2, (MyCrystal.acell[2]+dz)**2)
                        distkj = valx + valy + valz
                        if distkj<BondTable[MyCrystal.typat[katom]][MyCrystal.typat[jatom]]:
                            #print("distance = ",distkj, "on ",BondTable[MyCrystal.typat[katom]][MyCrystal.typat[jatom]])
                            lines[katom].append(jatom)
                            lines[jatom].append(katom)
                            
    return lines
    

        

def main(argv):
    umdp.headerumd()
    UMDname='output.umd.dat'
    Nsteps = 1
    InputFile = ''
    header = ''
    numCells=5
    try:
        opts, arg = getopt.getopt(argv,"hf:s:l:i:n:",["fUMDfile","sSampling_Frequency","lMaxLength","iInputFile","nNumCells"])
    except getopt.GetoptError:
        print ('speciation.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -i <InputFile> -numCells')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('computation of bonding map')
            print ('bond.py -f <UMD_filename> -s <Sampling_Frequency> -l <MaxLength> -c <Cations> -a <Anions> -m <MinLife>  -i <InputFile> -r <Rings>')
            print ('  default values: -f output.umd.dat -s 1 -l 3.0 -m 0 -r 1 -n 5')
            print (' the input file contains the bond lengths for the different atom pairs. \n the values overwrite the option -l')
            print (' rings = 1 default, polymerization, all anions and cations bond to each other; rings = 0 only individual cation-anion groups')
            sys.exit()
        elif opt in ("-l","--lmaxlength"):
            maxlength = float(arg)    
            header = header + "-l=" + str(maxlength)
        elif opt in ("-f", "--fUMDfile"):
            UMDname = str(arg)
            header = header + 'FILE: -f=' + UMDname
            maxlength=None
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

    MyCrystal = cr.Lattice()
    AllSnapshots = [cr.Lattice]
    (MyCrystal,AllSnapshots,TimeStep)=umdp.read_xcart(UMDname)       


    
    if maxlength==None and len(InputFile)>0 :
        BondTable = read_inputfile(InputFile,MyCrystal)
    elif maxlength!=None :
        BondTable = [[maxlength for _ in range(MyCrystal.ntypat)] for _ in range(MyCrystal.ntypat)]


    M=max([max(Bondlengths) for Bondlengths in BondTable])

    if MyCrystal.acell[0]/numCells<M or MyCrystal.acell[1]/numCells<M or MyCrystal.acell[2]/numCells<M:
        print('WARNING : dimension of the cell smaller than the greatest bondlength.')
    
        

    natom=MyCrystal.natom
    FileAll=UMDname+'.bondingII.dat'
    print ('Bondings will be written in ',FileAll,' file')
    fa=open(FileAll,'w')
    ff=open(UMDname,'r')
    
    for i in range(20):    
        line=ff.readline()
        fa.write(line)
        
    ff.close()
       
    with concurrent.futures.ProcessPoolExecutor() as executor :
        WriteBondingRed=partial(WriteBonding,MyCrystal=MyCrystal,BondTable=BondTable,timestep=TimeStep,natom=natom,numCells=numCells)
        Lines=list(executor.map(WriteBondingRed,AllSnapshots,[step for step in range(0,len(AllSnapshots),Nsteps)]))

        #print(Lines)
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
