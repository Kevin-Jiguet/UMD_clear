#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 10:04:14 2023

@author: kevin Jiguet
"""

import sys,getopt,os.path,itertools,math
import crystallography as cr
import umd_process as umdp


def compute_angles(Bonds,MyCrystal,AllSnapshots,internatoms,externatoms,fa,TimeStep):
    acell=MyCrystal.acell
    TotMean=0
    TotNangles=0
    for step in range(len(Bonds)) :
        print("step ",step," on",len(Bonds))
        SnapshotAngles=[]
        SnapshotBonds=Bonds[step]
        MySnapshot=AllSnapshots[step]
        fa.write('\nstep : '+str(step)+'\ntime :'+str(step*TimeStep)+" fs \n")
        line='Angles : '
        Nangles=0
        for internAt in SnapshotBonds :
            for iexternAt in SnapshotBonds[internAt]:
                for jexternAt in SnapshotBonds[internAt]:
                    if iexternAt<jexternAt :    
                        Nangles+=1
                        
                        xEi,yEi,zEi=MySnapshot.atoms[iexternAt].xcart[0],MySnapshot.atoms[iexternAt].xcart[1],MySnapshot.atoms[iexternAt].xcart[2]
                        xEj,yEj,zEj=MySnapshot.atoms[jexternAt].xcart[0],MySnapshot.atoms[jexternAt].xcart[1],MySnapshot.atoms[jexternAt].xcart[2]
                        xI,yI,zI=MySnapshot.atoms[internAt].xcart[0],MySnapshot.atoms[internAt].xcart[1],MySnapshot.atoms[internAt].xcart[2]
                        
                        dx = xEi-xEj
                        dy = yEi-yEj
                        dz = zEi-zEj
                        
                        valx = min(dx**2, (acell[0]-dx)**2, (acell[0]+dx)**2)
                        valy = min(dy**2, (acell[1]-dy)**2, (acell[1]+dy)**2)
                        valz = min(dz**2, (acell[2]-dz)**2, (acell[2]+dz)**2)
                        distEiEj = valx + valy + valz               
                        
                        dx = xI-xEj
                        dy = yI-yEj
                        dz = zI-zEj                        
                        
                        valx = min(dx**2, (acell[0]-dx)**2, (acell[0]+dx)**2)
                        valy = min(dy**2, (acell[1]-dy)**2, (acell[1]+dy)**2)
                        valz = min(dz**2, (acell[2]-dz)**2, (acell[2]+dz)**2)
                        distIEj = valx + valy + valz               
                        
                        dx = xEi-xI
                        dy = yEi-yI
                        dz = zEi-zI                        
                        
                        valx = min(dx**2, (acell[0]-dx)**2, (acell[0]+dx)**2)
                        valy = min(dy**2, (acell[1]-dy)**2, (acell[1]+dy)**2)
                        valz = min(dz**2, (acell[2]-dz)**2, (acell[2]+dz)**2)
                        distIEi = valx + valy + valz               
        
                        angle=math.acos((distIEi+distIEj-distEiEj)/(2*math.sqrt(distIEi*distIEj)))
                        
                        line+="\t"+str(angle)+" ("+str(iexternAt)+"-"+str(internAt)+"-"+str(jexternAt)+")"
                        
                        SnapshotAngles.append(angle)
        fa.write(line)
        fa.write("\n")

        TotNangles+=Nangles
        if Nangles!=0:
            mean=sum(SnapshotAngles)/Nangles                  
            msd=0
            TotMean+=sum(SnapshotAngles)
            for angle in SnapshotAngles:
                msd+=(angle-mean)**2
            msd/=mean
            msd=math.sqrt(msd)
        
            fa.write("\t"+"mean : "+str(mean)+"\tmsd : "+str(msd)+"\n")
        fa.write("end of step")

    if(TotNangles!=0):
        TotMean/=TotNangles

    fa.write("\n\n Mean of all angles : "+str(TotMean))
    fa.close()        

def read_bonds_angles(BondFile,s,centralatoms,externatoms):
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
        if atomslist!=[] and atomslist[0].isdigit() and int(atomslist[0]) in centralatoms:
            bondslist=[]
            for atom in atomslist[1:]:
                if int(atom) in externatoms :
                    bondslist.append(int(atom))
            if(len(bondslist)>1):
                SnapshotBondingTable[int(atomslist[0])]=bondslist
#                print(bondslist)
        if(atomslist!=[] and atomslist[0]=='end'):
            Bonds.append(SnapshotBondingTable)
            SnapshotBondingTable={}
            
    return Bonds


def main(argv):
    umdp.headerumd()
    Nsteps = 1
    header = ''
    try:
        opts, arg = getopt.getopt(argv,"hf:b:s:c:e:m:",["fUMDFile","bBondFile","sSampling_Frequency","cCentral_Atoms","eExternAtoms"])
    except getopt.GetoptError:
        print ('Angles.py -f <UMD_filename> -b <Bond_filename> -s <Sampling_Frequency> -c <CentralAtom> -e <ExterAtoms>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Angles.py program to compute angles between atoms')
            print ('Angles.py -f <UMD_filename> -b <BondFile> -c <innerAtoms> -e <externAtoms>')
            print (' the bond file contains the bonds for each atom. It is computed using the script Bond.py')
            sys.exit()
        elif opt in ("-f", "--fUMDFile"):
            UMDname = str(arg)
            header = header + 'FILE: -f=' + UMDname
        elif opt in ("-b","--bBondFile"):
            BondFile = str(arg)
            header = header + ' -b=' + BondFile
        elif opt in ("-s","--sNsteps"):
            Nsteps = int(arg)
            header = header + ' -s=' + arg
            print('Will sample the data every ',Nsteps,' steps')        
        elif opt in ("-c","--cCentralAtom"):
            cen=str(arg)
        elif opt in ("-e","--eExternAtom"):
            ext=str(arg)
            header=header+"\t angles : "+ext+"-"+cen+"-"+ext

    MyCrystal = cr.Lattice()
    AllSnapshots = [cr.Lattice]
    (MyCrystal,AllSnapshots,TimeStep)=umdp.read_xcart(UMDname)
    
    inneratoms=[]
    outeratoms=[]
    
    for iatom in range(MyCrystal.natom):
        elem=MyCrystal.elements[MyCrystal.typat[iatom]]
        if elem==cen:
            inneratoms.append(iatom)
        if elem==ext:
            outeratoms.append(iatom)
    
        
    Bonds=read_bonds_angles(BondFile,Nsteps,inneratoms,outeratoms)    
    
    File=str(UMDname)+".angles.dat"
    print(File)
    
    fa=open(File,'w')
    fa.write(header)
    compute_angles(Bonds,MyCrystal,AllSnapshots,inneratoms,outeratoms,fa,TimeStep)    


if __name__=='__main__':
   main(sys.argv[1:])
