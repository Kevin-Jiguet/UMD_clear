#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 10:04:14 2023

@author: k
"""

import sys,getopt,math
import crystallography as cr
import umd_processes_fast as umdpf


def compute_angles(Bonds,MyCrystal,AllSnapshots,CentMin,CentMax,OutMin,OutMax,fa,TimeStep,Nsteps):
    acell=MyCrystal.acell
    TotMean=0
    TotNangles=0
    for step in range(len(Bonds)) :
        print("step ",step," on",len(Bonds))
        SnapshotAngles=[]
        SnapshotBonds=Bonds[step]
        MySnapshot=AllSnapshots[step]
        fa.write('\nstep : '+str(step*Nsteps)+'\ntime :'+str(step*TimeStep*Nsteps)+" fs \n")
        line='Angles : \n'
        Nangles=0
        for internAt in SnapshotBonds :
            if internAt<=CentMax and internAt>=CentMin:
                for iexternAt in SnapshotBonds[internAt]:
                    if iexternAt<=OutMax and iexternAt>=OutMin:
                        for jexternAt in SnapshotBonds[internAt]:
                            if jexternAt<=OutMax and jexternAt>=OutMin:
                                if iexternAt<jexternAt :    
                                    Nangles+=1
                        
                                    xEi,yEi,zEi=MySnapshot[iexternAt][0],MySnapshot[iexternAt][1],MySnapshot[iexternAt][2]
                                    xEj,yEj,zEj=MySnapshot[jexternAt][0],MySnapshot[jexternAt][1],MySnapshot[jexternAt][2]
                                    xI,yI,zI=MySnapshot[internAt][0],MySnapshot[internAt][1],MySnapshot[internAt][2]
                        
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
        
                                    angle=math.acos((distIEi+distIEj-distEiEj)/(2*math.sqrt(distIEi*distIEj)))/math.pi*180
                                    
                                    line+=str(angle)+"\t("+str(iexternAt)+"-"+str(internAt)+"-"+str(jexternAt)+") \n"
                        
                                    SnapshotAngles.append(angle)
        fa.write(line)

        TotNangles+=Nangles
        if Nangles!=0:
            mean=sum(SnapshotAngles)/Nangles                  
            msd=0
            TotMean+=sum(SnapshotAngles)
            for angle in SnapshotAngles:
                msd+=(angle-mean)**2
            msd/=mean
            msd=math.sqrt(msd)
        
            fa.write("mean : "+str(mean)+"\tmsd : "+str(msd)+"\n")
        fa.write("end of step\n")

    if(TotNangles!=0):
        TotMean/=TotNangles

    fa.write("\n\n Mean of all angles : "+str(TotMean))
    fa.close()        

def main(argv):
    umdpf.headerumd()
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
            header = header + 'FILE: -f=' + UMDname+'\t'
        elif opt in ("-b","--bBondFile"):
            BondFile = str(arg)
            header = header + ' -b=' + BondFile+'\t'
        elif opt in ("-s","--sNsteps"):
            Nsteps = int(arg)
            header = header + ' -s=' + arg
            print('Will sample the data every ',Nsteps,' steps')        
        elif opt in ("-c","--cCentralAtom"):
            cen=str(arg)
        elif opt in ("-e","--eExternAtom"):
            ext=str(arg)
            header=header+"\t angles : "+ext+"-"+cen+"-"+ext + "  (degrees)"

    MyCrystal = cr.Lattice()
    AllSnapshots = [cr.Lattice]
        
    CentMin,CentMax,AdjMin,AdjMax,MyCrystal,Bonds,TimeStep = umdpf.read_bonds(BondFile,cen,ext,Nsteps,"Dictionary")

    (MyCrystal,AllSnapshots,TimeStep,length)=umdpf.read_values(UMDname,"xcart","lists",Nsteps)


    File=str(UMDname)+".angles_umdpf.dat"
    
    fa=open(File,'w')
    fa.write(header)
    compute_angles(Bonds,MyCrystal,AllSnapshots,CentMin,CentMax,AdjMin,AdjMax,fa,TimeStep,Nsteps)    
    
    print("Angles computed successfully. File created under the name ",File)

if __name__=='__main__':
   main(sys.argv[1:])
