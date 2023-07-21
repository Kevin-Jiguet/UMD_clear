#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 07:43:57 2023

@author: Kevin Jiguet-Covex
"""
import sys, getopt, os.path, itertools, math
import numpy as np

def main(argv):
    
    filename = ""
    logname = ""
    allname = ""
    press = False
    
    try:
        opts,args = getopt.getopt(argv,"f:l:a:",["fFile","lLog","aAllpress"])
    except getopt.GetoptError:
        print("-f : LAMMPS file ; -l : log file ; -a : allpress file")
        sys.exit()
    for opt,arg in opts : 
        if opt in "-f":
            filename = str(arg)
        elif opt in "-l":
            logname = str(arg)
        elif opt in "-a":
            allname = str(arg)
            press = True
    print("Will use the file <"+filename+"> for positions") 
    print("Will use the file <"+logname+"> for thermodynamic data")
    print("Will use the file <"+allname+"> for stress tensor values\n") 
        
    fl = open(logname,'r')
    fa = open(filename+".umd.dat",'w')    
    
    natom = 0
    types =  []
    elements = []
    ntypat = 0
    typat = []
    t0 = 0
    style=""
    
    
    while True : 
        line=fl.readline()
        if not line : 
            break
        line = line.strip()
        
        if line == "Reading data file ...":
            line=fl.readline().strip().split()
            acellx = float(line[7].strip('('))-float(line[3].strip('('))
            acelly = float(line[8].strip('('))-float(line[4].strip('('))
            acellz = float(line[9].strip(')'))-float(line[5].strip(')'))
        elif line == "# SIMULATION":
            break
        else :
            line = line.split()
            if line[0] == "timestep":
                line = fl.readline().strip().split()
                timestepLog = float(line[1])*1000
            elif line[0] == "thermo_style":
                style = line[1]
                

    while True :
        line=fl.readline()
        if not line : 
            break
        line = line.strip().split()
        if len(line)>0 and line[0]=="group":
            elements.append(line[1].split("_")[0])
            line = fl.readline().strip().split()
            types.append(int(line[0]))
            natom+=types[-1]
            typat+=[ntypat for _ in range(types[-1])]
            ntypat+=1
        if len(line)>2 and line[0] == "Per" and line[1] == "MPI" and line[2] == "rank":
            if style == "one":
                Variables = fl.readline().strip().split()
                IEindex = None
                Enthindex = None
                Tempindex = None
                Pressindex = None
                for i in range(len(Variables)) : 
                    v=Variables[i]
                    if v == "Temp":
                        Tempindex = i
                    elif v == "Press":
                        Pressindex = i
                    elif v == "TotEng":
                        IEindex = i
                    elif v =="Enthalpy":
                        Enthindex = i
            break

    
    ff=open(filename,'r')
    t0=None
    while True :
        l = ff.readline()
        if not l :
            break
        line = l.strip().split()
        if len(line)>1 and line[1]=="Timestep:":
            if t0 == None :
                t0=float(line[2])
            else : 
                timestepFile = (float(line[2]) - t0)*timestepLog
                break
    ff.close()

    fa.write("natom "+str(natom)+"\n")
    fa.write("ntypat "+str(ntypat)+"\n")
    string = ""
    for n in types :
        string+=str(n)+" "
    string+="\n"
    fa.write("types "+string)
    string=""
    for el in elements :
        string+=el+" "
    string+="\n"
    fa.write("elements "+string)
    string=""
    for x in typat :
        string+=" "+str(x)
    string+="\n\n"
    fa.write("typat "+string)




    IE = '0'
    Enth = '0'
    Press = '0'
    Temp = '0'

    ff=open(filename,'r')
    if press :
        fall=open(allname,'r')
    
    while True :
        line = ff.readline()
        if not line :
            break
        line = line.strip().split()
        if line[0] == "Atoms.":
            string="timestep "+str(timestepFile)+" fs\n"
            fa.write(string)
            string = "time "+str(float(line[2])*timestepLog)+" fs\n"
            fa.write(string)
            logline = fl.readline().strip().split()
            
            if style == "one":         
                while logline[0] != line[2]:
                    logline = fl.readline().strip().split()
                    if not logline : 
                        break
                if Tempindex != None :
                    Temp = logline[Tempindex]
                if Pressindex != None :
                    Press = str(float(logline[Pressindex])/10000)
                if Enthindex != None :
                    Enth = logline[Enthindex]
                if IEindex != None :
                    IE = logline[IEindex]

            elif style == "multi":
                while logline[2] != line[2]:
                    logline = fl.readline().strip().split()
                    if not logline : 
                        break
                logline = fl.readline().strip().split()
                IE = logline[2]
                Temp = logline[8]
                logline = fl.readline().strip().split()
                logline = fl.readline().strip().split()
                logline = fl.readline().strip().split()
                Press = str(float(logline[8])/10000)
            
            rprimVecs = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
            cellVecs = [np.array([1.0,0.0,0.0])*acellx,np.array([0.0,1.0,0.0])*acellx,np.array([0.0,0.0,1.0])*acellx]
                
            rprimVecsstr=["  ".join(str(el) for el in vec) for vec in rprimVecs]
            cellVecsstr=["  ".join(str(el) for el in vec) for vec in cellVecs]
                
            acellLine = "acell "+str(acellx)+" "+str(acelly)+" "+str(acellz)+" A\n"
            rLines ="rprim_a "+rprimVecsstr[0] +" \nrprim_b "+ rprimVecsstr[1] + " \nrprim_c " +rprimVecsstr[2]+ " \n"
            rdLines ="rprimd_a "+cellVecsstr[0] +" A\nrprimd_b "+ cellVecsstr[1] + " A\nrprimd_c " +cellVecsstr[2]+ " A\n"
            
            string = "Internal Energy "+IE+" eV\nEnthalpy "+Enth+" eV\nTemperature "+Temp+" K\nPressure "+Press+" GPa\n"
            
            if press :
                allLine = fall.readline().strip().split()
            
                while allLine[0]!=line[2]:
                    allLine = fall.readline().strip().split()
                    if not allLine :
                        break
            
                stressline = "StressTensor "+str(float(allLine[1])/10000)+" "+str(float(allLine[2])/10000)+" "+str(float(allLine[3])/10000)+" "+str(float(allLine[4])/10000)+" "+str(float(allLine[5])/10000)+" "+str(float(allLine[6])/10000)+" GPa\n"

                    
            
        
            fa.write(string)
            if press :
                fa.write(stressline)
            fa.write(acellLine)
            fa.write(rLines)
            fa.write(rdLines)            
            
            
            string = "atoms: reduced*3 cartesian*3(A) abs.diff.*3(A)\n"
            fa.write(string)
            for at in range(natom):
                line = ff.readline()
                line = line.strip().split()
                x,y,z = float(line[1]),float(line[2]),float(line[3])
                redx,redy,redz = x/acellx,y/acelly,z/acellz
                string = str(redx)+" "+str(redy)+" "+str(redz)+" "+str(math.fmod(acellx+x,acellx))+" "+str(math.fmod(acelly+y,acelly))+" "+str(math.fmod(acellz+z,acellz))+" "+str(x)+" "+str(y)+" "+str(z)
                fa.write(string+"\n")
            
            fa.write("\n")
    
    ff.close()
    fl.close()
    if press :
        fall.close()
    fa.close()
    print("umd file created under the name <"+filename+".umd.dat>")
            
        

if __name__ =="__main__":
    main(sys.argv[1:])