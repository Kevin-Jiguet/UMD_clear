#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 11:00:43 2023

@author: Kevin Jiguet-Covex
"""
import umd_process as umdp
import sys, getopt


def main(argv):
    umdp.headerumd()
    UMDname = ''
    try:
        opts, arg = getopt.getopt(argv,"hf:",['fUMDfile'])
    except getopt.GetoptError:
        print ('UMD_to_LAMMPS.py -f <umdfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('UMD_to_LAMMPS.py to convert a .umd file into a LAMMPS-like file')
            sys.exit()
        elif opt in ("-f",'fUMDfile'):
            UMDname = str(arg)
        
    fa = open(UMDname[:-8]+'.lammps','w')#creating the lammps file
    ff = open(UMDname,'r')
    #We read the header to extract some information about the system, such as the number of atoms of each type
    line = ff.readline().strip().split()
    natom = int(line[1])
    ff.readline()
    ff.readline()
    ff.readline()
    line = ff.readline().strip().split()
    typat = [line[i] for i in range(1,natom+1)]
    ff.readline()
    line = ff.readline().strip().split()
    timestep = line[1]
    #We fill the file with the atoms coordinates
    while True :
        line = ff.readline()
        if not line : 
            break
        l = line.strip().split()
        if len(l)>0:
            if l[0]=='time':
                fa.write(str(natom)+'\nAtoms. TimeStep '+ str(int(float(l[1])/float(timestep)))+'\n')
                print('Converting TimeStep ',l[1])
            elif l[0]=='atoms:':
                for atom in range(natom):
                    newstring = typat[atom]+' '
                    line = ff.readline().strip().split()
                    newstring+=line[3]+' '+line[4]+' '+line[5]+'\n'#Only the xcart (cartesian) ones
                    fa.write(newstring)
                    
    
    fa.close()
    ff.close()   
        
    print('LAMMPS file successfully created under the name <'+UMDname[:-8]+'.lammps'+'>')
        

if __name__ == "__main__":
   main(sys.argv[1:])
