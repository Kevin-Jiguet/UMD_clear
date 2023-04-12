#!/usr/bin/env python3

###
##AUTHORS: RAZVAN CARACAS
###

import sys,getopt,os
import crystallography as cr
import umd_process as umdp

def print_xyz(MyCrystal,AllSnapshots,UMDname,firststep,iterstep):
    xyzfile = UMDname[:-7] + 'xyz'
    ff = open(xyzfile,'w')
    for istep in range(firststep,len(AllSnapshots),iterstep):
        string = str(MyCrystal.natom) + '\n' + UMDname +'\n'
        ff.write(string)
        for iatom in range(MyCrystal.natom):
            string = MyCrystal.elements[MyCrystal.typat[iatom]] + ' ' + str(AllSnapshots[istep].atoms[iatom].xcart[0]) + ' ' + str(AllSnapshots[istep].atoms[iatom].xcart[1]) + ' ' + str(AllSnapshots[istep].atoms[iatom].xcart[2]) +'\n'
            ff.write(string)
    ff.close()

def main(argv):
    iterstep = 1
    firststep = 0
    UMDname = 'output.umd.dat'
    UMDname = ''
    umdp.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:i:n:")
    except getopt.GetoptError:
        print ('umd2xyz.py -f <umdfile> -i <InitialStep> -s <Sampling_Frequency>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('umd2xyz.py program to write an xyz file with the atomic trajectories from the umd file')
            print ('umd2xyz.py -f <umdfile> -i <InitialStep> -s <Sampling_Frequency>')
            print (' default values: -f output.umd.dat -i 0 -s 1')
            sys.exit()
        elif opt in ("-f"):
            UMDname = str(arg)
        elif opt in ("-i"):
            firststep = int(arg)
        elif opt in ("-n"):
            iterstep = int(arg)
    if (os.path.isfile(UMDname)):
        print('The first ',firststep, 'timesteps will be discarded')
        print('The XYZ file contains every ',iterstep,' timesteps')
        MyCrystal = cr.Lattice()
        AllSnapshots = [cr.Lattice]
        (MyCrystal,AllSnapshots,TimeStep)=umdp.readumd(UMDname)
        print_xyz(MyCrystal,AllSnapshots,UMDname,firststep,iterstep)
    else:
        print ('the umdfile ',UMDname,' does not exist')
        sys.exit()

if __name__ == "__main__":
   main(sys.argv[1:])
