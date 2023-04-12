#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS, ANAIS KOBSCH
###

import sys,getopt,os.path
import crystallography as cr
import umd_process as umdp


def msd(MyCrystal,AllSnapshots,TimeStep,hh,vv,ballistic,umdfile):
    
    niter = len(AllSnapshots)
    print('Nombre de pas = ', niter-ballistic+1)
    msd = [[0.0 for x in range(niter-ballistic+1)] for x in range(MyCrystal.ntypat + 1 + MyCrystal.natom)]
    weight = [0 for x in range(niter)]
    ref = [0.0 for x in range(3)]
    currZ = 0
    for iatom in range(MyCrystal.natom):                        #compute the MSDs for each individual atom
        currZ = MyCrystal.typat[iatom]
        for ii in range(ballistic+hh,int(niter/2),hh):          #initial time t
            for kk in range(3): 
                ref[kk]=AllSnapshots[ii].atoms[iatom].absxcart[kk]     #coordinates of the reference atom at time t
            for jj in range(ballistic,int(niter/2),vv):         #new time tao
                msd[iatom][jj] = msd[iatom][jj] + (AllSnapshots[ii+jj].atoms[iatom].absxcart[0]-ref[0])**2 +(AllSnapshots[ii+jj].atoms[iatom].absxcart[1]-ref[1])**2 +(AllSnapshots[ii+jj].atoms[iatom].absxcart[2]-ref[2])**2       #distance between t+tao - t
                weight[jj]=weight[jj]+1.0                         #normalization to the number of t times (larger tao's are counted fewer times)

    for jj in range(ballistic,int(niter/2),vv):
        weight[jj]=weight[jj]/float(MyCrystal.natom)
                                                                #this wieght carries also the number of atoms with it


    for iatom in range(MyCrystal.natom):                        #combine the individual MSDs for the different atoms per atom type
        currZ = MyCrystal.typat[iatom]
        for jj in range(ballistic,int(niter/2),vv):
            msd[MyCrystal.natom+currZ][jj] = msd[MyCrystal.natom+currZ][jj] + msd[iatom][jj]

    msdfile = umdfile[:-8] + '.msd_all.dat'
    f = open(msdfile,'w')
    string='time_(fs)\t'
    for iatom in range(MyCrystal.natom):                        #prints the header with all the checmial elements, indexed by the current number of the atom
        string = string + MyCrystal.elements[MyCrystal.typat[iatom]] + str(iatom) + '\t'
    for itypat in range(MyCrystal.ntypat):
        string=string + MyCrystal.elements[itypat] + '\t'
    string = string + '\n'
    f.write(string)
    for ii in range(int(niter/2)):                              # all possible times, max is total_simulation_time / 2
        if (weight[ii]>0):                                      # further check to be sure that tao's have been counted
            string = str(float(ii)*TimeStep)
            for iatom in range(MyCrystal.natom):                #first writes the individual atomic MSDs
                string = string + '\t' + str(msd[iatom][ii]/float(weight[ii]))       #double division with natom as it is counted in weight
            for jj in range(MyCrystal.ntypat):                  #now, write the cummulate MSDs by atomic type
                string = string + '\t' + str(msd[MyCrystal.natom+jj][ii]/(float(MyCrystal.types[jj])*float(weight[ii])))
            string = string + '\n'
            f.write(string)
    print ('MSDs printed in file ',msdfile)

def main(argv):
    zz = 1
    vv = 1
    ballistic = 0
    TimeStep = 1
    umdp.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:z:v:b:",["fumdfile","zHorizontalJump","vVerticalJump","bBallistic"])
    except getopt.GetoptError:
        print ('msd_sparse.py -f <XYZ_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Program to compute the mean square displacements first of all the individual atoms and and then averaged by atomic types. Usage: ')
            print ('msd_sparse.py -f <UMD_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic>')
            print ('UMD_filename = input file with the trajectory, in UMD format.')
            print (' As the MSD is measured with a sliding window of various size up to half the trajectory\'s length, ')
            print (' we can accelerate the calculation using a selected reduced sampling ')
            print ('HorizontalJump = discretization for the start of the sampling window.')
            print ('VerticalJump = discretization for the length of the sampling window.')
            print ('Ballistic = estimation of the ballistic part of the trajectory. Default is 0. Typical values of 100 are sufficient.')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            print ('I will use the ',umdfile,' for input')
        elif opt in ("-z", "--zHorizontalJump"):
            zz = int(arg)
            if zz<=0 :
                print('ERROR : horizontal jump -z has to be strictly positive.')
                sys.exit()
        elif opt in ("-v", "--vVerticalJump"):
            vv = int(arg)
            if vv<=0 :
                print('ERROR : vertical jump -v has to be strictly positive.')                
                sys.exit()
        elif opt in ("-b", "--bBallistic"):
            ballistic = int(arg)
    if (os.path.isfile(umdfile)):
        MyCrystal = cr.Lattice()
        AllSnapshots = [cr.Lattice]
        (MyCrystal,AllSnapshots,TimeStep)=umdp.read_absxcart(umdfile)
        print ('Number of atoms of each type is ',MyCrystal.types)
        msd(MyCrystal,AllSnapshots,TimeStep,zz,vv,ballistic,umdfile)
    else:
        print ('umd file ',umdfile,'does not exist')
        sys.exit()

 

if __name__ == "__main__":
   main(sys.argv[1:])

