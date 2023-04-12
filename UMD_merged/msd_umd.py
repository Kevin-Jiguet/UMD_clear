#!/usr/bin/env python3
###
##AUTHORS: RAZVAN CARACAS, ANAIS KOBSCH
###

import sys,getopt,os.path
import crystallography as cr
import umd_process as umdp


def msd(MyCrystal,AllSnapshots,TimeStep,hh,vv,ballistic,umdfile): 
    
    niter = len(AllSnapshots)
    msd = [[0.0 for x in range(niter-ballistic+1)] for x in range(MyCrystal.ntypat+1)]
    weight = [0 for x in range(niter)]
    ref = [0.0 for x in range(3)]
    currZ = 0
    for iatom in range(MyCrystal.natom):
        currZ = MyCrystal.typat[iatom]
        for ii in range(ballistic+hh,int(niter/2),hh):          #initial time t
            for kk in range(3): 
                ref[kk]=AllSnapshots[ii].atoms[iatom].absxcart[kk]     #coorindates of the reference atom at time t
            for jj in range(ballistic,int(niter/2),vv):         #new time tao
                msd[currZ][jj] = msd[currZ][jj] + (AllSnapshots[ii+jj].atoms[iatom].absxcart[0]-ref[0])**2 +(AllSnapshots[ii+jj].atoms[iatom].absxcart[1]-ref[1])**2 +(AllSnapshots[ii+jj].atoms[iatom].absxcart[2]-ref[2])**2       #distance between t+tao - t
                weight[jj]=weight[jj]+1                           #normalization to the number of t times (larger tao's are counted fewer times) mutiplied by the number of atoms .
                                                                #when applyging the final normalization this number has to be divided by the total number of atoms to give only the number of t's:
                                                                # weight = sum_no.atoms sum_t's (0 ... tao)
                                                                # normalization is no. of taos
                                                                # which is this weight divided by no.atoms
    msdfile = umdfile[:-8] + '.msd.dat'
    f = open(msdfile,'w')
    string='time_(fs)\t'
    for itypat in range(MyCrystal.ntypat):
        string=string + MyCrystal.elements[itypat] + '\t'
    string = string + '\n'
    f.write(string)
    for ii in range(int(niter/2)):                              # all possible times, max is total_simulation_time / 2
        if (weight[ii]>0):                                      # further check to be sure that tao's have been counted
            string = str(float(ii)*TimeStep)
            for jj in range(MyCrystal.ntypat):
                string = string + '\t' + str(msd[jj][ii]/(float(MyCrystal.types[jj])*float(weight[ii])/float(MyCrystal.natom)))
            string = string + '\n'
            f.write(string)
    print ('MSDs printed in file ',msdfile)

def main(argv):
    hh = 1
    vv = 1
    ballistic = 0
    TimeStep = 1
    umdp.headerumd()
    try:
        opts, arg = getopt.getopt(argv,"hf:z:v:b:",["fumdfile","zHorizontalJump","vVerticalJump","bBallistic"])
    except getopt.GetoptError:
        print ('msd_umd.py -f <XYZ_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('Program to compute the Mean Square Displacement. Usage: ')
            print ('msd_umd.py -f <UMD_filename> -z <HorizontalJump> -v <VerticalJump> -b <Ballistic>')
            print ('UMD_filename = input file with the trajectory, in UMD format.')
            print (' As the MSD is measured with a sliding window of various size up to half the trajectory\'s length, ')
            print (' we can accelerate the calculation using a selected reduced sampling ')
            print ('HorizontalJump = discretization for the start of the sampling window.')
            print ('VerticalJump = discretization for the length of the sampling window.')
            print ('Ballistic = estimation of the ballistic part of the trajectory. Default is 0. Typical values of 100 are sufficient.')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            print ('I will use the ',umdfile,' for input','\n')
        elif opt in ("-z", "--zHorizontalJump"):
            hh = int(arg)
            if hh<=0 :
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
        msd(MyCrystal,AllSnapshots,TimeStep,hh,vv,ballistic,umdfile)
    else:
        print ('umd file ',umdfile,'does not exist')
        sys.exit()
 

if __name__ == "__main__":
   main(sys.argv[1:])
