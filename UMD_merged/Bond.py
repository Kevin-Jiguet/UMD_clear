import sys,getopt,os.path
import crystallography as cr
import umd_process as umdp

#This script calculates which atom is bound with each one. It creates a .bonding.dat file containing all the relevant information. This file is then used by speciation_alt.py.

def read_inputfile(InputFile,MyCrystal):
    BondTable = [[0.0 for _ in range(MyCrystal.ntypat)] for _ in range(MyCrystal.ntypat)]
    nolines = 0
    with open(InputFile) as ff:
        for line in ff:
            line=line.strip()
            entry=line.split()
            if (len(entry)==3):
                nolines +=1
                for ii in range(MyCrystal.ntypat):
                    if MyCrystal.elements[ii]==entry[0]:
                        for jj in range(MyCrystal.ntypat):
                            if MyCrystal.elements[jj]==entry[1]:
                                BondTable[ii][jj]=float(entry[2])*float(entry[2])
                                BondTable[jj][ii]=float(entry[2])*float(entry[2])
    return BondTable


def WriteBonding2(MyCrystal,Mysnapshot,BondTable,step,timestep,File,natoms):
    header = 'time '+str(step*timestep)+' fs\nstep '+str(step)+'\n'
    File.write(header)
#    print(MyCrystal.typat[1])
    for iatom in range(natoms):
        [x1,y1,z1]=Mysnapshot.atoms[iatom].xcart
        typeatom1=MyCrystal.typat[iatom]
        line=str(iatom)+'\t'
        for jatom in range(natoms):
            if (jatom!=iatom):
                [x2,y2,z2]=Mysnapshot.atoms[jatom].xcart
                typeatom2=MyCrystal.typat[jatom]
                dx,dy,dz=(x1-x2),(y1-y2),(z1-z2)
                
                valx = min(dx**2, (MyCrystal.acell[0]-dx)**2, (MyCrystal.acell[0]+dx)**2)
                valy = min(dy**2, (MyCrystal.acell[1]-dy)**2, (MyCrystal.acell[1]+dy)**2)
                valz = min(dz**2, (MyCrystal.acell[2]-dz)**2, (MyCrystal.acell[2]+dz)**2)
                
                distsquared=valx+valy+valz
                
                if(distsquared < BondTable[typeatom1][typeatom2]):
                    line+=str(jatom)+'\t'
        line+='\n'
        File.write(line)
    File.write('end\n')
    
def WriteBonding(MyCrystal,MySnapshot,BondTable,step,timestep,File,natoms):
#    print(MyCrystal.typat[1])
    lines=[[at] for at in range(natoms)]
    for iatom in range(natoms):
        for jatom in range(iatom+1,natoms):
                dx = MySnapshot.atoms[jatom].xcart[0] - MySnapshot.atoms[iatom].xcart[0]
                dy = MySnapshot.atoms[jatom].xcart[1] - MySnapshot.atoms[iatom].xcart[1]
                dz = MySnapshot.atoms[jatom].xcart[2] - MySnapshot.atoms[iatom].xcart[2]
                
                valx = min(dx**2, (MyCrystal.acell[0]-dx)**2, (MyCrystal.acell[0]+dx)**2)
                valy = min(dy**2, (MyCrystal.acell[1]-dy)**2, (MyCrystal.acell[1]+dy)**2)
                valz = min(dz**2, (MyCrystal.acell[2]-dz)**2, (MyCrystal.acell[2]+dz)**2)
                distij = valx + valy + valz
                    
                if (distij < BondTable[MyCrystal.typat[iatom]][MyCrystal.typat[jatom]]): #if the distance between the two atoms is below the cutoff distance between the two atoms
                        lines[iatom].append(jatom)
                        lines[jatom].append(iatom)
        
    header = 'time '+str(step*timestep)+' fs\nstep '+str(step)+'\n'
    File.write(header)
    for line in lines :
        l=''
        for atom in line :
                l+=str(atom)+'\t'
        l+='\n'
        File.write(l)
    File.write('end\n')

def main(argv):
    umdp.headerumd()
    UMDname='output.umd.dat'
    Nsteps = 1
    InputFile = ''
    header = ''
    try:
        opts, arg = getopt.getopt(argv,"hf:s:i:",["fUMDfile","sSampling_Frequency","iInputFile"])
    except getopt.GetoptError:
        print ('Bond.py -f <UMD_filename> -s <Sampling_Frequency> -i <InputFile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('computation of bonding map')
            print ('bond.py -f <UMD_filename> -s <Sampling_Frequency>  -i <InputFile> ')
            print ('default values: -f output.umd.dat -s 1')
            print ('the UMD file contains the position and nature of each atom at each time.'
            print ('the input file contains the bond lengths for the different atom pairs.')
            sys.exit()
        elif opt in ("-f", "--fUMDfile"):
            UMDname = str(arg)
            header = header + 'FILE: -f=' + UMDname
        elif opt in ("-s","--sNsteps"):
            Nsteps = int(arg)
            header = header + ' -s=' + arg
            print('Will sample the MD trajectory every ',Nsteps,' steps')
        elif opt in ("-i", "--iInputFile"):
            InputFile = str(arg)
            header = header + ' -i=' + InputFile
            print ('Bonding cutoffs to be read from file ',InputFile)
        

    if not (os.path.isfile(UMDname)):
        print ('the UMD files ',UMDname,' does not exist')            
        sys.exit()


    MyCrystal = cr.Lattice()
    AllSnapshots = [cr.Lattice]
    (MyCrystal,AllSnapshots,TimeStep)=umdp.read_xcart(UMDname)       

    if len(InputFile)>0:
        BondTable = read_inputfile(InputFile,MyCrystal)

 
    natom=MyCrystal.natom
    FileAll=UMDname+'.bonding.dat'
    print ('Bondings will be written in ',FileAll,' file')
    fa=open(FileAll,'w')
    ff=open(UMDname,'r')
    
    for i in range(20):    
        line=ff.readline()
        fa.write(line)
        
    ff.close()
    
    
    for step in range(0,len(AllSnapshots),Nsteps):
        print('step',step,'on',len(AllSnapshots))
        WriteBonding(MyCrystal,AllSnapshots[step],BondTable,step,TimeStep,fa,natom)            
    fa.close()

if __name__ == "__main__":
    main(sys.argv[1:])
