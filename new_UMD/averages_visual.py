#!/usr/bin/python
import numpy
import sys
import getopt
import os
import subprocess
import matplotlib.pyplot as plt


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def main(argv):
    FileName = ''
    Pattern = ''
    data = []
    SkipSteps=0
    average = 0 
    stdev = 0 
    variance = 0 
    d=0
    anchor=''
    try:
        opts, arg = getopt.getopt(argv,"hf:p:s:",["fFileName","pPattern","sSkipSteps"])
    except getopt.GetoptError:
        print('average.py -f <FileName> -p <Pattern> -s <SkipSteps> ')
        sys.exit(d)
    for opt, arg in opts:
        if opt == '-h':
            print('average.py program to extract and average numerical values')
            print('average.py -f <FileName> -p <Pattern> -s <SkipSteps>')
            sys.exit()
        elif opt in ("-f", "--fFileName"):
            FileName = str(arg)
        elif opt in ("-p", "--pPattern"):
            Pattern = str(arg)
            anchor=Pattern.split()[0]
        elif opt in ("-s", "--sSkipSteps"):
            SkipSteps = int(arg)
    if(anchor==''):
        print('No parameter given ; no mean calculated. Use -p to specify a parameter.')
        sys.exit()
    if (os.path.isfile(FileName)): 
        try :
            patterns=subprocess.check_output(['grep',Pattern,FileName])
        except subprocess.CalledProcessError :
            print('Parameter ',Pattern,' not found.')
            sys.exit()
        patternsstr=patterns.decode()
        greps=patternsstr.split('\n')
        for isteps in range(SkipSteps+1,len(greps)):
            elems=greps[isteps].split()
            for ii in range(len(elems)):
                if elems[ii] == anchor:
                    for jj in range(ii+1,len(elems)):
                        if is_number(elems[jj]):
                            data.append(float(elems[jj]))
                            break
        average = sum(data)/len(data)
        for ii in range(len(data)):
            variance = variance + (data[ii]-average)**2
        variance = variance/len(data)
        stdev = numpy.sqrt(variance)
        plt.plot(data)
        plt.show()
        print('Averages over ',len(data),' ares: mean = ',average,' variance = ', variance, ' stdev = ', stdev)
    else:
        print('No input file or file ',FileName,'does not exist')
        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])

