#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 14:54:47 2023

@author: Kevin Jiguet-Covex
"""

import sys,getopt,os.path
import crystallography as cr
import umd_processes_fast as umdpf
import time
import numpy as np
import concurrent.futures
from functools import partial
import ctypes
from os.path import join
from scipy.fftpack import dct, fftfreq
import math
import matplotlib.pyplot as plt
import platform

current_path=os.path.abspath(__file__)#For all this to work, the file c_autocorrelation_vib.so must be in the same directory than this script
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u

OS = platform.system()
LibraryName =""

if OS == "Linux":
    LibraryName = 'c_vibr_clusters.so'
elif OS == "Windows":
    LibraryName = 'c_vibr_clusters.dll'
elif OS == "Darwin":
    LibraryName = 'c_vibr_clusters.dylib'


clust_lib = ctypes.cdll.LoadLibrary(join(path_new, LibraryName))

clust_lib.correlation.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                  ctypes.POINTER(ctypes.c_double),
                                  ctypes.POINTER(ctypes.c_double), 
                                  ctypes.c_int, 
                                  ctypes.c_int, 
                                  np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C'), 
                                  np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C'), 
                                  np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C'), 
                                  np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C'),
                                  np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C'),
                                  np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C'),
                                  ctypes.c_int,
                                  ctypes.c_double]
clust_lib.correlation.restype = None
clust_lib.free_memory.argtypes= [ctypes.POINTER(ctypes.c_double)]
clust_lib.free_memory.restype = None


def nan_in(list):
    for i in list :
        if math.isnan(i):
            return True
    return False

def correlation_C(PosMatrix,VelMatrix,timestep,numsteps,masses,natom):
    maxtau = numsteps//2
    Vel_autoc_X = np.zeros(maxtau*natom,dtype=np.float64)
    Vel_autoc_Y = np.zeros(maxtau*natom,dtype=np.float64)
    Vel_autoc_Z = np.zeros(maxtau*natom,dtype=np.float64)
    Rad_autoc = np.zeros(maxtau*(natom-1),dtype=np.float64)
    Theta_autoc = np.zeros(maxtau*(natom-1),dtype=np.float64)
    Phi_autoc = np.zeros(maxtau*(natom-1),dtype=np.float64)

    massesP = (ctypes.c_double * len(masses))(*masses)
    PosMatrixP = (ctypes.c_double * len(PosMatrix))(*PosMatrix)
    VelMatrixP = (ctypes.c_double * len(VelMatrix))(*VelMatrix)
    
    clust_lib.correlation(PosMatrixP,VelMatrixP,massesP,natom,numsteps,Vel_autoc_X,Vel_autoc_Y,Vel_autoc_Z,Rad_autoc,Theta_autoc,Phi_autoc,len(PosMatrix),timestep)

    fft_correlation_VelX = []
    fft_correlation_VelY = []
    fft_correlation_VelZ = []
    fft_correlation_Rad = []
    fft_correlation_Theta = []
    fft_correlation_Phi = []
    
    fft_correlation_VelX.append(dct(Vel_autoc_X[:maxtau],1)/2.0 * 2.0 * timestep)
    fft_correlation_VelY.append(dct(Vel_autoc_Y[:maxtau],1)/2.0 * 2.0 * timestep)
    fft_correlation_VelZ.append(dct(Vel_autoc_Z[:maxtau],1)/2.0 * 2.0 * timestep)

    
    for at in range(natom-1):
        fft_correlation_VelX.append(dct(Vel_autoc_X[(at+1)*maxtau:(at+2)*maxtau],1)/2.0 * 2.0 * timestep)
        fft_correlation_VelY.append(dct(Vel_autoc_Y[(at+1)*maxtau:(at+2)*maxtau],1)/2.0 * 2.0 * timestep)
        fft_correlation_VelZ.append(dct(Vel_autoc_Z[(at+1)*maxtau:(at+2)*maxtau],1)/2.0 * 2.0 * timestep)
        fft_correlation_Rad.append(dct(Rad_autoc[at*maxtau:(at+1)*maxtau],1)/2.0 * 2.0 * timestep)
        fft_correlation_Theta.append(dct(Theta_autoc[at*maxtau:(at+1)*maxtau],1)/2.0 * 2.0 * timestep)
        fft_correlation_Phi.append(dct(Phi_autoc[at*maxtau:(at+1)*maxtau],1)/2.0 * 2.0 * timestep) 

    freq = fftfreq(2*int(numsteps/2),d=timestep)[:int(numsteps/2)]
    return [Vel_autoc_X,Vel_autoc_Y,Vel_autoc_Z,Rad_autoc,Theta_autoc,Phi_autoc],[fft_correlation_VelX,fft_correlation_VelY,fft_correlation_VelZ,fft_correlation_Rad,fft_correlation_Theta,fft_correlation_Phi,freq]


def correlation(TimeMatrix,timestep):
  
    # TimeMatrix should be in matrix format
    # entry1 entry2 entry3 ..
    # 0 1 2 3 4 5 6 
    # 1 1 2 3 4 5 6
    #.....
#    print("TimeMatrix=",TimeMatrix)
    
    nostep = len(TimeMatrix)
    noentries = len(TimeMatrix[1])
    maxtau = int(nostep / 2)
    autocorrelation = np.empty((maxtau,noentries))
    fft_correlation = np.empty((maxtau,noentries)) 
    temp = 1.0/np.arange(nostep,nostep-maxtau,-1)
    normalization = np.diag(temp)        
    
    for ientry in range(noentries): 

        temp1 = np.correlate(TimeMatrix[:,ientry],TimeMatrix[:,ientry],mode='full')[len(TimeMatrix[:,ientry])-1:] #same as that of [len(TimeMatrix[ientry])-1:]
        #although it's fast, it's not normalized
        autocorrelation[:,ientry] = np.matmul(normalization,temp1[0:maxtau]) 
        fft_correlation[:,ientry] = dct(autocorrelation[:,ientry],1)/2.0 * 2.0 * timestep #this gives exactly the same answer as above and hope you know what is time shifting in discrete fourier transform(that is the reason to have np.abs)

#    print("autoc=",autocorrelation[0])
    # we calculate frequency
    freq = fftfreq(2*maxtau,d=timestep)[:maxtau] 
#    print(fft_correlation)       
    return autocorrelation,fft_correlation,freq                	

def relative(ListPos,ListVel,Cluster,life,MyCrystal):
    
  #  print(len(ListPos),len(ListVel))
    ListPosM = np.array([ListPos[i]*MyCrystal.masses[MyCrystal.typat[Cluster[i]]] for i in range(len(Cluster))])
    ListVelM = np.array([ListVel[i]*MyCrystal.masses[MyCrystal.typat[Cluster[i]]] for i in range(len(Cluster))])
 #   print(np.shape(ListPosM),np.shape(ListVelM))
    
    totMass=sum([MyCrystal.masses[MyCrystal.typat[atom]] for atom in Cluster])
    
    Pos_center=np.sum(ListPosM,0)/totMass
    Vel_center=np.sum(ListVelM,0)/totMass
    
#    print(np.shape(Pos_center))
    
    relativeP = ListPos-np.array([Pos_center for _ in range(len(Cluster))])
    relativeV = ListVel-np.array([Vel_center for _ in range(len(Cluster))])    
    

#    print(Cluster)
#    print(relativeP[0])
#    print(relativeV[0])
#    sys.exit()
    
    return relativeP,relativeV
    
def main(argv):
    start=time.time()
    maxsize=10
    minsize=0
    umdfile=''
    popfile=''
    minlife=0
    centralatom=''
    outeratom=''
    temperature=5000
    try:
        opts, args = getopt.getopt(argv,"hf:p:c:o:t:s:S:T:",["fumdfile","pPopulation","cCentralatom","oOuteratom","tMintime","sMinsize","SMaxsize","TTemperature"])
    except getopt.GetoptError:
        print ('msd_umd.py -f <umdfilename> -p <populfilename> -c <centralatom> -o <outeratom> -t <mintime> -s <minsize> -S <maxsize> -T <temperature>')
        sys.exit(2)
    print("opts=",opts)

    for o in opts:
        opt=o[0]
        arg=o[1]
        if opt in ('-h', "--help"):
            print ('vibr_clusters_fast_umd.py program to compute the atomic velocity self-correlation of atomic clusters')
            print ('and extract relevant properties: vibrational spectrum and self-correlation')
            print ('vibr_clusters_fast_umd.py -f <umdfilename> -p <populfilename> -c <centralatom> -o <outeratom> -t <mintime> -s <minsize> -S <maxsize> -T <temperature>')
            print ('the program needs an UMD file and a population file')
            sys.exit()
        elif opt in ("-f", "--fumdfile"):
            umdfile = str(arg)
            print ('I will use the ',umdfile,' for input','\n')
        elif opt in ("-p", "--pPopulation"):
            popfile = str(arg)
            print ('I will use the ',popfile,' to determine the population','\n')
        elif opt in ("-c", "--cCentralatom"):
            centralatom = str(arg)
        elif opt in ("-o", "--oOuteratom"):
            outeratom = str(arg)
        elif opt in ("-t", "--tMintime"):
            minlife=float(arg)
        elif opt in ("-S", "--SMaxsize"):
            maxsize=float(arg)
        elif opt in ("-s", "--sMinsize"):
            minsize=float(arg)
        elif opt in ("-t","--tTemperature"):
            temperature = float(arg)

    Boltzmann=8.6173303*10**(-5) #eV
    Avogadro=6.022140857*(10**(23)) # mol-1
    au_angstrom_squre = 1.0 / Avogadro * (10**(-3)) * (10**10) / (1.602176634 * 10**-19) # unit is eV
    kB_T = temperature * Boltzmann
                
    MyCrystal,Values,TimeStep,length = umdpf.read_values(umdfile, "everything")
    
    Centralatoms = []
    Outeratoms = []
    AllAtoms=[]
    ff=open(umdfile,'r')    
    
    for iatom in range(MyCrystal.natom):
        if MyCrystal.elements[MyCrystal.typat[iatom]]==centralatom:
            Centralatoms.append(iatom)    
            AllAtoms.append(iatom)#contains the list with the index of the central atoms from the 0 ... natom
        if MyCrystal.elements[MyCrystal.typat[iatom]]==outeratom:
            Outeratoms.append(iatom) 
            AllAtoms.append(iatom)#contains the list with the index of the coordinating atoms from the 0 ... natom
        
    print("searching for clusters containing only the following atoms : ", AllAtoms)
    
    listVel=[[] for _ in range(MyCrystal.natom)]
    listPos=[[] for _ in range(MyCrystal.natom)]

    
    for step in range(length):
        for at in AllAtoms : 
            listPos[at]+=[Values[step][12*at+6],Values[step][12*at+7],Values[step][12*at+8]]
            listVel[at]+=[Values[step][12*at+9],Values[step][12*at+10],Values[step][12*at+11]]

    dicoClusters={}
    Names=[]
    dicoNames = {}
        
    fp = open(popfile,'r')   
    fp.readline()
    fp.readline()
    while True : 
        line=fp.readline()
        if not line : 
            break
        if len(line)>1 :
            entry=line.strip().split()
            if float(entry[3])>=minlife and minsize<=len(entry[4:])<=maxsize:
                Cluster="["+(line.split('[')[-1]).strip("\n")
                speciesName = ''
                clusterType = [0 for _ in range(MyCrystal.ntypat)]
                Atoms = eval(Cluster)
                for at in Atoms : 
                    clusterType[MyCrystal.typat[at]]+=1
                for indice in range(len(clusterType)) :
                    if clusterType[indice] > 0 :
                        speciesName += MyCrystal.elements[indice]+'_'+str(clusterType[indice])+" "
                        
                ti = int(entry[1])
                tf = int(entry[2])
                 
                if speciesName in dicoNames.keys():
                    dicoNames[speciesName]['Times'].append([ti,tf])
                    dicoNames[speciesName]['Clusters'].append(Atoms)
                else :
                    masses=[]
                    for atom in Atoms:
                        clusterType[MyCrystal.typat[atom]]+=1                    
                        masses.append(MyCrystal.masses[MyCrystal.typat[atom]])
                    dicoNames[speciesName]={'Times':[[ti,tf]],'Clusters':[Atoms],'Masses':masses, 'Correlations':[], 'Spectra':[], "Average Correlations":[],"Average Spectra":[]}
    
    for Name in dicoNames :

        masses = dicoNames[Name]['Masses']
        Clusters = dicoNames[Name]['Clusters']
        Times = dicoNames[Name]['Times']
        maxtime = 0
        maxindex = 0
        numatoms = len(Clusters[0])

        for index in range(len(Times)):
            Atoms = Clusters[index]
            [ti,tf] = Times[index]
            if(tf-ti+1)>maxtime :
                maxtime = tf-ti+1
                maxindex=index
            PosList = []
            VelList = []
            i=0;
            for atom in Atoms : 
                i+=1;
                PosList+=listPos[atom][3*ti:3*(tf+1)]
                VelList+=listVel[atom][3*ti:3*(tf+1)]
                
            Correlations,Spectra=correlation_C(PosList,VelList,TimeStep,tf-ti+1,masses,numatoms)

            dicoNames[Name]['Correlations'].append(Correlations)
            dicoNames[Name]['Spectra'].append(Spectra)

        frequencies = dicoNames[Name]['Spectra'][maxindex][-1]
        MaxTau = int((dicoNames[Name]['Times'][maxindex][1]-dicoNames[Name]['Times'][maxindex][0]+1)/2)

        num_freqbins = len(frequencies)

        NormFreq = np.zeros(num_freqbins)
        NormCorr = np.zeros(MaxTau)

        Avg_fft_Vel_central = np.zeros(num_freqbins)
        Avg_fft_Vel_outer = np.zeros(num_freqbins)
        Avg_fft_Rad = np.zeros(num_freqbins)
        Avg_fft_Theta = np.zeros(num_freqbins)
        Avg_fft_Phi = np.zeros(num_freqbins)
        Avg_corr_Vel_central = np.zeros(MaxTau)
        Avg_corr_Vel_outer = np.zeros(MaxTau)
        Avg_corr_Rad = np.zeros(MaxTau)
        Avg_corr_Theta = np.zeros(MaxTau)
        Avg_corr_Phi = np.zeros(MaxTau)
    
        for t in range(len(dicoNames[Name]['Times'])):
    #        print("t=",t)
            spectrum = dicoNames[Name]['Spectra'][t]
            correlation = dicoNames[Name]['Correlations'][t]
            N=len(spectrum[0][0])                
            maxtau = int(len(correlation[0])/numatoms)

            for i in range(maxtau):

                NormCorr[i]+=1
                Avg_corr_Vel_central[i]+=correlation[0][i]+correlation[1][i]+correlation[2][i]
                    
                for iat in range(1,numatoms):

                    Avg_corr_Vel_outer[i]+=correlation[0][iat*maxtau+i]+correlation[1][iat*maxtau+i]+correlation[2][iat*maxtau+i] 
                    Avg_corr_Rad[i]+=correlation[3][(iat-1)*maxtau+i] 
                    Avg_corr_Theta[i]+=correlation[4][(iat-1)*maxtau+i]
                    Avg_corr_Phi[i]+=correlation[5][(iat-1)*maxtau+i]
                        
            for k in range(N):
                    
                freqbin = int(k/N*num_freqbins + 0.5)
                NormFreq[freqbin]+=1
                Avg_fft_Vel_central[freqbin]+=(spectrum[0][0][k]+spectrum[1][0][k]+spectrum[2][0][k])/3
                
                for iat in range(numatoms-1) :
                
                    Avg_fft_Vel_outer[freqbin]+=(spectrum[0][iat+1][k]+spectrum[1][iat+1][k]+spectrum[2][iat+1][k])/3
                    Avg_fft_Rad[freqbin]+=spectrum[3][iat][k]
                    Avg_fft_Theta[freqbin]+=spectrum[4][iat][k]
                    Avg_fft_Phi[freqbin]+=spectrum[5][iat][k]
            
            
        for k in range(num_freqbins):#Normalization of spectra
            
            Avg_fft_Vel_central[k]*=MyCrystal.masses[MyCrystal.elements.index(centralatom)]* 2 * au_angstrom_squre / kB_T/(NormFreq[k])
            Avg_fft_Vel_outer[k]*=MyCrystal.masses[MyCrystal.elements.index(outeratom)]* 2 * au_angstrom_squre / kB_T/((numatoms-1)*NormFreq[k])
            Avg_fft_Rad[k]*=MyCrystal.masses[MyCrystal.elements.index(outeratom)]* 2 * au_angstrom_squre / kB_T/((numatoms-1)*NormFreq[k])
            Avg_fft_Theta[k]*=MyCrystal.masses[MyCrystal.elements.index(outeratom)]* 2 * au_angstrom_squre / kB_T/((numatoms-1)*NormFreq[k])
            Avg_fft_Phi[k]*=MyCrystal.masses[MyCrystal.elements.index(outeratom)]* 2 * au_angstrom_squre / kB_T/((numatoms-1)*NormFreq[k])
                    
        for i in range(MaxTau):#Normalization of correlations
                
            Avg_corr_Vel_central[i]/=NormCorr[i]/MyCrystal.masses[MyCrystal.elements.index(centralatom)]
            Avg_corr_Vel_outer[i]/=(numatoms-1)*NormCorr[i]/MyCrystal.masses[MyCrystal.elements.index(outeratom)]
            Avg_corr_Rad[i]/=(numatoms-1)*NormCorr[i]/MyCrystal.masses[MyCrystal.elements.index(outeratom)]
            Avg_corr_Theta[i]/=(numatoms-1)*NormCorr[i]/MyCrystal.masses[MyCrystal.elements.index(outeratom)]
            Avg_corr_Phi[i]/=(numatoms-1)*NormCorr[i]/MyCrystal.masses[MyCrystal.elements.index(outeratom)]
            
     #   print(Avg_fft_Vel_central)
        dicoNames[Name]["Average Spectra"]=[Avg_fft_Vel_central,Avg_fft_Vel_outer,Avg_fft_Rad,(Avg_fft_Theta+Avg_fft_Phi)/2,frequencies]
        dicoNames[Name]["Average Correlations"]=[Avg_corr_Vel_central,Avg_corr_Vel_outer,Avg_corr_Rad,(Avg_corr_Theta+Avg_corr_Phi)/2]
    
#        print(dicoNames[Name]["Average Spectra"][4],dicoNames[Name]["Average Correlations"][4])
        
    specfilename = umdfile[:-8] + '.vibr_clust_fast.dat'
    nf = open(specfilename,'w')
    vaffilename = umdfile[:-8] + '.vels_clust_fast.scf.dat'
    nv = open(vaffilename,'w')
    
    header = "vibrational spectrum of "+centralatom+" (central atom) and "+outeratom+" (coordinating atoms)"
    headerVel = "average correlation of atoms per species : "+centralatom+" (central atom) and "+ outeratom+" (coordinating atoms)"
    
    nf.write(header+'\n')
    nv.write(headerVel+'\n')
    
    for species in dicoNames.keys() :
        headerstring = '\n'+species+'\nFrequency(cm^-1)\t'
        headerstringVel = '\n'+species+'\ntime(fs)\tVel_AutoCorr_Cen(fs)\tVel_AutoCorr_Out(fs)'
        headerstring = headerstring + 'VDos_Central \t VDos_Outer \t VDosRad_Outer \t VDosAng_Outer'
                
        nf.write(headerstring+'\n')
        nv.write(headerstringVel+'\n')
        
        frequencies = dicoNames[species]['Average Spectra'][-1]   
        Ang_fft_correlation = dicoNames[species]['Average Spectra'][3]        
        Rad_fft_correlation = dicoNames[species]['Average Spectra'][2]        
        Cen_fft_correlation = dicoNames[species]['Average Spectra'][0]
        Out_fft_correlation = dicoNames[species]['Average Spectra'][1]
        
        avg_Cen_correlation = dicoNames[species]['Average Correlations'][0]
        avg_Out_correlation = dicoNames[species]['Average Correlations'][1]



        for ii in range(len(frequencies)):
            string=str(frequencies[ii]*33356.41)
            stringVel = str(ii*TimeStep)
            string=string+'\t'+str(Cen_fft_correlation[ii])+'\t'+str(Out_fft_correlation[ii])+'\t'+str(Rad_fft_correlation[ii])+'\t'+str(Ang_fft_correlation[ii])
            stringVel = stringVel +'\t' + str(avg_Cen_correlation[ii]) + '\t' +str(avg_Out_correlation[ii])
            nf.write(string+'\n')
            nv.write(stringVel+'\n')
    nf.close()
    nv.close()
#    print("total time : ",time.time()-start," time calc :",tautf-taut,"time merge :", timemergef-timemerge," timemat :",timematf-timemat," timelec :", timelecf-timelec, "timelecf =",timelecb-timelec)
     
    print("Vibrational spectrum file created under the name "+specfilename)        
    print("Correlations values file created under the name "+vaffilename)        


    sys.exit()
if __name__ == "__main__":
    main(sys.argv[1:])
