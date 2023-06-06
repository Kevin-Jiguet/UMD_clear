#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 10:32:24 2023

@author: Kevin Jiguet-Covex
"""
import sys
import numpy as np
import os
import ctypes
from os.path import join

current_path=os.path.abspath(__file__)#For all this to work, the file c_gofr.so must be in the same directory than gofr_umd
path_split=current_path.split('/')
path_red=path_split[1:-1]
path_new=''
for u in path_red:
    path_new+='/'+u

acorr_lib = ctypes.cdll.LoadLibrary(join(path_new, 'c_correlate_H.so'))
acorr_lib.compute_autocorrelation.argtype = [ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.c_int,ctypes.c_int,ctypes.c_int, ctypes.c_int]
acorr_lib.compute_autocorrelation.restype = None

def autocorrel(Tab,firststep,originstep,length):
    
    if(Tab.ndim==1):
        Tab=np.array([Tab])
    elif(Tab.ndim!=2):
        print("The input must be a 1- or 2-dimensional numpy array")
        sys.exit()
        
    tablength=np.shape(Tab)[1]
    tabwidth=np.shape(Tab)[0]
    print(tabwidth)
    ListTab=[Tab[i,:] for i in range(tabwidth)]
    print(ListTab)
    numborigin=0
    ListTabautocorr=[]
    ListTabautocorrsigma=[]
    ListNorm=[]
    ListNumborigin=[]
    
    for tab in ListTab :
        tabautocorr=np.zeros(length)
        tabautocorrsigma=np.zeros(length)
    
        tabp=tab.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        tabautocorrp=tabautocorr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        tabautocorrsigmap=tabautocorrsigma.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        numborigin = tablength-length
        
        acorr_lib.compute_autocorrelation(tabp,tabautocorrp,tabautocorrsigmap,originstep,length,firststep,tablength)
    
        normalization=np.mean(tab[np.arange(firststep,tablength-length,originstep)]**2)
        
        ListTabautocorr.append(tabautocorr.copy())
        ListNorm.append(normalization)
        ListTabautocorrsigma.append(tabautocorrsigma.copy())
        ListNumborigin.append(numborigin)

    return ListTabautocorr, ListTabautocorrsigma, ListNumborigin, ListNorm
