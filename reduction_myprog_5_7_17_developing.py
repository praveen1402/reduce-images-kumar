#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 23:17:50 2017

@author: praveen
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import pylab
from scipy.integrate import trapz, simps
from mpl_toolkits.mplot3d import Axes3D
import time
import os
import io   
import anatools as ana  
from scipy.ndimage.filters import gaussian_filter
#%%

start_time = time.time()
Signal = np.zeros((900,900))
Back = np.zeros((900,900))
prefix= '/Users/praveen/Desktop/Asdf1/'
list = os.listdir('/Users/praveen/Desktop/Asdf1/')
number_files = len(list)
print('the number of files in', prefix, 'are',number_files)
A = np.zeros((1,1000))
k = 0
l = 0
TOF5 = [] 


# This part print the files from the alphanumeric order

fnames = os.listdir('/Users/praveen/Desktop/Asdf1/')
#fnames = os.listdir('/Users/praveen/Documents/Study/Pravin work/Concerning Fel/h5files/sideband/Run_036/rawdata/')
for j in range(4,len(list)):#len(list)):#File number in the folder
    def foo(path):
        fnames[j]   
#    print(fnames[j])
    with h5py.File(prefix+fnames[j], 'r') as f:
        for i in range(0,100):#shot numbers
             
             TOF=f['digitizer/channel1'][(i)]#shot number in the file
             sum_0 = np.sum(TOF[0:5050]/5050)
             sum1 = np.sum(sum_0-TOF[5051:9000])
             Img = f['vmi/andor'][(i)]#VMI Image Matrix:Look at the Camera size(910x910)
             BP = f['Background_Period'][(...)]#Background Period
             BNO = f['bunches'][(i)]#Bunch Number
             TOF5 = sum_0-TOF[5000:8000]
             if BNO%BP!=0:# definetion for background: If BNO%BP==0, then it is background
#             if np.max(TOF5[1127:1190])>=10:
                 k +=1
                 Signal += Img
             else:
                 l +=1
                 Back += Img
#

print(k,l)

plt.matshow(Signal)
plt.pause(1)
plt.colorbar()

plt.matshow(Back)
plt.pause(1)
plt.colorbar()

plt.matshow(Signal-Back)
plt.pause(1)
plt.colorbar()


#%%
# This lines are for writing the data to a .hf file

with h5py.File('/Users/praveen/Desktop/Asdf1/reduced123.h5', 'w') as f:
     f.create_dataset("vmi1/andor1",  data=np.absolute(Signal/k-Back/l))

     
#%%
     
with h5py.File('/Users/praveen/Desktop/Asdf1/reduced123.h5', 'r') as f:
     Img11 = f['vmi1/andor1'][()]

     
     
plt.matshow(Img11)
plt.pause(1)
plt.colorbar()
#




