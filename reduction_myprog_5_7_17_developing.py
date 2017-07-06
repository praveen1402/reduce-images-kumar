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
##A = np.zeros((1,1000))
Img1 = np.zeros((900,900))
Img2 = np.zeros((900,900))
prefix= '/Users/praveen/Desktop/Asdf1/'
list = os.listdir('/Users/praveen/Desktop/Asdf1/')
#prefix= '/Users/praveen/Documents/Study/Pravin work/Concerning Fel/h5files/sideband/Run_036/rawdata/'
#list = os.listdir('/Users/praveen/Documents/Study/Pravin work/Concerning Fel/h5files/sideband/Run_036/rawdata/') # dir is your directory path
number_files = len(list)
print('the number of files in', prefix, 'are',number_files)

A = np.zeros((1,1000))
k = 0
l = 0
Inten = []
maximaTOF = []
SumTOFMAT = np.zeros((1,1000))
SumTOFMATtra = np.zeros((1000,1))
SumTOFMATmul = np.zeros((1000,1000))
SumITOFMAT = np.zeros((1,1000))
SumITOFMATtra = np.zeros((1000,1))
SumITOFMATmul = np.zeros((1000,1000))
SumInten = 0.0
SumIntenSq = 0.0
TOF6 = np.zeros((1,3000))
TOF5 = [] 


# This part print the files from the alphanumeric order
fnames = os.listdir('/Users/praveen/Desktop/Asdf1/')
#fnames = os.listdir('/Users/praveen/Documents/Study/Pravin work/Concerning Fel/h5files/sideband/Run_036/rawdata/')
for j in range(7,8):#len(list)):#File number in the folder
    def foo(path):
        fnames[j]   
#    print(fnames[j])
    with h5py.File(prefix+fnames[j], 'r') as f:
        for i in range(18,19):#shot numbers
             
             TOF=f['digitizer/channel1'][(i)]#shot number in the file
             sum_0 = np.sum(TOF[0:5050]/5050)
             sum1 = np.sum(sum_0-TOF[5051:9000])
             Img = f['vmi/andor'][(i)]
             BP = f['Background_Period'][(...)]
#             print(sum1)
             TOF5 = sum_0-TOF[5000:8000]
             print(np.sum(TOF))
             plt.hist(TOF,60000)
#             print(np.sum(TOF5))
#             TOF6 += TOF5.reshape(1,3000)
#             TOF5 = TOF6.reshape(3000,)
#             TOF7 = TOF5[1000:2000]                
#             plt.plot(TOF5)
#             plt.pause(0.5)
#             Inten.append(np.sum(TOF5[1127:1190])
##             Intensity = f['photon_diagnostics/FEL02/I0_monitor/iom_sh_a'][(i)]
#             if np.max(TOF5[1127:1190])>=3.2:
#                 k +=1
#                 Img1 += Img
##                 TOF5 = sum_0-TOF[5000:8000]
###             print(np.sum(TOF5))
##                 TOF6 += TOF5.reshape(1,3000)
#                 
#             else:
#                 l +=1
#                 Img2 = Img
#                 TOF5 = sum_0-TOF[5000:8000]
#             print(np.sum(TOF5))
#                 TOF6 += TOF5.reshape(1,3000)

#print(np.sum(TOF6))
#print(Inten)
#        print(BP)
#        print(k,l)  
#print(np.max(TOF5),np.min(TOF5),np.mean(TOF5))               
#Img3 = Img1-(k/l)*Img2
#                 
##
#plt.matshow(Img1)
#plt.matshow(Img2)
#plt.matshow(Img1*(1/k)-Img2*(1/l))

#%%
with h5py.File('/Users/praveen/Desktop/Asdf1/reduced.h5', 'w') as f:
     f.create_dataset("vmi/andor",  data=Img1-Img2)
# This lines are for writing the data to a .hf file
#%%
#number_files = len(list)
#print('the number of files in', prefix, 'are',number_files)
#for j in range(1,len(list)):#File number in the folder
#    def foo(path):
#        fnames[j]
#    print(fnames[j])
     
with h5py.File('/Users/praveen/Desktop/Asdf1/reduced3.h5', 'r') as f:
     Img = f['vmi/andor'][()]
#print(Img.shape)         
#plt.matshow(Img)
#plt.colorbar() 
plt.plot(np.diagonal(Img1-Img2))
#filename = '/Users/praveen/Desktop/Asdf1/reduced.h5'
#key = 'vmi/andor'
#blur = 0
#  
#
#with h5py.File(filename, 'r') as f:
#     img = f['vmi/andor'][...]
#     img2 = img/img
#origin = ana.Hist(gaussian_filter(img, blur).T)  # shape=(x,z)
#
#plt.figure(figsize=(8, 8))
#plt.pcolormesh(origin.hist.T)
#plt.colorbar() 

#A = np.diagonal(Img1-(67/33)*Img2)
#plt.plot(A)
#pylab.xlim([0,910])
#pylab.ylim([0,910]) 
#plt.colorbar()
#plt.title('Covariance Map:Cov(Y,X)')
##
#plt.matshow(CovIXY)
#pylab.xlim([0,1000])
#pylab.ylim([0,1000]) 
#plt.colorbar()
#plt.title('Intensity Covariance Map:Cov(Y,I)*Cov(X,I)/Cov(I,I)')
####
#plt.matshow(PCovXY)
#pylab.xlim([0,1000])
#pylab.ylim([0,1000]) 
#plt.colorbar()
#plt.title('Partial Covariance Map:Cov(Y,X)-Cov(Y,I)*Cov(X,I)/Cov(I,I)') 

#A = np.diagonal(PCovXY)
#
#
#    
#plt.plot(A)    








