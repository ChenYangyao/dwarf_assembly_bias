# coding: utf-8

# Copyright (C) 2025 Hui-Jie Hu (huhuijienao@gmail.com) - All Rights Reserved
# You may use, distribute and modify this code under the MIT license. We kindly
# request you to give credit to the original author(s) of this code, and cite 
# the paper, Qi Guo et al. Nature Astronomy 4, 246–251 (2020) if you use this 
# code in your research.

# In[15]:

import numpy as np
from math import *
import matplotlib
import csv
import time

Gc=4.3e-6 #(km/s)^2*kpc/M_sun
h=0.7
rhocri=2.8e2*h**2 #M_sun/kpc^3
cctn = 13
pc=3.0857e16

def Mhaldrive(Rmax,Vmax,Ac=None,Rac=None): #re kpc Vmax km/s
    c00 = cctn #13 primary
    Mc=log(np.power(Vmax,2)*Rmax/Gc,10)
    if Ac == None:
        Ac = 0.01
    if Rac == None:
        Rac = 0.005
    c1 = 1
    while c1 >= Ac :
        c0 = c00
        rhos=200.0/3.0*rhocri*c0**3/(log(1+c0)-c0/(1+c0))
        xs_max = 1000 #Rs max 1000kpc
        xs_min = 0
        xs = 500
        while (xs_max-xs_min) > Rac:
            MM = log(4*pi*rhos*(xs**3)*(log((xs+Rmax)/xs)-Rmax/(xs+Rmax)),10)
            if MM < Mc:
                xs_min = xs
            else:
                xs_max = xs
            xs = (xs_min+xs_max)/2
        Rs = xs
        Mhalo = log(4.0/3.0*pi*rhocri*200.0*c0**3*np.power(Rs,3),10)
        c00 = 10**(0.905-0.101*(Mhalo-12+log(h,10)))
        c1 = abs(c00 - c0)
    return Mhalo

def Mhaldrive_Burk(Rmax,Vmax,Rac=None): #re kpc Vmax km/s
    Mc=np.log10(np.power(Vmax,2)*Rmax/Gc)
    r00 = 0.65 #primary
    r0 = r00
    if Rac == None:
        Rac = 0.0005
    r1 = 100
    while r1 > Rac:
        rho0 = 10**Mc/(2*pi*(r0**3)*(np.log((r0+Rmax)/r0)-np.arctan(Rmax/r0)+0.5*np.log(1+(Rmax/r0)**2)))
        Rvir = (3/(800*pi*rhocri)*10**(np.log10(r0)/0.58-0.66/0.58+11))**(1/3)
        Mvir = np.log10(2*pi*rho0*(r0**3)*(log((r0+Rvir)/r0)-np.arctan(Rvir/r0)+0.5*log(1+(Rvir/r0)**2)))
        if abs(10**(0.66+0.58*(Mvir-11))-r0) <= r1:
            r1 = abs(10**(0.66+0.58*(Mvir-11))-r0)
            r0 = 10**(0.66+0.58*(Mvir-11))
        else:
            break
    return Mvir

def Mhaldrive_Pseudo(Rmax,Vmax,Rac=None): #re kpc Vmax km/s
    Mc=np.log10(np.power(Vmax,2)*Rmax/Gc)
    r00 = 1 #primary
    r0 = r00
    rho0 = 10**(3.8175/0.47-0.41/0.47*np.log10(r0))
    if Rac == None:
        Rac = 0.0005
    rma = 15
    rmi = 0.001
    r1 = 100
    while r1 > Rac:
        MM = np.log10(4*pi*rho0*(r0**3)*(Rmax/r0-np.arctan(Rmax/r0)))
        if MM > Mc:
            rma = r0
        else:
            rmi = r0
        r0 = (rma+rmi)/2
        rho0 = 10**(3.8175/0.47-0.41/0.47*np.log10(r0))
        r1 = rma-rmi
    
    Mrh = 800*pi/3*rhocri
    Rvma = 1000
    Rvmi = Rmax
    Rvir = 10*Rmax
    R1 = 10
    while R1 >10*Rac:
        MrhM = 4*pi*rho0*(r0**3)*(Rvir/r0-np.arctan(Rvir/r0))/Rvir**3
        if MrhM > Mrh:
            Rvmi = Rvir
        else:
            Rvma = Rvir
        Rvir = (Rvma+Rvmi)/2
        R1 = Rvma-Rvmi
    Mvir = np.log10(800*pi/3*rhocri*Rvir**3)
    return Mvir