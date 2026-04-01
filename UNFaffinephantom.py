import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit, fsolve
import warnings
from scipy.special import erf, erfc
from scipy.stats import t
import time
import os

Dres = np.loadtxt('/home/fernando/ownCloud/Papers/SegmentalOrder/NMR_NRP-NRS.dat',dtype=float)[:,0]

data = np.loadtxt('/home/fernando/ownCloud/Papers/SegmentalOrder/NRP-NRS.dat',dtype='float',delimiter=',')

Gc = data[:,0]

Ge = data[:,2]

rho = 0.92 #g/cm3
T   = 298 #K
R   = 8.3144 #J/(mol K)
Ms  = 131 #g/mol
Na  = 6.022*np.float_power(10,23) #Avogadro constant
kb  = 1.380649*np.float_power(10,-23) #Boltzmann constant J/K

Ncaf = rho*R*T/(Gc*Ms)
Ncph = 0.5*rho*R*T/(Gc*Ms)
Ne = (1/(np.sqrt(6)))*rho*R*T/(Ge*Ms)

m=[]
Dstatk=[]
print('AFFINE GOING ON')
#######################################################################################################################3
for i in range(len(Dres)):
    m += [1/Ncaf[i]+1/Ne[i]]
    Dstatk += [(5/3)*Dres[i]/(m[i]*1.4)]

dirdest  = '/home/fernando/ownCloud/Papers/SegmentalOrder/Affine'
if not os.path.exists(dirdest):
        os.makedirs(dirdest)
Acav = np.ones(len(Dres))

print('Dstat/k')
print(Dstatk)
print('')

print('Nc')
print(Ncaf)

print('Ne')
print(Ne)

print('')

info = np.array([Ncaf,m,Dstatk,Acav]).T
f = open(dirdest+'/affine.dat','a')
np.savetxt(f,info,delimiter=',')

f.close()

###################################################################################################################################
m=[]
Dstatk=[]
print('PHANTOM GOING ON')
for i in range(len(Dres)):
    m += [0.5*(1/(Ncph[i])+1/Ne[i])]
    Dstatk += [(5/3)*Dres[i]/(m[i]*1.4)]

dirdest  = '/home/fernando/ownCloud/Papers/SegmentalOrder/Phantom'
if not os.path.exists(dirdest):
        os.makedirs(dirdest)
Acav = 0.5*np.ones(len(Dres))

print('Dstat/k')
print(Dstatk)
print('')

print('Nc')
print(Ncph)

print('Ne')
print(Ne)

print('')

info = np.asarray([Ncph,m,Dstatk,Acav]).T
f = open(dirdest+'/phantom.dat','a')
np.savetxt(f,info,delimiter=',')

f.close()