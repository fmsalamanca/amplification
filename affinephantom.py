import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit, fsolve
import warnings
from scipy.special import erf, erfc
from scipy.stats import t
import time
import os

dirID = 11
path  = '/home/fernando/ownCloud/Papers/Fillers/Mechprop/CICLOS/'+str(dirID)+'_LS28649_ciclos.is_tens_RawData/'
#path  = '/home/fernando/Papers/Inhomogeneities_peroxide/Samples/'+str(dirID)+'DCP/'

COAN   = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/COAN.dat')
COAN   = COAN[dirID-1]

phi   = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/phiCB.dat')
phi   = phi[dirID-1]

#phieff  = (1 + COAN*0.0181)*phi/1.59

Dres = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/Dresmax.dat',dtype=float)
Dres = Dres[dirID-1]

data = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/NRCBMAX.dat',dtype='float',delimiter=',')
gcdeff = data[dirID-1,0]
gcdefferr = data[dirID-1,1]
gedeff = data[dirID-1,2]
gedefferr = data[dirID-1,3]
dd = data[dirID-1,4]
dderr = data[dirID-1,5]
A = data[dirID-1,6]
Aerr = data[dirID-1,7]
B = data[dirID-1,8]
Berr = data[dirID-1,9]

rho = 0.92 #g/cm3
T   = 298 #K
R   = 8.3144 #J/(mol K)
Ms  = 131 #g/mol
Na  = 6.022*np.float_power(10,23) #Avogadro constant
kb  = 1.380649*np.float_power(10,-23) #Boltzmann constant J/K

Ncunf = rho*R*T/(0.282*Ms)
Neunf  = (1/np.sqrt(6))*rho*R*T/(0.132*Ms)

#######################################################################################################################3
if True:
    Dstatk = (5/3)*208/(1/Ncunf+1/Neunf)
    print('Dstat/k')
    print(Dstatk/1.4)
    print('')

    Ne  = (1/np.sqrt(6))*rho*R*T/(gedeff*Ms)
    Ncinv = (5/3)*Dres/Dstatk - 1/Ne
    Nc = 1/Ncinv

    nuc   = Na*rho/(Nc*Ms)
    nue   = Na*rho/(Ne*Ms)
    gcnet = rho*R*T/(Nc*Ms)
    nufr  = (gcdeff - gcnet*(1-phi))/(kb*T*phi)

    print('nufr')
    print(nufr)
    print('Nc')
    print(Nc)
    print('nuc')
    print(nuc)
    print('Ne')
    print(Ne)
    print('nue')
    print(nue)
    print('')

    dirdest  = '/home/fernando/ownCloud/Papers/Fillers/Affine'
    if not os.path.exists(dirdest):
            os.makedirs(dirdest)
    Ncerr = Neerr = nufrerr = 0
    Acav = 1
    info = np.asarray([[gcdeff,gcdefferr,gedeff,gedefferr,dd,dderr,A,Aerr,B,Berr,Nc,Ncerr,Ne,Neerr,nufr,nufrerr,Acav]])
    f = open(dirdest+'/'+str(dirID)+'-NRCBMAXAF.dat','a')
    np.savetxt(f,info,delimiter=',')

    f.close()

###################################################################################################################################

Ncunf = 0.5*Ncunf
Dstatk = 2*(5/3)*208/(1/Ncunf+1/Neunf)
print('Dstat/k')
print(Dstatk/1.4)
print('')

Ne  = (1/np.sqrt(6))*rho*R*T/(gedeff*Ms)
Ncinv = 2*(5/3)*Dres/Dstatk - 1/Ne
Nc = 1/Ncinv

nuc   = Na*rho/(Nc*Ms)
nue   = Na*rho/(Ne*Ms)
gcnet = 0.5*rho*R*T/(Nc*Ms)
nufr  = (gcdeff - gcnet*(1-phi))/(kb*T*phi)

print('nufr')
print(nufr)
print('Nc')
print(Nc)
print('nuc')
print(nuc)
print('Ne')
print(Ne)
print('nue')
print(nue)

dirdest  = '/home/fernando/ownCloud/Papers/Fillers/Phantom'
if not os.path.exists(dirdest):
        os.makedirs(dirdest)
Ncerr = Neerr = nufrerr = 0
Acav = 0.5
info = np.asarray([[gcdeff,gcdefferr,gedeff,gedefferr,dd,dderr,A,Aerr,B,Berr,Nc,Ncerr,Ne,Neerr,nufr,nufrerr,Acav]])
f = open(dirdest+'/'+str(dirID)+'-NRCBMAXPH.dat','a')
np.savetxt(f,info,delimiter=',')

f.close()