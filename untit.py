import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit, fsolve
import warnings
from scipy.special import erf, erfc
from scipy.stats import t
import time
from tabulate import tabulate
import os

gcdeff = 0.282
gedeff = 0.132


dirID = 22
path  = '/home/fernando/ownCloud/Papers/Fillers/Mechprop/CICLOS/'+str(dirID)+'_LS28649_ciclos.is_tens_RawData/'
#path  = '/home/fernando/Papers/Inhomogeneities_peroxide/Samples/'+str(dirID)+'DCP/'
COAN   = np.loadtxt('COAN.dat')
COAN   = COAN[dirID-1]
phi   = np.loadtxt('phiCB.dat')
phi   = phi[dirID-1]
phieff  = (1 + COAN*0.0181)*phi/1.59
ad    = 1/(1-phieff)
aas    = 1+0.5*phieff+2.2*phieff*phieff
amp   = aas*ad



Dres = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/Mechprop/params.dat',dtype=float)[:,0]



def Ac(Nc,Ne):
    a = (3*Ne)/(2*Nc)
    return 0.5 + (np.sqrt(a/np.pi)*np.exp(-a))/erf(np.sqrt(a))

rho = 0.92 #g/cm3
T   = 298 #K
R   = 8.3144 #J/(mol K)
Ms  = 131 #g/mol
Na  = 6.022*np.float_power(10,23) #Avogadro constant
kb  = 1.380649*np.float_power(10,-23) #Boltzmann constant J/K

Ne    = []
Nc    = []
nue   = []
nuc   = []
gcnet = []
nufr  = []

def n(Nc,Ne):
    return Nc*(1-Ac(Nc,Ne))/(2*Ac(Nc,Ne))

z = 0.387 #from revisiting segmental order
Dstatk = 10918*1.4 #Hz, from UF sample


Ne  = (1/np.sqrt(6))*rho*R*T/(gedeff*Ms)
def func(Nc):
    return (1.213)*(Dres[dirID-1]/Dstatk) - (2*z/Nc)*np.sqrt((Nc+2*n(Nc,Ne))/(0.8*Ne))*np.arctan(Nc/(2*np.sqrt(Nc*n(Nc,Ne)+n(Nc,Ne)**2)))

nue = Na*rho/(Ne*Ms)

Nc    = fsolve(func,x0=1) #Ncc is an array, so we select the component zero which is a float

nuc   = Na*rho/(Nc[0]*Ms)
gcnet = Ac(Nc[0],Ne)*rho*R*T/(Nc[0]*Ms)
nufr  = (gcdeff - gcnet*(1-phieff))/(kb*T*phieff)

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
print('Ac')
print(Ac(Nc,Ne))