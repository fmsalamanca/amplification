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
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
            'figure.figsize': (16,9),
            'axes.labelsize': 'xx-large',
            'axes.titlesize': 'xx-large',
            'xtick.labelsize': 'xx-large',
            'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

dirID = 21
path  = '/home/fernando/ownCloud/Papers/Fillers/Mechprop/CICLOS/'+str(dirID)+'_LS28649_ciclos.is_tens_RawData/'

def Ac(Nc,Ne):
    a = (3*Ne)/(2*Nc)
    return 0.5 + (np.sqrt(a/np.pi)*np.exp(-a))/erf(np.sqrt(a))

def n(Nc,Ne):
    return Nc*(1-Ac(Nc,Ne))/(2*Ac(Nc,Ne))

rho = 0.92 #g/cm3
T   = 298 #K
R   = 8.3144 #J/(mol K)
Ms  = 131 #g/mol
Na  = 6.022*np.float_power(10,23) #Avogadro constant
kb  = 1.380649*np.float_power(10,-23) #Boltzmann constant J/K

z = 0.387 #From revisiting segmental order...

Dres = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/Dresmax.dat',dtype=float)
Dres = Dres[dirID-1]
Dstatk = 10918*1.4

data     = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/NRCB.dat',dtype='float',delimiter=',')

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

COAN   = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/COAN.dat')
COAN   = COAN[dirID-1]
phi   = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/phiCB.dat')
phi   = phi[dirID-1]
phieff  = (1 + COAN*0.0181)*phi/1.59

Ne  = (1/np.sqrt(6))*rho*R*T/(gedeff*Ms)

def func(Nc):
    return (1.213)*(Dres/Dstatk) - (2*z/Nc)*np.sqrt((Nc+2*n(Nc,Ne))/(0.8*Ne))*np.arctan(Nc/(2*np.sqrt(Nc*n(Nc,Ne)+n(Nc,Ne)**2)))

nue = Na*rho/(Ne*Ms)

Nc    = fsolve(func,x0=1)[0] #Ncc is an array, so we select the component zero which is a float
nuc   = Na*rho/(Nc*Ms)
gcnet = Ac(Nc,Ne)*rho*R*T/(Nc*Ms)
nufr  = (gcdeff - gcnet*(1-phieff))/(kb*T*phieff)

n = Nc*(1-Ac(Nc,Ne))/(2*Ac(Nc,Ne))
m = (2/Nc)*z*np.sqrt((Nc+2*n)/(0.8*Ne))*np.arctan(Nc/(2*np.sqrt(Nc*n+n**2)))

print('')
if True:
    print('nuc')
    print(nuc)
    print('nue')
    print(nue)
    print('Nc')
    print(Nc)
    print('Ne')
    print(Ne)
    print('nufr')
    print(nufr)
    print('Ac')
    print(Ac(Nc=Nc,Ne=Ne))
#saving = bool(input('Hit enter if you do not want to save the data. '))
saving = True
if saving:
    dirdest  = path+'Iterator MAX'

    info = np.asarray([[gcdeff,gcdefferr,gedeff,gedefferr,dd,dderr,A,Aerr,B,Berr,Nc,0,Ne,0,nufr,0,Ac(Nc=Nc,Ne=Ne)]])
    f = open(dirdest+'/dataMAX.dat','a')
    np.savetxt(f,info,delimiter=',')

    f.close()