import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit, fsolve
import warnings
from scipy.special import erf, erfc
from scipy.stats import t
import time
import os

gcdeff = float(input('Gc (MPa) = '))
gedeff = float(input('Ge (MPa) = '))
Dres = float(input('Dres (avg) in Hz = '))

rho = 0.92 #g/cm3
T   = 353 #K
R   = 8.3144 #J/(mol K)
Ms  = 131 #g/mol
Na  = 6.022*np.float_power(10,23) #Avogadro constant
kb  = 1.380649*np.float_power(10,-23) #Boltzmann constant J/K

#######################################################################################################################3

Ncaf = rho*R*T/(gcdeff*Ms)
Ne  = (1/np.sqrt(6))*rho*R*T/(gedeff*Ms)
Dstatk = (5/3)*Dres/(1/Ncaf+1/Ne)

m = 1/Ncaf+1/Ne
print('AFFINE MODEL')
print('Dstat/k')
print(Dstatk/1.4)
print('')
print('Nc')
print(Ncaf)
print('Ne')
print(Ne)
print('')
print('m')
print(m)
print('')
###################################################################################################################################

Ncph = 0.5*Ncaf
Dstatk = 2*(5/3)*Dres/(1/Ncph+1/Ne)
print('PHANTOM MODEL')
print('Dstat/k')
print(Dstatk/1.4)
print('')

m = 0.5*(1/Ncph+1/Ne)

print('Nc')
print(Ncph)
print('Ne')
print(Ne)
print('')
print('m')
print(m)