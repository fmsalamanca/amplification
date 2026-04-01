import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import os.path
import time as t
import shutil
import subprocess
from scipy import stats
import inspect
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
            'figure.figsize': (16,9),
            'axes.labelsize': 'xx-large',
            'axes.titlesize': 'xx-large',
            'xtick.labelsize': 'xx-large',
            'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

Dstatk = 6300
z = 0.387
Ac = 0.67
Ne = 25
dist = np.loadtxt('bestdistribution_22.dat')
pop = dist[:,0]
Dres = dist[:,1]
print(np.sum(pop[1]*Dres))
print('NORMALIZED')

maxdres = Dres[Dres.argmax()]
maxdres = pop[Dres.argmax()]

fig1 = plt.figure()
plt.plot(pop,Dres,'k')
plt.xlim(left=0,right=1)
plt.ylim(bottom=0)
plt.xlabel(r'$D_{res} (kHz)$')
plt.ylabel('Frequency')
plt.text(0.75*plt.xlim()[1],0.8*plt.ylim()[1],r'$D^*_{res} = %.3f\,\,Hz$' % (maxdres*1000))
plt.title(r'$D_{res}$'+' distribution')
plt.show()

m = (1.213/Dstatk)*Dres
print(np.sum(pop[1]*m))
print('NOT NORMALIZED')

fig2 = plt.figure()
plt.plot(pop,m,'b')
plt.xlim(left=0,right=1)
plt.ylim(bottom=0)
plt.xlabel(r'$m$')
plt.ylabel('Frequency')
plt.title(r'$m$'+' distribution')
plt.show()


distrib = stats.norm
params = distrib.fit()
prefac = 2*z/(0.8*Ne*Ac)
postfac = np.pi/2 - 2*np.sqrt(((1-Ac)/2*Ac)+((1-Ac)/2*Ac)**2)
Nc = prefac**2*(1/m**2)*postfac**2

fig3 = plt.figure()
plt.plot(pop,Nc,'g')
plt.xlim(left=0,right=1)
plt.ylim(bottom=0)
plt.xlabel(r'$m$')
plt.ylabel('Frequency')
plt.title(r'$N_c$'+' distribution')
plt.show()