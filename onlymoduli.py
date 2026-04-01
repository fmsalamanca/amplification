import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit, fsolve
import warnings
from scipy.special import erf, erfc
from scipy.stats import t
import time
import os
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
            'figure.figsize': (16,9),
            'axes.labelsize': 'xx-large',
            'axes.titlesize': 'xx-large',
            'xtick.labelsize': 'xx-large',
            'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)

path = '/home/fernando/ownCloud/Papers/SegmentalOrder/TensileTest/NRP/'
ID = '2'
dirID = '.is_tens_RawData/'

deff = float(input('Defects (%) = '))

def stressfit(l,dd,gc,ge): #the most useful form. Note that here dd = d**2. C gives flexibility for initial stress after corrections in strain axis
    denom  = 1 - (dd)*((l**2) + 2/l - 3)
    factor = l - 1/(l**2)
    term1 = (1 - dd)*factor/(denom**2)
    term2 = dd*factor/denom
    term3 = -l**float(-1 - 1) + l**float(1/2 - 1)
    return gc*(term1 - term2) + 2*(ge/1)*term3

def stressfitNOdd(l,gc,ge): #the most useful form. Note that here dd = d**2. C gives flexibility for initial stress after corrections in strain axis
    factor = l - 1/(l**2)
    term3 = -l**float(-1 - 1) + l**float(1/2 - 1)
    return gc*(factor) + 2*(ge/1)*term3

df_list = []

for files in os.listdir(path+ID+dirID):
    if files.endswith('.csv'):
        df = pd.read_csv(path+ID+dirID+files,header=1,dtype=np.float64,sep=';',encoding='latin-1') 
        df_list.append(df)

reps     = len(df_list)
strain   = {}
stress   = {}

if True: #true if csv coming from our INSTRON
    for i in range(reps):
        #y axis
        stress["stress{0}".format(i)] = df_list[i].values[:,3]
        print('Added stress ',i+1)
        #x axis
        strain["strain{0}".format(i)] = df_list[i].values[:,4]+1
        print('Added strain ',i+1)

fig1 = plt.figure()
for i in range(reps):
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
plt.legend()
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Full Stress-Strain curves')
plt.show()

pleased   = False
sec_value = True     

strainselected = {}
strainaux      = {}
stressselected = {}
stressaux      = {}
maxstresses    = []
strainmax      = []
for i in range(reps):
    strain["strain{0}".format(i)]=strain["strain{0}".format(i)] - strain["strain{0}".format(i)][0]+1#adjust the first point as close as possible to 1
    maxstresses += [max(stress['stress{0}'.format(i)])]
    stress["stress{0}".format(i)]=stress["stress{0}".format(i)] - stress["stress{0}".format(i)][0]#adjust the first point as close as possible to 0
    strainmax += [strain["strain{0}".format(i)][len(strain["strain{0}".format(i)])-1]]

if True:
    lmin = 1
    lmax = float(input('Select the maximum value for unitary deformation to avoid regions where cristalization induced by deformation are more important than elastic behaviour: '))
    lmax = lmax*np.ones(reps)

for i in range(reps):
    strainaux["strain{0}".format(i)]=strain["strain{0}".format(i)][strain["strain{0}".format(i)]<lmax[i]]
    strainselected["strain{0}".format(i)]=strainaux["strain{0}".format(i)][strainaux["strain{0}".format(i)]>lmin]
    stressaux["stress{0}".format(i)]=stress["stress{0}".format(i)][strain["strain{0}".format(i)]<lmax[i]]
    stressselected["stress{0}".format(i)]=stressaux["stress{0}".format(i)][strainaux["strain{0}".format(i)]>lmin]

fig4 = plt.figure()
for i in range(reps):
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
    plt.vlines([lmin,lmax[i]],ymin=0,ymax=5,linestyles='dashed')

plt.xticks(lmax,np.round(lmax,1))
plt.legend()
plt.xlabel(r'$\lambda$')
plt.xlim(left=1)
plt.ylim(bottom=0)
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Stress-Strain curves fully corrected')
#plt.show()

del strainaux
del stressaux
del df_list

ddd = []
gcdefff = []
gedefff = []
RRlist = []
fig5 = plt.figure()
for i in range(reps):
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1)) 
    popt, pcov = curve_fit(stressfit,xdata=strainselected["strain{0}".format(i)],ydata=stressselected["stress{0}".format(i)],bounds=((0,0,0),(0.1,np.inf,np.inf)))
    RR         = 1-np.sum(np.sqrt(np.diag(pcov)))
    plt.plot(strain["strain{0}".format(i)],stressfit(l=strain["strain{0}".format(i)],dd=popt[0],gc=popt[1],ge=popt[2]),label='HSH fitting for sample {0} with R squared = {1}'.format(i+1,np.round(RR,5)),linestyle='dashed')
    plt.vlines([lmin,lmax[i]],ymin=0,ymax=5,linestyles='dashed')
    plt.xticks([lmin,lmax[i]],[lmin,np.round(lmax[i],1)])
    ddd.append(popt[0])
    gcdefff.append(popt[1]/((1-deff/100)))
    gedefff.append(popt[2]/((1-deff/100)))
    RRlist.append(RR)
plt.legend()
plt.xlim(left=1)
plt.ylim(bottom=0,top=max(maxstresses)+1)

plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Fitting data: Gc = '+str(np.round(np.average(gcdefff),2))+' Ge = '+str(np.round(np.average(gedefff),2))+' '+r'$\delta^2 = $'+str(np.round(np.average(ddd),2)))
#plt.show()

print('')

gcdeff = np.average(gcdefff)
gedeff = np.average(gedefff)
dd     = np.average(ddd)

if True:
    print('Gc WITH DEFECTS')
    print(gcdeff)
    print('Ge WITH DEFECTS')
    print(gedeff)
    print('delta2')
    print(dd)
    print('')
    print(np.average(RRlist))