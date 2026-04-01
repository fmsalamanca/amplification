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

dirID = 'NR_rCBVL'
path  = '/home/fernando/ownCloud/Papers/ICB/Mechprop/'+str(dirID)+'_ciclos.is_tens_RawData/'
#path  = '/home/fernando/Papers/Inhomogeneities_peroxide/Samples/'+str(dirID)+'DCP/'


def stressfit(l,dd,gc,ge): #simplest form from Heinrich, Sommer and Helmis model
    denom  = 1 - (dd)*((l**2) + 2/l - 3)
    factor = l - 1/(l**2)
    term1 = (1 - dd)*factor/(denom**2)
    term2 = dd*factor/denom
    term3 = -l**float(-1 - 1) + l**float(1/2 - 1)
    return gc*(term1 - term2) + 2*(ge/1)*term3

df_list = []
for files in os.listdir(path):
    if files.endswith('.csv'):
        df = pd.read_csv(path+files,header=1,dtype=np.float64,sep=';') 
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

if False: #true if csv coming from our hands
    for i in range(reps):
        #y axis
        stress["stress{0}".format(i)] = df_list[i].values[:,0]
        print('Added stress ',i+1)
        #x axis
        strain["strain{0}".format(i)] = df_list[i].values[:,1]+1
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
if sec_value:
    while not pleased:
        precycling = input('Did you see precycling on the curves? Hit Enter if not. ')
        if precycling: #true if precycling
            points = int(input('How many points do you want to remove? '))
            fig2 = plt.figure()
            for i in range(reps):
                    #y axis
                    stress["stress{0}".format(i)] = df_list[i].values[points:len(df_list[i]),3]
                    print('Precycling removed from stress ',i+1)
                    #x axis
                    strain["strain{0}".format(i)] = df_list[i].values[points:len(df_list[i]),4]+1
                    print('Precycling removed from strain ',i+1)
                    plt.scatter(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
            plt.legend()
            plt.xlabel(r'$\lambda$')
            plt.ylabel(r'$\sigma$'+' (MPa)')
            plt.title('Stress-Strain curves without precycling')
            plt.show()
            pleased = input('Was precycling removed correctly? Hit Enter if not. ')
        else:
            pleased = True

if not sec_value:
    for i in range(reps):
                    #y axis
        stress["stress{0}".format(i)] = df_list[i].values[2900:len(df_list[i]),3]
        print('Precycling removed from stress ',i+1)
                    #x axis
        strain["strain{0}".format(i)] = df_list[i].values[2900:len(df_list[i]),4]+1
        print('Precycling removed from strain ',i+1)
        

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
if False: #true to create a csv file with the experimental curve corrected
    pd.DataFrame(np.array([strain["strain0"],stress["stress0"]]).T).to_csv('/home/fernando/Papers/Fillers/Mechprop/CICLOS/curvas_exp/'+"curva_exp"+str(dirID)+".csv",header=None)
    print('File created')

fig3 = plt.figure()
for i in range(reps):
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
plt.legend()
plt.xlim(left=1)
plt.ylim(bottom=0)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Stress-Strain curves corrected')
plt.show()


lmin = 1
lmax = 0.521*np.array(strainmax)

if False:
    lmin = float(input('Select the minimum value for unitary deformation to avoid regions where precycle effect is more important than elastic behaviour: '))
    lmax = float(input('Select the maximum value for unitary deformation to avoid regions where cristalization induced by deformation are more important than elastic behaviour: '))

for i in range(len(lmax)):
    print('Optimal interval: '+'('+str(lmin)+','+str(np.round(lmax[i],3))+')')

for i in range(reps):
    strainaux["strain{0}".format(i)]=strain["strain{0}".format(i)][strain["strain{0}".format(i)]<lmax[i]]
    strainselected["strain{0}".format(i)]=strainaux["strain{0}".format(i)][strainaux["strain{0}".format(i)]>lmin]
    stressaux["stress{0}".format(i)]=stress["stress{0}".format(i)][strain["strain{0}".format(i)]<lmax[i]]
    stressselected["stress{0}".format(i)]=stressaux["stress{0}".format(i)][strainaux["strain{0}".format(i)]>lmin]

fig4 = plt.figure()
for i in range(reps):
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
    plt.vlines([lmax[i]],ymin=0,ymax=5,linestyles='dashed')
plt.xticks(lmax,np.round(lmax,1))
plt.legend()
plt.xlabel(r'$\lambda$')
plt.xlim(left=1)
plt.ylim(bottom=0)
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Stress-Strain curves fully corrected')
plt.show()

del strainaux
del stressaux
del df_list

Dres = float(input('Write here Dres measured by DQ-NMR in Hz. '))
deff = float(input('Write here the fraction of defects measured by DQ-NMR in %. '))
Dstatk = float(input('Write here single-chain static frequency for the network in Hz. '))
Dstatk = Dstatk*1.4
dd     = []
gcdeff = []
gedeff = []
RRlist = []

fig5 = plt.figure()
for i in range(reps):
    popt, pcov = curve_fit(stressfit,xdata=strainselected["strain{0}".format(i)],ydata=stressselected["stress{0}".format(i)],bounds=((0,0,0),(0.1,np.inf,np.inf)))
    RR         = 1-np.sum(np.sqrt(np.diag(pcov)))
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
    plt.plot(strain["strain{0}".format(i)],stressfit(l=strain["strain{0}".format(i)],dd=popt[0],gc=popt[1],ge=popt[2]),label='HSH fitting for sample {0} with R squared = {1}'.format(i+1,np.round(RR,5)),linestyle='dashed')
    dd.append(popt[0])
    gcdeff.append(popt[1]/((1-deff/100)))
    gedeff.append(popt[2]/((1-deff/100)))
    RRlist.append(RR)
plt.legend()

plt.xlim(left=1)
plt.ylim(bottom=0,top=max(maxstresses)+1)
plt.vlines([lmax[i]],ymin=0,ymax=5,linestyles='dashed')
plt.xticks([lmax[i]],[np.round(lmax[i],1)])
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Fitting data: Gc = '+str(np.round(np.average(gcdeff),2))+' Ge = '+str(np.round(np.average(gedeff),2))+' '+r'$\delta^2 = $'+str(np.round(np.average(dd),2)))
plt.show()

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

def n(Nc,Ne):
    return Nc*(1-Ac(Nc,Ne))/(2*Ac(Nc,Ne))

z = 0.387 #from revisiting segmental order

for i in range(reps):
    Nee  = (1/np.sqrt(6))*rho*R*T/(gedeff[i]*Ms)
    def func1(Nc):
        aux = (1-Ac(Nc,Nee))/(2*Ac(Nc,Nee))
        return (1.213)*(Dres/Dstatk) - ((2*z)/np.sqrt(Ac(Nc,Nee)*0.8*Nee*Nc))*(np.pi/2-2*np.sqrt(aux+aux**2))

    def func(Nc):
        return (1.213)*(Dres/Dstatk) - (2*z/Nc)*np.sqrt((Nc+2*n(Nc,Nee))/(0.8*Nee))*np.arctan(Nc/(2*np.sqrt(Nc*n(Nc,Nee)+n(Nc,Nee)**2)))

    nuee = Na*rho/(Nee*Ms)
    Ne.append(Nee)
    nue.append(nuee)

    Ncc    = fsolve(func,x0=1) #Ncc is an array, so we select the component zero which is a float
    Nc.append(Ncc[0])
    nucc   = Na*rho/(Ncc[0]*Ms)
    nuc.append(nucc)

if True:
    x = np.arange(0,1000)
    plt.plot(func(x),label='exact form')
    plt.axhline(y=0)
    plt.plot(func1(x),label='Taylor form')
    plt.ylabel('k-m(Nc)')
    plt.xlabel('Nc')
    plt.title('Comparing segmental order vector parameter approaches')
    plt.legend()
    plt.show()

Ne    = np.asarray(Ne)
Nc    = np.asarray(Nc)
nue   = np.asarray(nue)
nuc   = np.asarray(nuc)

gcdefferr  = t.ppf(0.95,df=len(gcdeff)-1)*np.std(gcdeff)/np.sqrt(len(gcdeff))
gedefferr  = t.ppf(0.95,df=len(gedeff)-1)*np.std(gedeff)/np.sqrt(len(gedeff))
dderr  = t.ppf(0.95,df=len(dd)-1)*np.std(dd)/np.sqrt(len(dd))

Ncerr  = t.ppf(0.95,df=len(Nc)-1)*np.std(Nc)/np.sqrt(len(Nc))
Neerr  = t.ppf(0.95,df=len(Ne)-1)*np.std(Ne)/np.sqrt(len(Ne))
nucerr = t.ppf(0.95,df=len(nuc)-1)*np.std(nuc)/np.sqrt(len(nuc))
nueerr = t.ppf(0.95,df=len(nue)-1)*np.std(nue)/np.sqrt(len(nue))

Acerr  = t.ppf(0.95,df=len(Ac(Nc,Ne)-1)*np.std(Ac(Nc,Ne))/np.sqrt(len(Ac(Nc,Ne))))

gcdeff = np.average(gcdeff)
gedeff = np.average(gedeff)
dd   = np.average(dd)

Acav = np.average(Ac(Nc,Ne))
Nc   = np.average(Nc)
Ne   = np.average(Ne)
nuc  = np.average(nuc)
nue  = np.average(nue)

if True:
    print('Gc')
    print(gcdeff,' +- ',gcdefferr)
    print('Ge')
    print(gedeff,' +- ',gedefferr)
    print('delta2')
    print(dd,' +- ',dderr)
    print('Nc')
    print(Nc,' +- ',Ncerr)
    print('Ne')
    print(Ne,' +- ',Neerr)
    print('Ac')
    print(Acav,' +- ',Acerr)
    print('nuc')
    print(nuc,' +- ',nucerr)
    print('nue')
    print(nue,' +- ',nueerr)

dirdest  = path+'Analysis'
if not os.path.exists(dirdest):
    os.makedirs(dirdest)
else:
    raise('You should change the directory name for this analysis')
fig1.savefig(dirdest+'/Figure1.pdf')
if precycling:
    fig2.savefig(dirdest+'/Figure2.pdf')
fig3.savefig(dirdest+'/Figure3.pdf')
fig4.savefig(dirdest+'/Figure4.pdf')
fig5.savefig(dirdest+'/Figure5.pdf')

info = np.asarray([[gcdeff,gcdefferr,gedeff,gedefferr,dd,dderr,Nc,Ncerr,Ne,Neerr,Acav]])

np.savetxt(dirdest+'/data.dat',info,delimiter=',')