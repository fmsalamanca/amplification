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

#dirID = 'NR20'
#path  = '/home/fernando/ownCloud/Papers/Mc distribution/Datos de Antonio/Instron_DataFiles/'+str(dirID)+'.is_tens_RawData/'
#path  = '/home/fernando/ownCloud/Papers/TubeModelCNT/'+dirID+'/'
path = '/home/fernando/Nextcloud/Papers/SegmentalOrder/TensileTest/NRS/'
dirID = 'NRS37NOMasticado'
#dirID = '.is_tens_RawData/'

path = 'C:\\Users\\ferna\\Nextcloud\\Papers\\Entanglements\\Masticados-instron\\'+str(dirID)+'.is_tens_RawData\\'

auto = True
instron = True
print('')
print('Sample: '+str(dirID))
print('')
if auto:
    Dres = float(input('1H DQ-NMR test. Dres (Hz): '))

    deff = float(input('1H DQ-NMR test. Fraction of defects (%): '))

def stressfit(l,dd,gc,ge): #the most useful form. Note that here dd = d**2. C gives flexibility for initial stress after corrections in strain axis
    denom  = 1 - (dd)*((l**2) + 2/l - 3)
    factor = l - 1/(l**2)
    term1 = (1 - dd)*factor/(denom**2)
    term2 = dd*factor/denom
    term3 = -l**float(-1 - 1) + l**float(1/2 - 1)
    return gc*(term1 - term2) + 2*(ge/1)*term3

def stressfit_b(l,d,b,gc,ge,C): #too general because hardly ever b != 1
    denom  = 1 - (d**2)*((l**2) + 2/l - 3)
    factor = l - 1/(l**2)
    term1 = (1 - d**2)*factor/denom
    term2 = d**2*factor/(denom**2)
    term3 = -l**float(-b - 1) + l**float(b/2 - 1)
    return gc*(term1 - term2) + 2*(ge/b)*term3 + C

def stressfitNOdd(l,gc,ge): #the most useful form. Note that here dd = d**2. C gives flexibility for initial stress after corrections in strain axis
    factor = l - 1/(l**2)
    term3 = -l**float(-1 - 1) + l**float(1/2 - 1)
    return gc*(factor) + 2*(ge/1)*term3

df_list = []

for files in os.listdir(path):
    if files.endswith('.csv'):
        df = pd.read_csv(path+files,header=1,dtype=np.float64,sep=';',encoding='latin-1') 
        df_list.append(df)
'''
for files in os.listdir(path):
    if files.endswith('.csv'):
        df = pd.read_csv(path+files,dtype=np.float64,sep='\t') 
        df_list.append(df)
        print('Done!')

for files in os.listdir(path):
    if files.endswith('.xlsx'):
        df = pd.read_excel(path+files) 
        df_list.append(df)
'''
reps     = len(df_list)
strain   = {}
stress   = {}

if instron: #true if csv coming from our INSTRON
    for i in range(reps):
        #y axis
        stress["stress{0}".format(i)] = df_list[i].values[:,3]
        print('Added stress ',i+1)
        print(stress["stress{0}".format(i)])
        #x axis
        strain["strain{0}".format(i)] = df_list[i].values[:,4]+1
        print('Added strain ',i+1)
        print(strain["strain{0}".format(i)])
if not instron: #true if csv coming from our hands
    for i in range(reps):
        #y axis
        stress["stress{0}".format(i)] = df_list[i].values[:,1]
        print('Added stress ',i+1)
        #x axis
        strain["strain{0}".format(i)] = df_list[i].values[:,0]/100+1
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
        #precycling = False
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

if not precycling:
    points = 0
        

strainselected = {}
strainaux      = {}
stressselected = {}
stressaux      = {}
maxstresses    = []
strainmax      = []
for i in range(reps):
    strain["strain{0}".format(i)]=strain["strain{0}".format(i)] - strain["strain{0}".format(i)][0]+1#adjust the first point as close as possible to 1
    stress["stress{0}".format(i)]=stress["stress{0}".format(i)] - stress["stress{0}".format(i)][0]#adjust the first point as close as possible to 0
    maxstresses += [max(stress['stress{0}'.format(i)])]
    strainmax += [strain["strain{0}".format(i)][len(strain["strain{0}".format(i)])-1]]
if False: #true to create a csv file with the experimental curve corrected
    pd.DataFrame(np.array([strain["strain0"],stress["stress0"]]).T).to_csv('/home/fernando/Papers/Fillers/Mechprop/CICLOS/curvas_exp/'+"curva_exp"+str(dirID)+".csv",header=None)
    print('File created')
print('')
print(strainmax)
print(maxstresses)
print('')
fig3 = plt.figure()
for i in range(reps):
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
plt.legend()
plt.xlim(left=1)
plt.ylim(bottom=0)
plt.text(0.65*plt.xlim()[1],0.7*plt.ylim()[1],"%i points removed" % points,{'fontsize': 15})
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Stress-Strain curves corrected')
plt.show()

if not auto:
    lmin = 1
    lmax = float(input('Select the maximum value for unitary deformation to avoid regions where cristalization induced by deformation are more important than elastic behaviour: '))
    lmax = lmax*np.ones(reps)
else:
    lmin = 1 
    lmax = 0.521*np.array(strainmax)
    #lmax = 0.6*np.array(strainmax)
    if any(ele < 1 for ele in lmax):
        lmax = 1.2*np.ones(reps)
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
    plt.vlines([lmin,lmax[i]],ymin=0,ymax=5,linestyles='dashed')

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

#deff = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/Mechprop/params.dat',dtype=float)[:,2]
#deff = 5.96 # in %
ddplotter = []
ddd = []
gcplotter = []
gcdefff = []
geplotter = []
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
    ddplotter.append(popt[0])
    ddd.append(popt[0])
    gcplotter.append(popt[1])
    gcdefff.append(popt[1]/((1-deff/100)))
    geplotter.append(popt[2])
    gedefff.append(popt[2]/((1-deff/100)))
    RRlist.append(RR)
plt.legend()
plt.xlim(left=1)
plt.ylim(bottom=0,top=max(maxstresses)+1)

plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Fitting data: Gc = '+str(np.round(np.average(gcdefff),2))+' Ge = '+str(np.round(np.average(gedefff),2))+' '+r'$\delta^2 = $'+str(np.round(np.average(ddd),2)))
plt.show()
def Ac(Nc,Ne):
    a = (3*Ne)/(2*Nc)
    return 0.5 + (np.sqrt(a/np.pi)*np.exp(-a))/erf(np.sqrt(a))
rho = 0.92 #g/cm3
T   = 298 #K
R   = 8.3144 #J/(mol K)
Ms  = 131 #g/mol
Na  = 6.022*np.float_power(10,23) #Avogadro constant

Ne  = []
Nc  = []
nue = []
nuc = []

z = 0.387 #From revisiting segmental order...

#Dres = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/Mechprop/params.dat',dtype=float)[:,0]
#Dres = 150



for i in range(reps):
    Nee  = (1/np.sqrt(6))*rho*R*T/(gedefff[i]*Ms)
    def func(Nc):
        return Nc - rho*R*T*Ac(Nc,Nee)/(gcdefff[i]*Ms)
    nuee = Na*rho/(Nee*Ms)
    Ne.append(Nee)
    nue.append(nuee)

    Ncc  = fsolve(func,x0=1) #Ncc is an array, so we select the component zero which is a float
    Nc.append(Ncc[0])
    nucc = Na*rho/(Ncc[0]*Ms)
    nuc.append(nucc)

Ne     = np.asarray(Ne)
Nc     = np.asarray(Nc)
nue    = np.asarray(nue)
nuc    = np.asarray(nuc)
n      = Nc*(1-Ac(Nc,Ne))/(2*Ac(Nc,Ne))
m      = (2/Nc)*z*np.sqrt((Nc+2*n)/(0.8*Ne))*np.arctan(Nc/(2*np.sqrt(Nc*n+n**2)))

Dstatk = 1.213*Dres/(m)
Dstatk = Dstatk/1.4

gcdefferr  = t.ppf(0.95,df=len(gcdefff)-1)*np.std(gcdefff)/np.sqrt(len(gcdefff))
gedefferr  = t.ppf(0.95,df=len(gedefff)-1)*np.std(gedefff)/np.sqrt(len(gedefff))
dderr      = t.ppf(0.95,df=len(ddd)-1)*np.std(ddd)/np.sqrt(len(ddd))


Ncerr     = t.ppf(0.95,df=len(Nc)-1)*np.std(Nc)/np.sqrt(len(Nc))
Neerr     = t.ppf(0.95,df=len(Ne)-1)*np.std(Ne)/np.sqrt(len(Ne))
nucerr    = t.ppf(0.95,df=len(nuc)-1)*np.std(nuc)/np.sqrt(len(nuc))
nueerr    = t.ppf(0.95,df=len(nue)-1)*np.std(nue)/np.sqrt(len(nue))
nerr      = t.ppf(0.95,df=len(n)-1)*np.std(n)/np.sqrt(len(n))
zerr      = 0
Acerr     = t.ppf(0.95,df=len(Ac(Nc,Ne)-1)*np.std(Ac(Nc,Ne))/np.sqrt(len(Ac(Nc,Ne))))
Dstatkerr = t.ppf(0.95,df=len(Dstatk)-1)*np.std(Dstatk)/np.sqrt(len(Dstatk))

gcdeff = np.average(gcdefff)
gedeff = np.average(gedefff)
dd     = np.average(ddd)

Acav = np.average(Ac(Nc,Ne))
Nc   = np.average(Nc)
Ne   = np.average(Ne)
nuc  = np.average(nuc)
nue  = np.average(nue)
n    = np.average(n)
z    = np.average(z)

m    = np.average(m)

Dstatk = np.average(Dstatk)

fig6 = plt.figure()
for i in range(reps):
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1)) 
    plt.plot(strain["strain{0}".format(i)],stressfit(l=strain["strain{0}".format(i)],dd=ddplotter[i],gc=gcplotter[i],ge=geplotter[i]),label='HSH fitting for sample {0} with R squared = {1}'.format(i+1,np.round(RRlist[i],5)),linestyle='dashed')
    plt.vlines([lmin,lmax[i]],ymin=0,ymax=5,linestyles='dashed')
    plt.xticks([lmin,lmax[i]],[lmin,np.round(lmax[i],1)])
plt.legend()
plt.xlim(left=1)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Fitting data: Gc = '+str(np.round(np.average(gcdeff),2))+' Ge = '+str(np.round(np.average(gedeff),2))+' '+r'$\delta^2 = $'+str(np.round(np.average(dd),2)))
plt.ylim(bottom=0,top=max(maxstresses)+1)
plt.text(0.65*plt.xlim()[1],0.7*plt.ylim()[1],r"$D_{stat}/k = %i Hz$" % Dstatk,{'fontsize': 15})
plt.show()

print('')
if True:
    print('Gc WITH DEFECTS')
    print(gcdeff,' +- ',gcdefferr)
    print('Ge WITH DEFECTS')
    print(gedeff,' +- ',gedefferr)
    print('delta2')
    print(dd,' +- ',dderr)
    print('nuc')
    print(nuc,' +- ',nucerr)
    print('nue')
    print(nue,' +- ',nueerr)
    print('Nc')
    print(Nc,' +- ',Ncerr)
    print('Ne')
    print(Ne,' +- ',Neerr)
    print('z')
    print(z,' +- ',zerr)
    print('Ac')
    print(Acav,' +- ',Acerr)
    print('Dstatk')
    print(Dstatk,' +- ',Dstatkerr)
    print(' ')
    print(np.average(RRlist))
saving = bool(input('Hit enter if you do not want to save the data. '))
if saving:
    dirdest  = path+'Analysis-'+str(dirID)
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
    fig6.savefig(dirdest+'/Figure6.pdf')

    info = np.asarray([[gcdeff,gcdefferr,gedeff,gedefferr,dd,dderr,Nc,Ncerr,Ne,Neerr,Acav,Dstatk,m]])

    np.savetxt(dirdest+'/data'+str(dirID)+'.dat',info,delimiter=',')
else:
    print('Retry.')
