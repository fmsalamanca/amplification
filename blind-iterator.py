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

dirID = 'NR-S3RCB'
dat = False
lmin  = 1
lmax1 = 1.5
lmax2 = 2
lmax3 = 3.25

#deff = 4.94
#Dres = 182

path  = '/home/fernando/Nextcloud/Papers/ICB/Mechprop-2023/'+str(dirID)+'PiroliticoPreciclos.is_tens_RawData/'
#path  = '/home/fernando/ownCloud/Papers/ICB/Mechprop/'+str(dirID)+'_ciclos.is_tens_RawData/'
path = 'C:\\Users\\Elastomeros Fernando\\Nextcloud\\Papers\\ICB\\Mechprop-2023\\'+str(dirID)+'Preciclos.is_tens_RawData\\'

#COAN   = 81
if dat:
    phi   = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/phiCB.dat')
    phi   = phi[dirID-1]
else:
    phi   = 0.195
ad    = 1/(1-phi)

auto = True
#geunf = float(input('Ge for the unfilled sample in MPa. Ge = '))
geunf = 0.1

df_list = []
for files in os.listdir(path):
    if files.endswith('.csv'):
        df = pd.read_csv(path+files,header=1,dtype=np.float64,sep=';',encoding='latin-1') 
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
        #precycling = input('Did you see precycling on the curves? Hit Enter if not. ')
        precycling = True
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
        stress["stress{0}".format(i)] = df_list[i].values[2200:len(df_list[i]),3]
        print('Precycling removed from stress ',i+1)
                    #x axis
        strain["strain{0}".format(i)] = df_list[i].values[2200:len(df_list[i]),4]+1
        print('Precycling removed from strain ',i+1)
        

strainselected1 = {}
strainselected2 = {}
strainselected3 = {}
strainaux       = {}
stressselected1 = {}
stressselected2 = {}
stressselected3 = {}
stressaux       = {}
maxstresses     = []
strainmax       = []
#warnings.filterwarnings('ignore',message='invalid value encountered in less')
for i in range(reps):
    strain["strain{0}".format(i)]=strain["strain{0}".format(i)] - strain["strain{0}".format(i)][0]+1#adjust the first point as close as possible to 1
    maxstresses += [max(stress['stress{0}'.format(i)])]
    stress["stress{0}".format(i)]=stress["stress{0}".format(i)] - stress["stress{0}".format(i)][0]#adjust the first point as close as possible to 0
    strainmax += [strain["strain{0}".format(i)][len(strain["strain{0}".format(i)])-1]]
if False: #true to create a csv file with the experimental curve corrected
    pd.DataFrame(np.array([strain["strain0"],stress["stress0"]]).T).to_csv('/home/fernando/ownCloud/Papers/Fillers/Mechprop/CICLOS/curvas_exp/'+"curva_exp"+str(dirID)+".csv",header=None)
    print('File created')

fig3 = plt.figure()
for i in range(reps):
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
plt.legend()
plt.xlim(left=1)
plt.text(0.65*plt.xlim()[1],0.7*plt.ylim()[1],"%i points removed" % points,{'fontsize': 15})
plt.ylim(bottom=0)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Stress-Strain curves corrected')
plt.show()



for i in range(reps):
    strainaux["strain{0}".format(i)]=strain["strain{0}".format(i)][strain["strain{0}".format(i)]<lmax1]
    strainselected1["strain{0}".format(i)]=strainaux["strain{0}".format(i)][strainaux["strain{0}".format(i)]>lmin]
    stressaux["stress{0}".format(i)]=stress["stress{0}".format(i)][strain["strain{0}".format(i)]<lmax1]
    stressselected1["stress{0}".format(i)]=stressaux["stress{0}".format(i)][strainaux["strain{0}".format(i)]>lmin]
    
    strainaux["strain{0}".format(i)]=strain["strain{0}".format(i)][strain["strain{0}".format(i)]<lmax2]
    strainselected2["strain{0}".format(i)]=strainaux["strain{0}".format(i)][strainaux["strain{0}".format(i)]>lmin]
    stressaux["stress{0}".format(i)]=stress["stress{0}".format(i)][strain["strain{0}".format(i)]<lmax2]
    stressselected2["stress{0}".format(i)]=stressaux["stress{0}".format(i)][strainaux["strain{0}".format(i)]>lmin]

    strainaux["strain{0}".format(i)]=strain["strain{0}".format(i)][strain["strain{0}".format(i)]<lmax3]
    strainselected3["strain{0}".format(i)]=strainaux["strain{0}".format(i)][strainaux["strain{0}".format(i)]>lmin]
    stressaux["stress{0}".format(i)]=stress["stress{0}".format(i)][strain["strain{0}".format(i)]<lmax3]
    stressselected3["stress{0}".format(i)]=stressaux["stress{0}".format(i)][strainaux["strain{0}".format(i)]>lmin]


fig4 = plt.figure()
for i in range(reps):
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
plt.legend()
plt.vlines([lmin,lmax1],ymin=0,ymax=5,linestyles='dashed')
plt.xticks([lmin,lmax1,lmax2,lmax3],[lmin,lmax1,lmax2,lmax3])
plt.vlines([lmin,lmax2],ymin=0,ymax=5,linestyles='dashed')
plt.vlines([lmin,lmax3],ymin=0,ymax=5,linestyles='dashed')
plt.xlabel(r'$\lambda$')
plt.xlim(left=1)
plt.ylim(bottom=0)
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Stress-Strain curves fully corrected')
plt.show()

#del strain
#del stress
#del strainaux
#del stressaux
#del df_list
if dat:
    deff = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/Mechprop/params.dat',dtype=float)[:,2]
    deff = deff[dirID-1]
else:
    deff = float(input('1H DQ-NMR test. Fraction of defects (%): '))
#################################################################################################################################################################################################################################

gcinit = []
A = []
B = []

RRlist = []
fig5 = plt.figure()
for i in range(reps):
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
    def func1(l,gc,A,B): #simplest form from Heinrich, Sommer and Helmis model
        epsilon = (l - 1)
        XH      = 5.2*phi**2 + 2.5*phi + 1
        XFF     = 1 + A * (1/(1 + B * np.sign(epsilon)*np.abs(epsilon)**float(1.2)))
        X       = XH*XFF
        return (X/ad)*(gc*(l-l**float(-2))+2*geunf*(-l**float(-2)+l**(-0.5)))
    popt, pcov = curve_fit(func1,xdata=strainselected1["strain{0}".format(i)],ydata=stressselected1["stress{0}".format(i)],bounds=((0,0,0),(np.inf,np.inf,np.inf)))
    RR         = 1-np.sum(np.sqrt(np.diag(pcov)))
    plt.plot(strain["strain{0}".format(i)],func1(l=strain["strain{0}".format(i)],gc=popt[0],A=popt[1],B=popt[2]),label='HSH fitting for sample {0} with R squared = {1}'.format(i+1,np.round(RR,5)),linestyle='dashed')
    gcinit.append(popt[0]/((1-deff/100)))
    A.append(popt[1])
    B.append(popt[2])
    RRlist.append(RR)
plt.legend()
plt.xlim(left=1,right=lmax1+0.1)
plt.ylim(bottom=0,top=2)
plt.vlines([lmin,lmax1],ymin=0,ymax=5,linestyles='dashed')
plt.xticks([lmin,lmax1],[lmin,np.round(lmax1,1)])
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('1st iteration: amplification quantification')
plt.text(1.2*plt.xlim()[0],0.5*plt.ylim()[1],r"$A = %f$" % (np.average(A)),{'fontsize': 15})
plt.text(1.2*plt.xlim()[0],0.4*plt.ylim()[1],r"$B = %f$" % (np.average(B)),{'fontsize': 15})
plt.text(1.2*plt.xlim()[0],0.3*plt.ylim()[1],r"$G_c = %f\,\,MPa}$" % (np.average(gcinit)),{'fontsize': 15})
plt.show()

print('')
print('A')
print(np.average(A))
print('B')
print(np.average(B))
print('gc 1st')
print(np.average(gcinit))
print('')
######################################################################################################################################################################################################################################

gcmeh = []
gemeh = []

RRlist = []
fig6 = plt.figure()
for i in range(reps):
    def func2(l,gc,ge):
        epsilon = (l - 1)
        XH      = 5.2*phi**2 + 2.5*phi + 1
        XFF     = 1 + A[i] * (1/(1 + B[i] * np.sign(epsilon)*np.abs(epsilon)**float(1.2)))
        X       = XH*XFF
        return (X/ad)*(gc*(l-l**float(-2))+2*ge*(-l**float(-2)+l**(-0.5)))
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
    popt, pcov = curve_fit(func2,xdata=strainselected2["strain{0}".format(i)],ydata=stressselected2["stress{0}".format(i)],bounds=((0,0),(np.inf,np.inf)))
    RR         = 1-np.sum(np.sqrt(np.diag(pcov)))
    plt.plot(strain["strain{0}".format(i)],func2(l=strain["strain{0}".format(i)],gc=popt[0],ge=popt[1]),label='HSH fitting for sample {0} with R squared = {1}'.format(i+1,np.round(RR,5)),linestyle='dashed')
    
    gcmeh.append(popt[0])
    gemeh.append(popt[1])
    RRlist.append(RR)
plt.legend()
plt.xlim(left=1,right=lmax2+0.1)
plt.ylim(bottom=0,top=2)
plt.vlines([lmin,lmax2],ymin=0,ymax=5,linestyles='dashed')
plt.xticks([lmin,lmax2],[lmin,np.round(lmax2,1)])
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('2nd iteration: Gc, Ge with delta=0')
plt.text(1.2*plt.xlim()[0],0.5*plt.ylim()[1],r"$G_c = %f\,\,MPa$" % (np.average(gcmeh)),{'fontsize': 15})
plt.text(1.2*plt.xlim()[0],0.4*plt.ylim()[1],r"$G_e = %f\,\,MPa$" % (np.average(gemeh)),{'fontsize': 15})
plt.show()

print('')
print('gcmeh')
print(np.average(gcmeh))
print('gemeh')
print(np.average(gemeh))
print('')
##########################################################################################################################################################################################################################

dd = []
gcit = []

RRlist = []
fig7 = plt.figure()
for i in range(reps):
    def func3(l,dd,gc): #simplest form from Heinrich, Sommer and Helmis model
        epsilon = (l - 1)
        XH      = 5.2*phi**2 + 2.5*phi + 1
        XFF     = 1 + A[i] * (1/(1 + B[i] * np.sign(epsilon)*np.abs(epsilon)**float(1.2)))
        X       = XH*XFF
        denom  = 1 - (dd)*((l**2) + 2/l - 3)
        factor = l - 1/(l**2)
        term1 = (1 - dd)*factor/(denom**2)
        term2 = dd*factor/denom
        term3 = -l**float(-1 - 1) + l**float(1/2 - 1)
        return (X/ad)*(gc*(term1 - term2) + 2*(gemeh[i]/1)*term3)
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
    popt, pcov = curve_fit(func3,xdata=strainselected3["strain{0}".format(i)],ydata=stressselected3["stress{0}".format(i)],bounds=((0,0),(0.1,np.inf)))
    RR         = 1-np.sum(np.sqrt(np.diag(pcov)))
    plt.plot(strain["strain{0}".format(i)],func3(l=strain["strain{0}".format(i)],dd=popt[0],gc=popt[1]),label='HSH fitting for sample {0} with R squared = {1}'.format(i+1,np.round(RR,5)),linestyle='dashed')
    
    dd.append(popt[0])
    gcit.append(popt[1])
    RRlist.append(RR)
plt.legend()
plt.xlim(left=1,right=lmax3+0.1)
plt.ylim(bottom=0,top=2)
plt.vlines([lmin,lmax3],ymin=0,ymax=5,linestyles='dashed')
plt.xticks([lmin,lmax3],[lmin,np.round(lmax3,1)])
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('3rd iteration: Gc and dd with Ge fixed')
plt.text(1.2*plt.xlim()[0],0.5*plt.ylim()[1],r"$\delta^2 = %f$" % (np.average(dd)),{'fontsize': 15})
plt.text(1.2*plt.xlim()[0],0.4*plt.ylim()[1],r"$G_c = %f\,\,MPa$" % (np.average(gcit)),{'fontsize': 15})
plt.show()

print('')
print('delta2')
print(np.average(dd))
print('gcit')
print(np.average(gcit))
print('')
##################################################################################################################################################################################################################################

gcdeff = []
gedeff = []

RRlist = []
fig8 = plt.figure()
for i in range(reps):
    def func4(l,gc,ge): #simplest form from Heinrich, Sommer and Helmis model
        epsilon = (l - 1)
        XH      = 5.2*phi**2 + 2.5*phi + 1
        XFF     = 1 + A[i] * (1/(1 + B[i] * np.sign(epsilon)*np.abs(epsilon)**float(1.2)))
        X       = XH*XFF
        denom  = 1 - (dd[i])*((l**2) + 2/l - 3)
        factor = l - 1/(l**2)
        term1 = (1 - dd[i])*factor/(denom**2)
        term2 = dd[i]*factor/denom
        term3 = -l**float(-1 - 1) + l**float(1/2 - 1)
        return (X/ad)*(gc*(term1 - term2) + 2*(ge/1)*term3)
    plt.plot(strain["strain{0}".format(i)],stress["stress{0}".format(i)],label='Sample {0}'.format(i+1))
    popt, pcov = curve_fit(func4,xdata=strainselected3["strain{0}".format(i)],ydata=stressselected3["stress{0}".format(i)],bounds=((0,0),(np.inf,np.inf)))
    RR         = 1-np.sum(np.sqrt(np.diag(pcov)))
    plt.plot(strain["strain{0}".format(i)],func4(l=strain["strain{0}".format(i)],gc=popt[0],ge=popt[1]),label='HSH fitting for sample {0} with R squared = {1}'.format(i+1,np.round(RR,5)),linestyle='dashed')

    gcdeff.append(popt[0]/((1-deff/100)))
    gedeff.append(popt[1]/((1-deff/100)))
    RRlist.append(RR)
plt.legend()
plt.xlim(left=1,right=lmax3+0.1)
plt.ylim(bottom=0,top=2)
plt.vlines([lmin,lmax3],ymin=0,ymax=5,linestyles='dashed')
plt.xticks([lmin,lmax3],[lmin,np.round(lmax3,1)])
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('4th iteration: Gc and Ge')
plt.text(1.2*plt.xlim()[0],0.5*plt.ylim()[1],r"$G_c = %f\,\.MPa$" % (np.average(gcdeff)),{'fontsize': 15})
plt.text(1.2*plt.xlim()[0],0.4*plt.ylim()[1],r"$G_e = %f\,\,MPa$" % (np.average(gedeff)),{'fontsize': 15})
plt.show()


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

Ne  = []
Nc  = []
nue = []
nuc = []
gcnet = []
nufr = []

z = 0.387 #From revisiting segmental order...

if dat:
    Dres = np.loadtxt('/home/fernando/ownCloud/Papers/Fillers/Dresmax.dat',dtype=float)
    Dres = Dres[dirID-1]
else:
    Dres = float(input('1H DQ-NMR test. Dres* (Hz): '))

Dstatk = 12319*1.4

for i in range(reps):
    Nee  = (1/np.sqrt(6))*rho*R*T/(gedeff[i]*Ms)

    #def func(Nc):
    #    return (1.213)*(Dres[dirID-1]/Dstatk) - (2*z/Nc)*np.sqrt((Nc+2*n(Nc,Nee))/(0.8*Nee))*np.arctan(Nc/(2*np.sqrt(Nc*n(Nc,Nee)+n(Nc,Nee)**2)))
    def func(Nc):
        return (1.213)*(Dres/Dstatk) - (2*z/Nc)*np.sqrt((Nc+2*n(Nc,Nee))/(0.8*Nee))*np.arctan(Nc/(2*np.sqrt(Nc*n(Nc,Nee)+n(Nc,Nee)**2)))

    nuee = Na*rho/(Nee*Ms)
    Ne.append(Nee)
    nue.append(nuee)

    Ncc    = fsolve(func,x0=1) #Ncc is an array, so we select the component zero which is a float
    Nc.append(Ncc[0])
    nucc   = Na*rho/(Ncc[0]*Ms)
    nuc.append(nucc)
    gcnett = Ac(Ncc[0],Nee)*rho*R*T/(Ncc[0]*Ms)
    gcnet.append(gcnett)
    nufrr  = (gcdeff[i] - gcnett*(1-phi))/(kb*T*phi)
    nufr.append(nufrr)

Ne    = np.asarray(Ne)
Nc    = np.asarray(Nc)
nue   = np.asarray(nue)
nuc   = np.asarray(nuc)
gcnet = np.asarray(gcnet)
nufr  = np.asarray(nufr)

n = Nc*(1-Ac(Nc,Ne))/(2*Ac(Nc,Ne))
m = (2/Nc)*z*np.sqrt((Nc+2*n)/(0.8*Ne))*np.arctan(Nc/(2*np.sqrt(Nc*n+n**2)))


gcdefferr  = t.ppf(0.95,df=len(gcdeff)-1)*np.std(gcdeff)/np.sqrt(len(gcdeff))
gedefferr  = t.ppf(0.95,df=len(gedeff)-1)*np.std(gedeff)/np.sqrt(len(gedeff))
dderr      = t.ppf(0.95,df=len(dd)-1)*np.std(dd)/np.sqrt(len(dd))
Aerr       = t.ppf(0.95,df=len(A)-1)*np.std(A)/np.sqrt(len(A))
Berr       = t.ppf(0.95,df=len(B)-1)*np.std(B)/np.sqrt(len(B))


Ncerr   = t.ppf(0.95,df=len(Nc)-1)*np.std(Nc)/np.sqrt(len(Nc))
Neerr   = t.ppf(0.95,df=len(Ne)-1)*np.std(Ne)/np.sqrt(len(Ne))
nucerr  = t.ppf(0.95,df=len(nuc)-1)*np.std(nuc)/np.sqrt(len(nuc))
nueerr  = t.ppf(0.95,df=len(nue)-1)*np.std(nue)/np.sqrt(len(nue))
nufrerr = t.ppf(0.95,df=len(nufr)-1)*np.std(nufr)/np.sqrt(len(nufr))

Acerr  = t.ppf(0.95,df=len(Ac(Nc,Ne)-1)*np.std(Ac(Nc,Ne))/np.sqrt(len(Ac(Nc,Ne))))

gcdeff = np.average(gcdeff)
gedeff = np.average(gedeff)
dd     = np.average(dd)
A      = np.average(A)
B      = np.average(B)

Acav = np.average(Ac(Nc,Ne))
Nc   = np.average(Nc)
Ne   = np.average(Ne)
nuc  = np.average(nuc)
nue  = np.average(nue)
nufr = np.average(nufr)

print('')
if True:
    print('Gc WITH DEFECTS')
    print(gcdeff,' +- ',gcdefferr)
    print('Ge WITH DEFECTS')
    print(gedeff,' +- ',gedefferr)
    print('delta2')
    print(dd,' +- ',dderr)
    print('A')
    print(A,' +- ',Aerr)
    print('B')
    print(B,' +- ',Berr)
    print('nuc')
    print(nuc,' +- ',nucerr)
    #print('Nc = {̣̣0} +/- {̣̣1}'.format(Nc,Ncerr))
    #time.sleep(1)
    print('nue')
    print(nue,' +- ',nueerr)
    #print('Ne = {̣̣0} +/- {̣̣1}'.format(Ne,Neerr))
    #time.sleep(1)
    print('Nc')
    print(Nc,' +- ',Ncerr)
    #print('nuc = {̣̣0} +/- {̣̣1}'.format(nuc,nucerr))
    #time.sleep(1)
    print('Ne')
    print(Ne,' +- ',Neerr)
    print('nufr')
    print(nufr,' +- ',nufrerr)
    print('Ac')
    print(Acav,' +- ',Acerr)
#print(tabulate([['Gc','Ge','Nc','Ne','nuc','nue','n','z','Sb'],np.round([gc,ge,Nc,Ne,nuc,nue,n,z,Sb],decimals=2),np.round([gcerr,geerr,Ncerr,Neerr,nucerr,nueerr,nerr,zerr,Sberr],decimals=2)],tablefmt="latex"))
    print(' ')
    print(np.average(RRlist))
saving = bool(input('Hit enter if you do not want to save the data. '))
if saving:
    dirdest  = path+'Iterator-BlindMAX'
    if not os.path.exists(dirdest):
        os.makedirs(dirdest)
    elif True:
        raise('You should change the directory name for this analysis')
    fig1.savefig(dirdest+'/Figure1.pdf')
    if precycling:
        fig2.savefig(dirdest+'/Figure2.pdf')
    fig3.savefig(dirdest+'/Figure3.pdf')
    fig4.savefig(dirdest+'/Figure4.pdf')
    fig5.savefig(dirdest+'/Figure5.pdf')
    fig6.savefig(dirdest+'/Figure6.pdf')
    fig7.savefig(dirdest+'/Figure7.pdf')
    fig8.savefig(dirdest+'/Figure8.pdf')

    info = np.asarray([[gcdeff,gcdefferr,gedeff,gedefferr,dd,dderr,A,Aerr,B,Berr,Nc,Ncerr,Ne,Neerr,nufr,nufrerr,Acav]])
    f = open(dirdest+'/datablind.dat','a')
    np.savetxt(f,info,delimiter=',')

    f.close()
    print('Saved!')