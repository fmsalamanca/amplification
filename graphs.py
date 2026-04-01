import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

dirID       = 14
path     = '/home/fernando/Papers/Fillers/Mechprop/CICLOS/'+'curvas_exp/'


def stressfit_C(l,dd,gc,ge,C): #the most useful form. Note that here dd = d**2. C gives flexibility for initial stress after corrections in strain axis
    denom  = 1 - (dd)*((l**2) + 2/l - 3)
    factor = l - 1/(l**2)
    term1 = (1 - dd)*factor/denom
    term2 = dd*factor/(denom**2)
    term3 = -l**float(-1 - 1) + l**float(1/2 - 1)
    return gc*(term1 - term2) + 2*(ge/1)*term3 + C

df = pd.read_csv(path+"curva_exp"+str(dirID)+'.csv',header=None,dtype=np.float64) 

strain = df.values[:,1]
stress = df.values[:,2]

strain0 = strain[strain<2]
strain1 = strain[strain<2.5]
strain2 = strain[strain<3]
strain3 = strain[strain<3.5]
strain4 = strain[strain<4]
strain5 = strain[strain<4.5]

#A88631. N990 50 phr. 2900 points removed
gc = [0.32,0.35,0.39,0.35,0.27]
ge = [0.39,0.30,0.26,0.34,0.67]
dd = [0.08,0.05,0.04,0.05,0.06]
C  = [0.06,0.08,0.10,0.05,-0.09]
RR = [0.9832,0.9948,0.9966,0.9950,0.9917]

if True:
    #A73012. N330 50 phr. 2475 points removed
    gc = [0.46,0.50,0.44,None,1.03]
    ge = [0.47,0.35,0.48,None,0.00]
    dd = [0.1,0.07,0.09,None,0.06]
    C  = [0.15,0.18,0.15,None,-0.65]
    RR = [0.9829,0.9939,0.9816,None,0.8792]

if True:
    #A66455 N234 50 phr. 1675 points removed
    gc = [0.53,0.63,0.68,0.70,0.77]
    ge = [0.44,0.27,0.18,0.12,0.00]
    dd = [0.10,0.10,0.10,0.10,0.10]
    C  = [0.23,0.26,0.29,0.31,0.32]
    RR = [0.8902,0.9700,0.9850,0.9910,0.9532]

if True:
    #A93501 N330 40 phr. 2850 points removed
    gc = [0.36,0.40,0.44,0.33,None,0.76]
    ge = [0.43,0.31,0.22,0.56,None,0.00]
    dd = [0.10,0.06,0.05,0.08,None,0.05]
    C  = [0.00,0.11,0.14,0.02,None,-0.39]
    RR = [0.9936,0.9964,0.9973,0.9863,None,0.8646]

fig = plt.figure()
plt.plot(strain,stress,label='Experimental curve',color='blue')
plt.plot(strain5,stressfit_C(l=strain5,dd=dd[5],gc=gc[5],ge=ge[5],C=C[5]),label='Fitting max 4.5, R^2 = {0}'.format(RR[5]),color='black')
#plt.plot(strain5,stressfit_C(l=strain4,dd=dd[4],gc=gc[4],ge=ge[4],C=C[4]),label='Fitting max 4, R^2 = {0}'.format(RR[4]),color='orange')
plt.plot(strain5,stressfit_C(l=strain5,dd=dd[3],gc=gc[3],ge=ge[3],C=C[3]),label='Fitting max 3.55, R^2 = {0}'.format(RR[3]),color='green')
plt.plot(strain5,stressfit_C(l=strain5,dd=dd[2],gc=gc[2],ge=ge[2],C=C[2]),label='Fitting max 3, R^2 = {0}'.format(RR[2]),color='red')
plt.plot(strain5,stressfit_C(l=strain5,dd=dd[1],gc=gc[1],ge=ge[1],C=C[1]),label='Fitting max 2.5, R^2 = {0}'.format(RR[1]),color='purple')
plt.plot(strain5,stressfit_C(l=strain5,dd=dd[0],gc=gc[0],ge=ge[0],C=C[0]),label='Fitting max 2, R^2 = {0}'.format(RR[0]),color='brown')

plt.legend()
plt.xticks([2,2.5,3,3.5,4,4.5],[2,2.5,3,3.5,4,4.5])
#plt.xticks([1.5,1.75,2,2.25,2.5],[1.5,1.75,2,2.25,2.5])
plt.vlines([2,2.5,3,3.5,4,4.5],ymin=0,ymax=5,linestyles='dashed')
#plt.vlines([1.5,1.75,2,2.25,2.5],ymin=0,ymax=5,linestyles='dashed')
plt.xlim(left=1)
plt.ylim(bottom=0)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Compound A-93501 (N330 with 40 phr)')
plt.show()