import numpy as np
import pandas as pd
import math
import h5py
import scipy as sp
import matplotlib.pyplot as plt

pval = 4 #desired ln p value
nfiles = 501 #number of files
step = 100 #timestepsize
xch = 1 #choosen x value
Bch = 1 #choosen B value
uch=0.03 
r=4
nu=1
pinj=1

c = 299792458
q0 = 1.602176634e-19
B0 = 1e-10
m0 = 1.672621923e-27


imp = pd.read_csv('extract/data.csv',skiprows=0, usecols=(1,2), delimiter=',')
imp = imp.to_numpy ()
t = imp [:,0]
file = imp [:,1]
f =[]
p = []
for n in range (0,nfiles):
    y = pd.read_csv('extract/'+file[n],skiprows=0, usecols=(0,8), delimiter=',')
    y = y.to_numpy ()
    ff = y [:,0]
    pp = y [:,1]
    f.append(ff)
    p.append(pp)
pdiff = []
for i in range (0,nfiles):
    for n in range (0,len(p[i])):
        pdif = p[i]-pval
        if (pdif[n]>0):
            pdiff.append(i)
            pdiff.append(n)
            pdiff.append(p[i][n])
            break
pdlen = len(pdiff)/3
pdlen=int(pdlen)
t=[]
n=pdiff[1]
p=np.exp(pdiff[2])
for i in range (0,pdlen):
    t.append(int(pdiff[3*i])*step)
fr=[]
#pr=[]
for i in range (0,nfiles):
    fr.append(f[i][n])
    #pr.append(p[i][n])
f=fr
#p=np.exp(pr[0])

'''eponenttial p
tanh in scat freq
power law'''


#theory
rg = m0*c/(q0*B0)
x = xch*rg
u = uch * c
f500=f[-1]
tau = 4/(3*nu*uch*uch)
fcalc=[0]
print (x,u)
c1=0.25*x/u+3*(r+1)/(4*(r-1))*np.log(p/pinj)*tau
c2=tau*c1/2
for i in range (1,nfiles):
    phi=0.5*(np.exp(2*c1*c1/c2)*sp.special.erfc(np.sqrt(c1*c1*c1/(2*t[i]*c2))+np.sqrt(c1*t[i]/(2*c2)))+sp.special.erfc(np.sqrt(c1*c1*c1/(2*t[i]*c2))-np.sqrt(c1*t[i]/(2*c2))))
    fi=f500*phi
    fcalc.append(fi)


injectionRate = r
injectionMomentum =pinj
UShock =u
compressionRatio = r
usedPValue = p

f000SteadyStateAna = 3 * injectionRate / (
    injectionMomentum * UShock) * compressionRatio / (compressionRatio - 1) * (
        usedPValue / injectionMomentum)**(-3 * compressionRatio /
                                          (compressionRatio - 1)) #erstes p0 zu viel? -> drury 992

# Use the final value of the simulation as the steady state value
f000SteadyState = f500

# Compute the mean acceleration time
scatteringFrequency = 1
tAcc = compressionRatio / (2 * UShock**2 * scatteringFrequency) * (
    compressionRatio + 1) / (compressionRatio - 1) * np.log(
        (1 + usedPValue**2) / (1 + injectionMomentum**2)) #why

# Analytical solution Drury
# gamma = np.sqrt(usedPValue**2 + 1)
# kappa =
# c1 = 3 / 4 * (compressionRatio - 1) / (compressionRatio + 1) * np.log(
#     usedPValue / injectionMomentum) * tau
c1 = tAcc
c2 = 1 / 3 * compressionRatio / UShock**4 * (compressionRatio**3 + 1) / (
    compressionRatio -
    1) * 1 / scatteringFrequency**2 * (1 / (usedPValue**2 + 1) - 1 /
                                       (injectionMomentum**2 + 1) + np.log(
                                           (usedPValue**2 + 1) /
                                           (injectionMomentum**2 + 1)))
anaSolDruryy = [0]
for i in range (1,nfiles):
    anaSolDrury = f000SteadyStateAna / 2 * (
    np.exp(2 * c1**2 / c2) * sp.special.erfc(
        np.sqrt(c1**3 / (2 * t[i] * c2)) + np.sqrt(c1 * t[i] /
                                                       (2 * c2))) +
    sp.special.erfc(
        np.sqrt(c1**3 / (2 * t[i] * c2)) - np.sqrt(c1 * t[i] /
                                                       (2 * c2))))
    anaSolDruryy.append(anaSolDrury)

anaSolDrury=anaSolDruryy
plt.plot (t,anaSolDrury,label='nils')


plt.title ('x='+str(xch)+' , p='+str(p))
plt.plot (t,f,label='sapphire')
plt.plot (t, fcalc,label='theory')
plt.legend ()
plt.show ()
