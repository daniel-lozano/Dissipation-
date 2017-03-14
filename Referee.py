import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sys



a=20.7989E5  #fermi
Vo=1.8E-6 # MeV
m=0.51     # MeV
hbar=197.327# MeV*fermi =[hbar*c]



E=np.linspace(0.01,0.99,100)
k=np.sqrt(2*m*Vo*(1-E))/hbar
kp=np.sqrt(2*m*Vo*(E))/hbar


def N_factor(k,kp):
    return 2*np.sqrt(k*kp/(k**2+kp**2))

def T_exact(k,E,Vo):
    
    termino= (np.sinh(k*a))**2
    
    return (4*E*(1-E)*Vo**2)/(4*E*(1-E)*Vo**2+termino*Vo**2)


T_wkb=np.exp(-2*k*a)
T_wkb_improved=np.exp(-2*k*a)*N_factor(k,kp)**4
T_Exact=T_exact(k,E,Vo)

plt.loglog(E,T_wkb,"k",label="$ T_{WKB} $",linewidth=1.0)
plt.loglog(E,T_wkb_improved,"b",label="$ T_{WKB_improved} $",linewidth=2.0)
plt.loglog(E,T_Exact,"r--",label="$ T_{Exact} $",linewidth=2.0)
plt.xlabel("$E/V_0$",size=20)
plt.ylabel("$ Transmission\ Coefficient $",size=20)
plt.legend(loc=2)

#plt.xlim(0.0,0.4)
#plt.ylim(0,0.000000001)

plt.savefig("Referee.png")
plt.savefig("Referee.pdf")

plt.show()

