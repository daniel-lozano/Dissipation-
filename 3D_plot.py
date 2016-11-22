import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import integrate
import sys
from mpl_toolkits.mplot3d import Axes3D


'''
%%%%%%%%%%%%%%%%%%%% Funciones utilizadas en el calculo del tiempo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

def integrand(x,k,eta,hbar,m,a):
    
    factor1= - 2*k*a*x
    
    factor2= - (eta/hbar)*pow(x*a,2)
    
    exp= np.exp(factor1+factor2)
    
    div= hbar*k + eta*x*a
    
    return a*m*exp/div


def expo(k,eta,hbar,m,a):
    
    factor1= - 2*k*a
    
    factor2= - (eta/hbar)*pow(a,2)
    
    exp= np.exp(factor1+factor2)
    
    return exp

def integrand1(xt,x,k,eta,hbar,m,a):
    
    factor1= - 2*k*x
    
    factor2= - (eta/hbar)*pow(x,2)
    
    cte=np.exp(+ 2*k*xt+ (eta/hbar)*pow(xt,2))
    
    exp= np.exp(factor1+factor2)
    
    div= hbar*k + eta*x
    
    return cte*m*exp/div


'''
%%%%%%%%%%%%%%%%%%%% Variables que serán utilizadas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''



Vo=1.8E-6 # MeV
m=0.51     # MeV
hbar=197.327# MeV*fermi =[hbar*c]
eta=float(sys.argv[1])



Ne=100
Nl=20

Eb=np.linspace(0.01,0.99,Ne)
Lb=np.linspace(1,200,Nl)*1.0E4




E=np.zeros((Nl,Ne))
L=np.zeros((Nl,Ne))
TimeD=np.zeros((Nl,Ne))
TimeND=np.zeros((Nl,Ne))


for i in range(Nl):
    for j in range(Ne):
        E[i][j]=Eb[j]
        L[i][j]=Lb[i]




##solving for E bellow


for j in range(Nl):
    
    for i in range(Ne):
        
        k=np.sqrt(2*m*Vo*(1-E[j][i]))/hbar
        
        kp=np.sqrt(2*m*Vo*(E[j][i]))/hbar
        
        N_factor= 2*np.sqrt(k*kp/(k**2+kp**2))
        
        funcion= lambda x: integrand(x,k,eta,hbar,m,L[j][i])
        
        resultados=integrate.quad(funcion,0,1,limit=100,limlst=96)
        
        TimeD[j][i]=(resultados[0])* (6.58*10**(-22)/hbar)*(N_factor**2)


for j in range(Nl):
    
    for i in range(Ne):
        
        k=np.sqrt(2*m*Vo*(1-E[j][i]))/hbar
        
        kp=np.sqrt(2*m*Vo*(E[j][i]))/hbar
        
        N_factor= 2*np.sqrt(k*kp/(k**2+kp**2))
        
        
        resultados=m*(1-np.exp(-2*k*L[j][i]))/(2*hbar*k)
        
        TimeND[j][i]=(resultados)* (6.58*10**(-22)/hbar)*(N_factor**2)






fig = plt.figure(1,figsize=(10,5))
ax = fig.add_subplot(121, projection='3d')
ax.plot_surface(E, L, TimeD, rstride=8, cstride=8, alpha=0.3)
#ax.plot_wireframe(E, L, TimeD, rstride=10, cstride=10)
cset = ax.contour(E, L, TimeD, zdir='z', offset=-100, cmap=cm.coolwarm)
cset = ax.contour(E, L, TimeD, zdir='x', offset=-40, cmap=cm.coolwarm)
cset = ax.contour(E, L, TimeD, zdir='y', offset=40, cmap=cm.coolwarm)


ax.set_title("With dissipation")

ax.set_xlabel("$ E/V_0  $")
    
ax.set_ylabel("$ L $")
                      
ax.set_zlabel("$ \tau  $")
                      
                      
                      
ax = fig.add_subplot(122, projection='3d')
#ax.plot_wireframe(E, L, TimeND, rstride=10, cstride=10)
ax.plot_surface(E, L, TimeND, rstride=8, cstride=8, alpha=0.3)
cset = ax.contour(E, L, TimeND, zdir='z', offset=-100, cmap=cm.coolwarm)
cset = ax.contour(E, L, TimeND, zdir='x', offset=-40, cmap=cm.coolwarm)
cset = ax.contour(E, L, TimeND, zdir='y', offset=40, cmap=cm.coolwarm)


ax.set_title("Without dissipation")

ax.set_xlabel("$ E/V_0  $")

ax.set_ylabel("$ L $")

ax.set_zlabel("$ \tau  $")

plt.savefig("3D_plot.png")
plt.show(1)

















