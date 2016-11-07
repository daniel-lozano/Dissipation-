import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import sys



'''
%%%%%%%%%%%%%%%%%%%% Funciones utilizadas en el calculo del tiempo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

def integrand(x,k,eta,hbar,m,a):
    
    factor1= - 2*k*a*x

    factor2= - (eta/hbar)*pow(x*a,2)

    exp= np.exp(factor1+factor2)

    div= hbar*k + eta*x*a

    return a*m*exp/div


def integrand1(xt,x,k,eta,hbar,m,a):
    
    factor1= - 2*k*x
    
    factor2= - (eta/hbar)*pow(x,2)
    
    cte=np.exp(+ 2*k*xt+ (eta/hbar)*pow(xt,2))
    
    exp= np.exp(factor1+factor2)
    
    div= hbar*k + eta*x
    
    return cte*m*exp/div



'''
%%%%%%%%%%%%%%%%%%%% Variables usadas en el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


#definiendo el espacio


print("el argumneto es",sys.argv[0])

#definiendo energias adimensionales
Ne=100
Eb=np.linspace(0.01,0.99,Ne)
Ea=np.linspace(1.01,2.0,Ne)

#Constantes usadas

a=20.7989E5  #fermi
Vo=1.8E-6 # MeV
m=0.51     # MeV
hbar=197.327# MeV*fermi =[hbar*c]

print("\nmass=",m, " MeV")
print("L=",a," fermi")
print("Vo=",Vo, " MeV")
print("hbar=",hbar," Mev*fermi\n")


'''
Tenemos que tener que eta*l << hbar*k, para satisfacerlo tomaremos el k mas pequeño posible correspondiente a E=0
'''


#eta=np.linspace(0,5E-11,2) # eV/Fermi
eta=[0,float(sys.argv[1])]
print("eta=",eta)


Tb=np.zeros(len(Eb))
Ta=np.zeros(len(Eb))

'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comenzando el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


##solving for E bellow

for j in range(len(eta)):

    for i in range(len(Eb)):
    
        k=np.sqrt(2*m*Vo*(1-Eb[i]))/hbar
        
        kp=np.sqrt(2*m*Vo*(Eb[i]))/hbar
        
        N_factor= 2*np.sqrt(k*kp/(k**2+kp**2))
    
        funcion= lambda x: integrand(x,k,eta[j],hbar,m,a)
    
        resultados=integrate.quad(funcion,0,1,limit=100,limlst=96)

        Tb[i]=(N_factor**2)*(resultados[0])* (6.58*10**(-22)/hbar)
        
        nombre="$ \eta= $"+str(eta[j]) +"$  MeV $"

        
    plt.plot(Eb,Tb,label=nombre)

##imrpimiento archivo

f=open("bellow_wkb.dat","w")

for i in range(len(Tb)):
    
    f.write(str(Eb[i])+str( )+str(Tb[i]))

f.close()


dE=np.zeros(len(Ea))


for i in range(len(dE)):
    dE[i]=eta[1]*np.sqrt(2*(Ea[i]-1)/(m*Vo))*a


print("minimun dissipation ",min(dE))

print("maximun dissipation ",max(dE))



##solving for E above


for j in range(len(eta)):
    
    for i in range(len(Ea)):
        
        k=np.sqrt(2*m*Vo*(Ea[i]-1))/hbar
        
        kp=np.sqrt(2*m*Vo*(Ea[i]))/hbar
        
        N_factor= 2*np.sqrt(k*kp/(k**2+kp**2))
        
        xturning=0#np.sqrt((Ea[i]-Vo)*m/eta[j])
        
        if(Ea[i]<1+max(dE) and eta[j]!=0):
            xturning=np.sqrt((Ea[i]-Vo)*m/eta[j])
        
        
        funcion1= lambda x: integrand1(xturning,x,k,eta[j],hbar,m,a)
        
        resultados=integrate.quad(funcion1,xturning,a,limit=100,limlst=96)
        
        Ta[i]=(N_factor**2)*(resultados[0])* (6.58*10**(-22)/hbar)
    
        nombre="$ \eta= $"+str(eta[j]) +"$  MeV $"
        
    
    plt.plot(Ea,Ta,label=nombre)









#plt.figure(figsize=(20,10))

plt.legend(loc=2)
plt.xlabel("$  E/V_0 $",size=20)
#plt.ylim(0,1E-15)
plt.ylabel(" $ Time $ ",size=20)
plt.savefig("Integral_disp_vel.png")

plt.show()









