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
%%%%%%%%%%%%%%%%%%%% Variables usadas en el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''


#definiendo el espacio


print("el argumneto es",sys.argv[0])


#Constantes usadas

a=20.7989E5  #fermi
Vo=1.8E-6 # MeV
m=0.51     # MeV
hbar=197.327# MeV*fermi =[hbar*c]

print("\nmass=",m, " MeV")
print("L=",a," fermi")
print("Vo=",Vo, " MeV")
print("hbar=",hbar," Mev*fermi\n")

#definiendo energias adimensionales

Ne=100
Eb=np.linspace(0.01,0.99,Ne)
Ea=np.linspace(1.01,2.0,Ne)


ENERGY=np.zeros(2*len(Eb))
ETA=np.zeros(2*len(Eb))


for i in range(len(Eb)):
    ENERGY[i]=Eb[i]
    ENERGY[i+len(Eb)]=Ea[i]




# Eta restrictions for the given condition \Delta E<<E

for i in range(len(ENERGY)):
    
    ETA[i]=ENERGY[i]*np.sqrt(m*Vo/(2*abs(1-ENERGY[i])))/a



plt.plot(ENERGY,ETA)
plt.xlabel("$ E/V_0   $",size=20)
plt.title("$ \eta_{max} $",size=20)
plt.savefig("eta_max.png")
plt.show(1)
plt.close()


'''
Tenemos que tener que eta*l << hbar*k, para satisfacerlo tomaremos el k mas pequeño posible correspondiente a E=0
'''



eta1=float(sys.argv[1])
eta=[0,eta1,eta1*2.0]
Min_E=np.zeros(len(eta))

print("eta=",eta)

for i in range(len(eta)):
    Min_E[i]=ENERGY[np.where(eta[i]<ETA)[0][0]]
    print(" eta=",eta[i]," Minimun Energy= ",ENERGY[np.where(eta[i]<ETA)[0][0]])





Tb=np.zeros(len(Eb))
Ta=np.zeros(len(Eb))

'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Comenzando el calculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

#calculando la maxima perdida de energia
dE=np.zeros(len(Ea))



for j in range(len(eta)-1):
    
    for i in range(len(dE)):
        
        dE[i]=(eta[j+1]*np.sqrt(2*(Vo-Vo*Eb[i])/(m))*a)/(Eb[i]*Vo)
        
        ETA_a=(1/a)*(Ea[i])*np.sqrt(m*Vo/2)/np.sqrt(Ea[i]-1)
        
        ETA_b=(1/a)*(Ea[i])*np.sqrt(m*Vo/2)/np.sqrt(Ea[i]-1)
    
    
    plt.plot(Eb,dE,label="$ \eta = $"+str(eta[j+1])+ "$ MeV/fermi $" )


plt.xlabel("$ E/V_0 $",size=20)
plt.ylabel("$  dE/E_b $",size=20)
plt.ylim(0,1)
plt.title("$ Energy\ percentage\ lost\  $")
plt.legend()
plt.savefig("Percent_energy.png")
plt.show(2)
plt.close()





##solving for E bellow

plt.figure(1,figsize=(15,8))

for j in range(len(eta)):

    for i in range(len(Eb)):
    
        k=np.sqrt(2*m*Vo*(1-Eb[i]))/hbar
        
        kp=np.sqrt(2*m*Vo*(Eb[i]))/hbar
        
        N_factor= 2*np.sqrt(k*kp/(k**2+kp**2))
    
        funcion= lambda x: integrand(x,k,eta[j],hbar,m,a)
    
        resultados=integrate.quad(funcion,0,1,limit=100,limlst=96)

        Tb[i]=(N_factor**2)*(resultados[0])* (6.58*10**(-22)/hbar)
        
        nombre="$ \eta= $"+str(eta[j]) +"$  MeV/fermi $"

    plt.subplot(121)
    plt.plot(Eb,Tb,label=nombre)





plt.legend(loc=2)
plt.xlabel("$  E/V_0 $",size=20)
plt.ylabel(" $ T\ [s] $ ",size=20)
plt.xlim(0.1486,1)
plt.title(" $ Dwell\ time\ $",size=20)


#plt.show()
#plt.close()




##imrpimiento archivo

f=open("bellow_wkb.dat","w")

for i in range(len(Tb)):
    
    f.write(str(Eb[i])+str( )+str(Tb[i]))

f.close()





'''

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

'''



#Calculando el tiempo traverso

T_trav=np.zeros(len(Ea))

for j in range(len(eta)):
    
    for i in range(len(Eb)):
        
        k=np.sqrt(2*m*Vo*(1-Eb[i]))/hbar
        
        T_trav[i]=(6.58*10**(-22)/hbar)*m*a/(hbar*k)
        
        if(j!=0):
            A=hbar*k/eta[j]
    
            T_trav[i]= (m/eta[j])*np.log( 1+ a/A)* (6.58*10**(-22)/hbar)
    plt.subplot(122)
    plt.plot(Eb,T_trav,label="$ \eta= $"+str(eta[j]) +"$  MeV/fermi $")



plt.xlim(0.1485,1)
plt.title("$ Traversal\ Time\ $")
plt.xlabel("$ E/V_0 $",size=20)
plt.xlim(0.1486,1)
plt.ylabel("$  Time\ [s] $",size=20)
plt.legend()

plt.savefig("Integral_disp_vel.png")
plt.show()
plt.close()




#Calculando el tiempo de mora de transmision

T_trans=np.zeros(len(Ea))
T_ref=np.zeros(len(Ea))

plt.figure(2,figsize=(15,8))


for j in range(len(eta)):
    
    for i in range(len(Eb)):
        
        k=np.sqrt(2*m*Vo*(1-Eb[i]))/hbar
        
        kp=np.sqrt(2*m*Vo*(Eb[i]))/hbar
        
        N_factor= 2*np.sqrt(k*kp/(k**2+kp**2))
        
        funcion= lambda x: integrand(x,k,eta[j],hbar,m,a)
        
        exponencial=expo(k,eta[j],hbar,m,a)
        
        resultados=integrate.quad(funcion,0,1,limit=100,limlst=96)
        
        T_trans[i]=(N_factor**2)*(resultados[0])* (6.58*10**(-22)/hbar)/(1.0/(1.0+ 1/exponencial))
        
        
        nombre="$ \eta= $"+str(eta[j]) +"$  MeV/fermi $"
    
    plt.subplot(121)
    plt.plot(Eb,T_trans,label=nombre)
    plt.title("$ Transmission\ Dwell\ Time\ $")
    plt.xlabel("$ E/V_0 $",size=20)
    plt.xlim(0.05,1)
    plt.ylabel("$  Time\ [s]\ $",size=20)
    plt.legend()



for j in range(len(eta)):
    
    for i in range(len(Eb)):
        
        k=np.sqrt(2*m*Vo*(1-Eb[i]))/hbar
        
        kp=np.sqrt(2*m*Vo*(Eb[i]))/hbar
        
        N_factor= 2*np.sqrt(k*kp/(k**2+kp**2))
        
        funcion= lambda x: integrand(x,k,eta[j],hbar,m,a)
        
        exponencial=expo(k,eta[j],hbar,m,a)
        
        resultados=integrate.quad(funcion,0,1,limit=100,limlst=96)
        
        T_ref[i]=(N_factor**2)*(resultados[0])* (6.58*10**(-22)/hbar)/((1.0/exponencial)/(1.0+ 1/exponencial))
        
        
        nombre="$ \eta= $"+str(eta[j]) +"$  MeV/fermi $"
    
    plt.subplot(122)
    plt.plot(Eb,T_ref,label=nombre)
    plt.title("$ Reflection\ Dwell\ Time\ $")
    plt.xlabel("$ E/V_0 $",size=20)
    plt.xlim(0.05,1)
    plt.ylabel("$  Time\ [s]\ $",size=20)
    plt.legend()





plt.savefig("trans_ref_time.png",dpi=400)
plt.show()
plt.close()














