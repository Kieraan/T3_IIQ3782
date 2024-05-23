from epcsaftpy import component, pcsaft
import numpy as np
from sgtpy.equilibrium import bubblePy,lle,lle_init
import matplotlib.pyplot as plt
Water = component('Water', ms = 1.2046817736, sigma = [2.7927, 10.11, -0.01775, -1.417, -0.01146], eps = 353.9449,
                 kappaAB = 0.045090, eAB = 2425.6714, sites = [0, 1, 1], Mw = 18.015)

Butanol = component('Butanol', ms = 4.2102, sigma = 3.0741 , eps = 219.92, 
               kappaAB = 0.006692, eAB = 1890.72, sites = [0, 1, 1], Mw = 74.123)

Ethanol = component('Ethanol', ms = 2.3827, sigma = 3.1771 , eps = 198.24, 
               kappaAB = 0.032384, eAB = 2653.4, sites = [0, 1, 1], Mw = 46.069)
mix1 = Water + Butanol
mix2 = Water + Ethanol
# adding the binary parameters
mix1.set_kijsaft(i = 0, j = 1, kij0 = -0.020102)
mix2.set_kijsaft(i = 0, j = 1, kij0 = -0.045)
eos1 = pcsaft(mix1)
eos2 = pcsaft(mix2)

# computing lle
T = 298.15 # K
P = 1E5 # Pa
z = np.array([0.7, 0.3])
x0, w0 = lle_init(z, T, P,eos1)
sol = lle(x0,w0,z,T,P,eos1)
xL11,xL21,_=sol
x0, w0 = lle_init(z, T, P,eos2)
sol = lle(x0,w0,z,T,P,eos2)
xL12,xL22,_=sol

#Stability function
def Stability (T,P,z,eos,rho0=None,Xass0=None):
    temp_aux  = eos.temperature_aux(T)
    nc = eos.nc
    rho, Xass0 = eos.density_aux(z, temp_aux, P, "L", rho0=rho0, Xass0=Xass0)
    rhoi = z*rho
    # muij
    muij_sol = eos.dmuad_aux(rhoi,temp_aux,Xass0)
    muij  = muij_sol[1]
    Stab_mat=np.zeros_like(muij)
    for i in range(nc):
        for j in range(nc):
            Stab_mat[i,j]=rho*muij[i,j]
    return np.linalg.det(Stab_mat)

x1=np.linspace(0.001,0.999,100)
x2=1-x1
z1=np.stack((x1,x2))
Stabi1=np.asarray([Stability(T,P,z1[:,i],eos1) for i,val in enumerate(z1[0,:])])
Stabi2=np.asarray([Stability(T,P,z1[:,i],eos2) for i,val in enumerate(z1[0,:])])
fig,ax=plt.subplots(2,1)
ax[0].plot(x1,Stabi1,"r-")
ax[0].plot(x1[Stabi1>0],Stabi1[Stabi1>0],'g--')
ax[0].plot(xL11[0],0,"kx")
ax[0].plot(xL21[0],0,"kx")
ax[0].set_ylim([-400,800])

ax[1].plot(x1,Stabi2,"g-")
ax[1].set_ylim([-400,800])
plt.show()

# plt.show()
