# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 12:32:45 2021

@author: icrma
"""

# %% strain driven shear coupling

from numpy.linalg import inv
import numpy as np
import matplotlib.pyplot as plt
import math as m

phi0 = .3
phi1 = 10
C11,C12,C44 = 169.3097, 87.2201, 41.0448 # moduli for copper
Cinv11,Cinv12,Cinv44 = 1/110, -0.34/110, 2*(1+0.34)/110 
P0 = np.zeros((3,3)) # applied load
Fgb12 = np.zeros((3,3)) # grain boundary shear
Fgb1 = np.identity((3)) # eigenstrain for grain 1
Fgb2 = np.identity((3)) # eigenstrain for grain 2
F0 = np.identity((3))
F0[1,1] = 0 # strain
Fgb12[0,1] = 1
a12 = 1
th_g = np.radians(30/2)
# yield at 1/8*phi0*fgb12^[0^(-1)]Cinv 

tfin = 5
dt = 0.0001
tt = int(tfin/dt)
V1 = np.zeros(tt+1)
V2 = np.zeros(tt+1)
h12 = np.zeros(tt+1)
time = np.zeros(tt+1)
strain = np.zeros(tt+1)
stress = np.zeros(tt+1)
gib1 = np.zeros(tt+1)

Vtot = 1.0
V1[0] = .3
V2[0] = Vtot - V1[0] 
v1dot= 0
v2dot = 0
dh12 = 0
dh21 = 0

C = np.array([[C11, C12, C12, 0, 0, 0],
              [C12, C11, C12, 0, 0, 0],
              [C12, C12, C11, 0, 0, 0],
              [0, 0, 0, C44, 0, 0],
              [0, 0, 0, 0, C44, 0],
              [0, 0, 0, 0, 0, C44]])

def CalSig(F,mat):
    # calculates sigma from a given strain
    
    sig = np.zeros((3,3))
    sig[0,0] = mat[0,0]*F[0,0] + mat[0,1]*F[1,1]  + mat[0,2]*F[2,2] + 2*(mat[0,3]*F[1,2] + mat[0,4]*F[2,0] + mat[0,5]*F[0,1])
    sig[1,1] = mat[1,0]*F[0,0] + mat[1,1]*F[1,1]  + mat[1,2]*F[2,2] + 2*(mat[1,3]*F[1,2] + mat[1,4]*F[2,0] + mat[1,5]*F[0,1])
    sig[2,2] = mat[2,0]*F[0,0] + mat[2,1]*F[1,1]  + mat[2,2]*F[2,2] + 2*(mat[2,3]*F[1,2] + mat[2,4]*F[2,0] + mat[2,5]*F[0,1])
    sig[1,2] = mat[3,0]*F[0,0] + mat[3,1]*F[1,1]  + mat[3,2]*F[2,2] + 2*(mat[3,3]*F[1,2] + mat[3,4]*F[2,0] + mat[3,5]*F[0,1])
    sig[0,2] = mat[4,0]*F[0,0] + mat[4,1]*F[1,1]  + mat[4,2]*F[2,2] + 2*(mat[4,3]*F[1,2] + mat[4,4]*F[2,0] + mat[4,5]*F[0,1])
    sig[0,1] = mat[5,0]*F[0,0] + mat[5,1]*F[1,1]  + mat[5,2]*F[2,2] + 2*(mat[5,3]*F[1,2] + mat[5,4]*F[2,0] + mat[5,5]*F[0,1])
    sig[1,0] = sig[0,1]
    sig[2,1] = sig[1,2]
    sig[0,2] = sig[2,0]
    return sig

def CalAstar(grain,v1,v2):
    Fstar = np.zeros((3,3))
    A = 0
    if(grain == 1):
        
        Fstar = (Vtot*F0 + v2*(Fgb1 - Fgb2))/Vtot
        F = Fstar - Fgb1
        P = CalSig(F,C)  
        for i in range(0,2):
            for j in range(0,2):
                A = A + 0.5*F[i,j]*P[i,j]
                
    elif(grain == 2):
        
        Fstar = (Vtot*F0 + v1*(Fgb2 - Fgb1))/Vtot
        F = Fstar - Fgb2
        P = CalSig(F,C)
        for i in range(0,2):
            for j in range(0,2):
                A = A + 0.5*F[i,j]*P[i,j] 
    return A,Fstar

def Caldadh(v1,v2):
    dh = 0
    
    P = CalSig(F0,C)
    for i in range(0,2):
            for j in range(0,2):
                dh = dh + (v1-v2)/Vtot*(2*Fgb12[i,j] + Fgb1[i,j] - Fgb2[i,j])*P[i,j] 
    return dh
# set up rotation matrix
theta = np.radians(0)
c, s = np.cos(theta), np.sin(theta)
rot2 = np.array(((c,0, s),(0,1,0), (-s,0, c)))
Fgb12 = np.matmul(np.matmul(rot2,Fgb12),np.transpose(rot2))

for t in range(tt):
    
    F0[0,1] = .17*np.sin(20*m.pi/(tt)*t) # applied strain
    
    A1,F1 = CalAstar(1,V1[t],V2[t]) # A* for grain 1
    A2,F2 = CalAstar(2,V1[t],V2[t]) # A* for grain 2
    dadh = Caldadh(V1[t],V2[t])
    
    dh12 = -(A1 - A2 + dadh + phi0)/phi1
    
    dh21 = -(-A1 + A2 - dadh + phi0)/phi1
    gib1[t] = A1 - A2 
    # stop if the gb interface moves more than 1/2 so that the volume does not become negative
    if (h12[t] < -1/2 or h12[t] > 1/2):
        dh12 = 0
        dh21 = 0
    
    if(dh12 < 0):
        dh12 = 0
    if(dh21 < 0):
        dh21 = 0
    # calculate eigenstrain 
    if(dh12 > 0):
        # V1 gets bigger
        Fgb1 = np.add(Fgb1, a12*dh12/V1[t]*dt*Fgb12)
        if(Fgb2[0,1] != 0):
            Fgb2 = np.add(Fgb2, -a12*dh12/V2[t]*dt*Fgb12)
        
    if(dh21 > 0):
        # V1 gets smaller 
        dh12 = -dh21 # set hpq = hqp via antisymmetry 
        Fgb2 = np.add(Fgb2, -a12*dh12/V2[t]*dt*Fgb12)
        if(Fgb1[0,1] != 0):
            Fgb1 = np.add(Fgb1, a12*dh12/V1[t]*dt*Fgb12)
    
    h12[t+1] = h12[t] + dt*dh12 # euler integration of interface position and volume
    v1dot = a12*dh12 
    v2dot = -a12*dh12
    
    V1[t+1] = V1[t] + v1dot*dt 
    V2[t+1] = V2[t] + v2dot*dt 
    
    time[t+1] = t*dt
    
    # stress strain curve
    temp = F1 - Fgb1
    P1 = CalSig(temp,C)
    temp = F2 - Fgb2
    P2 = CalSig(temp,C)
    
    stress[t+1] = (V1[t]*P1[0,1] + V2[t]*P2[0,1])/Vtot
    strain[t+1] = F0[0,1] 

fig = plt.figure()
ax = fig.add_subplot()
plt.xlabel('Time')
plt.ylabel('Change in h')
plt.title('Change in $h_{12}$')
plt.grid(alpha=.7,linestyle='-.')
plt.plot(time,h12)

plt.show()

fig1, ax1 = plt.subplots()
ax1.plot(time, V1)
ax1.plot(time, V2)
ax1.set_xlabel("Time")
ax1.set_ylabel("Volume")
ax1.set_title('Volume Change')
ax1.legend(['V1','V2'])
plt.grid(alpha=.7,linestyle='-.')

fig2, ax2 = plt.subplots()
ax2.plot(strain, stress)
ax2.set_xlabel("Strain")
ax2.set_ylabel("Stress")
ax2.set_title('Stress strain curve')
plt.grid(alpha=.7,linestyle='-.')
#plt.xlim([-1.5, 1.5])

fig1, ax1 = plt.subplots()
ax1.plot(time, gib1)
ax1.set_xlabel("Time")
ax1.set_ylabel("Forces")
ax1.set_title('Thermodynamic Forces')
plt.grid(alpha=.7,linestyle='-.')

