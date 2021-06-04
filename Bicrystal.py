# %% Stress Driven Motion

import numpy as np
import matplotlib.pyplot as plt
import math as m
# initualize values
phi0 = 0.1
phi1 = .02
Cinv11,Cinv12,Cinv44 = 1/110, -0.34/110, 2*(1+0.34)/110 # moduli for copper
P0 = np.zeros((3,3)) # applied load
Fgb12 = np.zeros((3,3)) # grain boundary shear
Fgb1 = np.zeros((3,3)) # eigenstrain for grain 1
Fgb2 = np.zeros((3,3)) # eigenstrain for grain 2
F0 = np.zeros((3,3)) # strain
F1 = np.zeros((3,3)) # strain for grain 1
F2 = np.zeros((3,3)) # strain for grain 2
Fgb12[0,1] = 0.6
a12 = 1

tfin = 1.5
dt = 0.00001
tt = int(tfin/dt)
V1 = np.zeros(tt+1)
V2 = np.zeros(tt+1)
h12 = np.zeros(tt+1)
time = np.zeros(tt+1)
strain = np.zeros(tt+1)
stress = np.zeros(tt+1)
gib1 = np.zeros(tt+1)
gib2 = np.zeros(tt+1)
gib1[0] = 0
gib2[0] = 0
V1[0] = 0.5
V2[0] = 0.5
v1dot= 0
v2dot = 0
dh12 = 0
dh21 = 0

def CalF(sig):
    Ff = np.zeros((3,3))
    Ff[0,0] = Cinv11*sig[0,0] + Cinv12*sig[1,1]  + Cinv12*sig[2,2]
    Ff[1,1] = Cinv12*sig[0,0] + Cinv11*sig[1,1]  + Cinv12*sig[2,2]
    Ff[2,2] = Cinv12*sig[0,0] + Cinv12*sig[1,1]  + Cinv11*sig[2,2]
    Ff[0,1] = 0.5*Cinv44*sig[0,1]
    Ff[1,2] = 0.5*Cinv44*sig[1,2]
    Ff[2,0] = 0.5*Cinv44*sig[2,0]
    Ff[1,0] = Ff[0,1]
    Ff[2,1] = Ff[1,2]
    Ff[0,2] = Ff[2,0]
    return Ff

def CalGibbs(grain):
    G = 0
    F = np.zeros((3,3))
    F = CalF(P0)
    if(grain == 1):
        for i in range(0,2):
            for j in range(0,2):
                G = G + 1/2*F[i,j]*P0[i,j] + P0[i,j]*Fgb1[i,j]
    elif(grain == 2):
        for i in range(0,2):
            for j in range(0,2):
                G = G + 1/2*F[i,j]*P0[i,j] + P0[i,j]*Fgb2[i,j]             
    return G

# set up rotation matrix
theta = np.radians(0)
c, s = np.cos(theta), np.sin(theta)
rot2 = np.array(((c,0, s),(0,1,0), (-s,0, c)))
Fgb12 = np.matmul(np.matmul(rot2,Fgb12),np.transpose(rot2))

for t in range(tt):
    
    P0[0,1] = .17*np.sin(8*m.pi/(tt)*t)
    
    G1 = CalGibbs(1)
    G2 = CalGibbs(2)
    dh12 = -(G1 - G2 - P0[0,1]*Fgb12[0,1] + phi0)/phi1
    dh21 = -(-G1 + G2 + P0[0,1]*Fgb12[0,1] + phi0)/phi1
    
    gib1[t+1] = Fgb1[0,1]
    gib2[t+1] = Fgb2[0,1]
    if (h12[t] < -1/2 or h12[t] > 1/2):
        dh12 = 0
        dh21 = 0

    if(dh12 < 0):
        dh12 = 0
    if(dh21 < 0):
        dh21 = 0
    if(dh12 > 0):
        # V1 gets bigger
        Fgb1 = np.add(Fgb1, a12*dh12/V1[t]*dt*Fgb12)
        if(Fgb2[0,1] != 0):
            Fgb2 = np.add(Fgb2, -a12*dh12/V2[t]*dt*Fgb12)
        
    if(dh21 > 0):
        # V1 gets smaller 
        dh12 = -dh21
        Fgb2 = np.add(Fgb2, -a12*dh12/V2[t]*dt*Fgb12)
        if(Fgb1[0,1] != 0):
            Fgb1 = np.add(Fgb1, a12*dh12/V1[t]*dt*Fgb12)
    
    h12[t+1] = h12[t] + dt*dh12 
    v1dot = a12*dh12
    v2dot = -a12*dh12
    
    V1[t+1] = V1[t] + v1dot*dt 
    V2[t+1] = V2[t] + v2dot*dt 
    
    time[t+1] = t*dt
    F0 = CalF(P0)
    stress[t+1] = P0[0,1]
    strain[t+1] = F0[0,1] + V1[t]*Fgb1[0,1] + V2[t]*Fgb2[0,1]
    
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

fig1, ax1 = plt.subplots()
ax1.plot(time, gib1)
ax1.plot(time, gib2)
ax1.set_xlabel("Time")
ax1.set_ylabel("Fgbp")
ax1.set_title('Eigenstrain Change')
ax1.legend(['Fgb1','Fgb2'])
plt.grid(alpha=.7,linestyle='-.')


# %% Strain Driven Motion
import numpy as np
import matplotlib.pyplot as plt
import math as m
# initualize values
phi0 = 0.3
phi1 = 0.5
C11,C12,C44 = 169.3097, 87.2201, 41.0448 # moduli for copper
Cinv11,Cinv12,Cinv44 = 1/110, -0.34/110, 2*(1+0.34)/110 
P0 = np.zeros((3)) # applied load
Fgb12 = np.identity((3)) # grain boundary shear
Fgb1 = np.identity((3)) # eigenstrain for grain 1
Fgb2 = np.identity((3)) # eigenstrain for grain 2
F0 = np.identity((3))
F0[1,1] = 0 # strain
Fgb12[0,2] = 1
a12 = 1
th_g = np.radians(30/2)
# yield at 1/8*phi0*fgb12^[0^(-1)]Cinv 

tfin = .5
dt = 0.00001
tt = int(tfin/dt)
V1 = np.zeros(tt+1)
V2 = np.zeros(tt+1)
h12 = np.zeros(tt+1)
time = np.zeros(tt+1)
strain = np.zeros(tt+1)
stress = np.zeros(tt+1)
gib1 = np.zeros(tt+1)

Vtot = 1.0
V1[0] = 0.5
V2[0] = Vtot - V1[0] 
v1dot= 0
v2dot = 0
dh12 = 0
dh21 = 0

def CalF(sig,th):
    # calculates strain from a given sigma
    Ff = np.identity((3))
    Ff[0,0] = Cinv11*sig[0,0] + Cinv12*sig[1,1]  + Cinv12*sig[2,2]
    Ff[1,1] = Cinv12*sig[0,0] + Cinv11*sig[1,1]  + Cinv12*sig[2,2]
    Ff[2,2] = Cinv12*sig[0,0] + Cinv12*sig[1,1]  + Cinv11*sig[2,2]
    Ff[0,1] = 0.5*Cinv44*sig[0,1]
    Ff[1,2] = 0.5*Cinv44*sig[1,2]
    Ff[2,0] = 0.5*Cinv44*sig[2,0]
    Ff[1,0] = Ff[0,1]
    Ff[2,1] = Ff[1,2]
    Ff[0,2] = Ff[2,0]
    return Ff
def CalSig(F,th):
    # calculates sigma from a given strain
    R = np.array([[np.cos(th)**2, np.sin(th)**2, 0, 0, 0, -2*np.sin(th)*np.cos(th)],
                  [np.sin(th)**2, np.cos(th)**2 ,0 ,0 ,0 ,2*np.sin(th)*np.cos(th)], 
                  [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, np.cos(th), np.sin(th), 0],
                  [0 ,0 ,0 ,-np.sin(th), np.cos(th) ,0],
                  [np.sin(th)*np.cos(th) ,-np.sin(th)*np.cos(th) ,0 ,0 ,0 ,2*np.cos(th)**2-1]])
    
    C = np.array([[C11, C12, C12, 0, 0, 0],
                  [C12, C11, C12, 0, 0, 0],
                  [C12, C12, C11, 0, 0, 0],
                  [0, 0, 0, C44, 0, 0],
                  [0, 0, 0, 0, C44, 0],
                  [0, 0, 0, 0, 0, C44]])
    CC = np.matmul(np.matmul(np.transpose(R),C),R)
    
    sig = np.zeros((3,3))
    sig[0,0] = CC[0,0]*F[0,0] + CC[0,1]*F[1,1]  + CC[0,2]*F[2,2] + 2.0*(CC[0,3]*F[1,2] + CC[0,4]*F[2,0] + CC[0,5]*F[0,1])
    sig[1,1] = CC[1,0]*F[0,0] + CC[1,1]*F[1,1]  + CC[1,2]*F[2,2] + 2.0*(CC[1,3]*F[1,2] + CC[1,4]*F[2,0] + CC[1,5]*F[0,1])
    sig[2,2] = CC[2,0]*F[0,0] + CC[2,1]*F[1,1]  + CC[2,2]*F[2,2] + 2.0*(CC[2,3]*F[1,2] + CC[2,4]*F[2,0] + CC[2,5]*F[0,1])
    sig[1,2] = CC[3,0]*F[0,0] + CC[3,1]*F[1,1]  + CC[3,2]*F[2,2] + 2.0*(CC[3,3]*F[1,2] + CC[3,4]*F[2,0] + CC[3,5]*F[0,1])
    sig[0,2] = CC[4,0]*F[0,0] + CC[4,1]*F[1,1]  + CC[4,2]*F[2,2] + 2.0*(CC[4,3]*F[1,2] + CC[4,4]*F[2,0] + CC[4,5]*F[0,1])
    sig[0,1] = CC[5,0]*F[0,0] + CC[5,1]*F[1,1]  + CC[5,2]*F[2,2] + 2.0*(CC[5,3]*F[1,2] + CC[5,4]*F[2,0] + CC[5,5]*F[0,1])
    sig[1,0] = sig[0,1]
    sig[2,1] = sig[1,2]
    sig[0,2] = sig[2,0]
    return sig

def CalA(grain, V):
    # calculates A* for a grain
    A = 0
    if(grain == 1):
        F1 = (Vtot*F0 + np.add(Fgb1*V,-Fgb2*V))/Vtot
        P = CalSig(F1,-th_g)
        for i in range(0,2):
            for j in range(0,2):
                A = A + 0.5*F1[i,j]*P[i,j]
    elif(grain == 2):
        F2 = (Vtot*F0 + np.add(Fgb2*V,-Fgb1*V))/Vtot
        P = CalSig(F2,th_g)
        for i in range(0,2):
            for j in range(0,2):
                A = A +  0.5*F2[i,j]*P[i,j]
    return A
    
def Caldadh(V1,V2,up):
    
    ddh = 0
    dadh = 0
    dh = 1e-1
    A1 = CalA(1,V1)
    A2 = CalA(2,V2)
    if(up == True): # boundary moves up / V1 gets smaller
        dV1 = -a12*dh
        dV2 = -dV1
        dFgb1 = a12*Fgb12*dh/V1
        dFgb2 = -a12*Fgb12*dh/V2
        
        F1 = (Vtot*F0 + np.add(dFgb1*dV2,-Fgb2*dV2))/Vtot
        F2 = (Vtot*F0 + np.add(dFgb2*dV1,-dFgb1*dV1))/Vtot
        P1 = CalSig(F1,-th_g)
        P2 = CalSig(F2,th_g)
        
        for i in range(0,2):
            for j in range(0,2):
                ddh = ddh + 0.5*F1[i,j]*P1[i,j]*V1/a12 + 0.5*F2[i,j]*P2[i,j]*V2/a12
    else:
        dV1 = a12*dh
        dV2 = -dV1
        dFgb1 = -a12*Fgb12*dh/V1
        dFgb2 = a12*Fgb12*dh/V2
        
        F1 = (Vtot*F0 + np.add(dFgb1*dV2,-Fgb2*dV2))/Vtot
        F2 = (Vtot*F0 + np.add(dFgb2*dV1,-dFgb1*dV1))/Vtot
        P1 = CalSig(F1,-th_g)
        P2 = CalSig(F2,th_g)
        
        for i in range(0,2):
            for j in range(0,2):
                ddh = ddh + 0.5*F1[i,j]*P1[i,j]*V1/a12 + 0.5*F2[i,j]*P2[i,j]*V2/a12
                
    dadh = (ddh - (A1+A2))/dh
    
    return dadh
# set up rotation matrix
theta = np.radians(0)
c, s = np.cos(theta), np.sin(theta)
rot2 = np.array(((c,0, s),(0,1,0), (-s,0, c)))
Fgb12 = np.matmul(np.matmul(rot2,Fgb12),np.transpose(rot2))

for t in range(tt):
    
    F0[0,1] = .8*np.sin(8*m.pi/(tt)*t) # applied strain
    
    A1 = CalA(1,V2[t]) # A* for grain 1
    A2 = CalA(2,V1[t]) # A* for grain 2
    # dA/dhpq
    dadh12 = Caldadh(V1[t],V2[t],True)
    dh12 = -(A1 - A2  + phi0)/phi1
    
    dadh21 = Caldadh(V1[t],V2[t],False)
    # dA/dhqp
    dh21 = -(-A1 + A2  + phi0)/phi1
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
    temp1 = (F0*Vtot + V2[t]*(Fgb1-Fgb2))/Vtot
    temp2 = (F0*Vtot + V1[t]*(Fgb2-Fgb1))/Vtot
    P1 = CalSig(temp1,-th_g)
    P2 = CalSig(temp2,th_g)
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

# %%  interface rotation - stress

import numpy as np
import matplotlib.pyplot as plt
import math as m

th = [0, 30, 45, 60, 90]

for i in range(5):
    
    phi0 = 0.1
    phi1 = .02
    Cinv11,Cinv12,Cinv44 = 1/110, -0.34/110, 2*(1+0.34)/110 # moduli for copper
    P0 = np.zeros((3,3)) # applied load
    Fgb12 = np.zeros((3,3)) # grain boundary shear
    Fgb1 = np.zeros((3,3)) # eigenstrain for grain 1
    Fgb2 = np.zeros((3,3)) # eigenstrain for grain 2
    F0 = np.zeros((3,3)) # strain
    F1 = np.zeros((3,3)) # strain for grain 1
    F2 = np.zeros((3,3)) # strain for grain 2
    Fgb12[0,1] = 0.6
    a12 = 1
    
    tfin = 1.5
    dt = 0.00001
    tt = int(tfin/dt)
    V1 = np.zeros(tt+1)
    V2 = np.zeros(tt+1)
    h12 = np.zeros(tt+1)
    time = np.zeros(tt+1)
    V1[0] = 0.5
    V2[0] = 0.5
    v1dot= 0
    v2dot = 0
    dh12 = 0
    dh21 = 0
    # set up rotation matrix
    theta = np.radians(th[i])
    c, s = np.cos(theta), np.sin(theta)
    rot2 = np.array(((c,0, s),(0,1,0), (-s,0, c)))
    Fgb12 = np.matmul(np.matmul(rot2,Fgb12),np.transpose(rot2))
    
    for t in range(tt):
        
        P0[0,1] = .25*np.sin(8*m.pi/(tt)*t)
        
        G1 = CalGibbs(1)
        G2 = CalGibbs(2)
        dh12 = -(-G1*np.sign(v1dot) - G2*np.sign(v2dot) - P0[0,1]*Fgb12[0,1] + phi0)/phi1
        dh21 = -(G1*np.sign(v1dot) + G2*np.sign(v2dot) + P0[0,1]*Fgb12[0,1] + phi0)/phi1
        
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
        
    plt.plot(time, h12)
plt.legend([r"$\theta = 0$",r"$\theta = \frac{\pi}{6}$",r"$\theta = \frac{\pi}{4}$",r"$\theta = \frac{\pi}{3}$",r"$\theta = \frac{\pi}{2}$"])
plt.title('Stress driven motion with differing angle')
plt.xlabel('Time')
plt.ylabel('Change in h')
plt.grid(alpha=.7,linestyle='-.')        
plt.show()

# %%  interface rotation - strain

import numpy as np
import matplotlib.pyplot as plt
import math as m

th = [0, 30, 45, 60, 90]

for i in range(5):
    
    # initualize values
    phi0 = 0.3
    phi1 = 0.5
    C11,C12,C44 = 169.3097, 87.2201, 41.0448 # moduli for copper
    Cinv11,Cinv12,Cinv44 = 1/110, -0.34/110, 2*(1+0.34)/110 
    P0 = np.zeros((3,3)) # applied load
    Fgb12 = np.zeros((3,3)) # grain boundary shear
    Fgb1 = np.zeros((3,3)) # eigenstrain for grain 1
    Fgb2 = np.zeros((3,3)) # eigenstrain for grain 2
    F0 = np.zeros((3,3)) # strain
    Fgb12[0,1] = 1
    Fgb12[1,0] = 1
    a12 = 1
    
    # yield at 1/8*phi0*fgb12^[0^(-1)]Cinv 
    
    tfin = .5
    dt = 0.0001
    tt = int(tfin/dt)
    V1 = np.zeros(tt+1)
    V2 = np.zeros(tt+1)
    h12 = np.zeros(tt+1)
    time = np.zeros(tt+1)
    V1[0] = 0.5
    V2[0] = 0.5
    v1dot= 0
    v2dot = 0
    dh12 = 0
    dh21 = 0
    # set up rotation matrix
    theta = np.radians(th[i])
    c, s = np.cos(theta), np.sin(theta)
    rot2 = np.array(((c,0, s),(0,1,0), (-s,0, c)))
    Fgb12 = np.matmul(np.matmul(rot2,Fgb12),np.transpose(rot2))
    
    for t in range(tt):
        
        F0[0,1] = .0008*np.sin(8*m.pi/(tt)*t) # applied strain
        
        A1 = CalA(1,V1[t]) # A* for grain 1
        A2 = CalA(2,V2[t]) # A* for grain 2
        # dA/dhpq
     

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
        
    plt.plot(time, h12)
plt.legend([r"$\theta = 0$",r"$\theta = \frac{\pi}{6}$",r"$\theta = \frac{\pi}{4}$",r"$\theta = \frac{\pi}{3}$",r"$\theta = \frac{\pi}{2}$"])
plt.title('Strain driven motion with differing angle')
plt.xlabel('Time')
plt.ylabel('Change in h')
plt.grid(alpha=.7,linestyle='-.')        
plt.show()
