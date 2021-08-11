from scipy import linalg
import numpy as np
import math as m
import matplotlib.pyplot as plt
np.set_printoptions(linewidth=2000)

C11,C12,C44 = 169.3097, 87.2201, 41.0448 # moduli for copper
Cinv11,Cinv12,Cinv44 = 1/110, -0.34/110, 2*(1+0.34)/110

N = 2 # Number of grains

F0 = np.identity((3)) # applied strain
F0[1,1] = .5 # strain

phi0 = .2
phi1 = 50
a12 = .5
Fgb12 = np.identity(3) # grain boundary shear
Fgb12[0,1] = .5
v1dot= 0
v2dot = 0
dh12 = 0
dh21 = 0
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

def CalSig(F,mat):
    # calculates sigma from a given strain
    CC = GetC(mat)
    x = np.zeros((6,1))
    sig = np.zeros((3,3))
    x[0] = F[0,0]
    x[1] = F[1,1]
    x[2] = F[2,2]
    x[3] = F[1,2]
    x[4] = F[2,0]
    x[5] = F[0,1]
       
    temp = np.matmul(CC,x)
   
    sig[0,0] = temp[0]
    sig[1,1] = temp[1]
    sig[2,2] = temp[2]
    sig[1,2] = temp[3]
    sig[2,0] = temp[4]
    sig[0,1] = temp[5]
   
    sig[1,0] = sig[0,1]
    sig[2,1] = sig[1,2]
    sig[2,0] = sig[0,2]
   
    #sig[0,0] = mat[0,0]*F[0,0] + mat[0,1]*F[1,1]  + mat[0,2]*F[2,2] + 2*(mat[0,3]*F[1,2] + mat[0,4]*F[2,0] + mat[0,5]*F[0,1])
    #sig[1,1] = mat[1,0]*F[0,0] + mat[1,1]*F[1,1]  + mat[1,2]*F[2,2] + 2*(mat[1,3]*F[1,2] + mat[1,4]*F[2,0] + mat[1,5]*F[0,1])
    #sig[2,2] = mat[2,0]*F[0,0] + mat[2,1]*F[1,1]  + mat[2,2]*F[2,2] + 2*(mat[2,3]*F[1,2] + mat[2,4]*F[2,0] + mat[2,5]*F[0,1])
    #sig[1,2] = mat[3,0]*F[0,0] + mat[3,1]*F[1,1]  + mat[3,2]*F[2,2] + 2*(mat[3,3]*F[1,2] + mat[3,4]*F[2,0] + mat[3,5]*F[0,1])
    #sig[0,2] = mat[4,0]*F[0,0] + mat[4,1]*F[1,1]  + mat[4,2]*F[2,2] + 2*(mat[4,3]*F[1,2] + mat[4,4]*F[2,0] + mat[4,5]*F[0,1])
    #sig[0,1] = mat[5,0]*F[0,0] + mat[5,1]*F[1,1]  + mat[5,2]*F[2,2] + 2*(mat[5,3]*F[1,2] + mat[5,4]*F[2,0] + mat[5,5]*F[0,1])
    #sig[1,0] = sig[0,1]
    #sig[2,1] = sig[1,2]
    #sig[2,0] = sig[0,2]
    return sig

def Init(th,Vtot):
   
    # initalize variables for Volume, elastic modulus tensor, and eigenstrain
   
    C = np.array([[C11, C12, C12, 0, 0, 0],
                  [C12, C11, C12, 0, 0, 0],
                  [C12, C12, C11, 0, 0, 0],
                  [0, 0, 0, C44, 0, 0],
                  [0, 0, 0, 0, C44, 0],
                  [0, 0, 0, 0, 0, C44]])
   
    CC = np.zeros((N, 6, 6))
    Fgb = np.zeros((N, 3, 3))
    V = np.zeros(N)
   
    for i in range(N):
       
        th_g = th[i]

        R = np.array([[np.cos(th_g)**2, np.sin(th_g)**2, 0, 0, 0, -2*np.sin(th_g)*np.cos(th_g)],
                  [np.sin(th_g)**2, np.cos(th_g)**2 ,0 ,0 ,0 ,2*np.sin(th_g)*np.cos(th_g)],
                  [0, 0, 1, 0, 0, 0],
                  [0, 0, 0, np.cos(th_g), np.sin(th_g), 0],
                  [0 ,0 ,0 ,-np.sin(th_g), np.cos(th_g) ,0],
                  [np.sin(th_g)*np.cos(th_g) ,-np.sin(th_g)*np.cos(th_g) ,0 ,0 ,0 ,2*np.cos(th_g)**2-1]])
   
       
        CC[i] = np.matmul(np.matmul(np.transpose(R),C),R)
        Fgb[i] = np.identity(3)
        V[i] = Vtot/N
   
    return CC,Fgb,V

def GetC(m):
   
    Cs = np.zeros((6,6))
    for i in range(0,3):
        for j in range(0,3):
            Cs[i,j] = m[i,j]
    for i in range(0,6):
        for j in range(3,6):
            Cs[i,j] = 2*m[i,j]
    for i in range(3,6):
        for j in range(0,3):
            Cs[i,j] = 2*m[i,j]
    return Cs
def GetA(x):
   
    A_aux = np.zeros(((N+1)*6,6 ))
    A_temp = np.zeros(( (N+1)*6,6 ))
    C_aux = np.zeros(( 6, 6 ))
    A_end = np.concatenate((np.tile(-1*np.identity(6),(N,1)),np.zeros((6,6))),axis=0)
    for i in range(N):
       
        C_aux = GetC(C[i])
        if(i == 0): # starting condition
            A_aux = np.concatenate((C_aux,np.tile(np.zeros((6,6)),(N-1,1)),v[i]*np.identity(6)),axis=0)
            A_temp = A_aux
        elif(N-i == 1): # ending condition
            A_aux = np.concatenate((np.tile(np.zeros((6,6)),(i,1)),C_aux,v[i]*np.identity(6)),axis=0)
            A_temp = np.concatenate((A_temp,A_aux),axis = 1)
        else:
            A_aux = np.concatenate((np.tile(np.zeros((6,6)),(i,1)),C_aux,(np.tile(np.zeros((6,6)),(N-i-1,1))),v[i]*np.identity(6)),axis=0)
            A_temp = np.concatenate((A_temp,A_aux),axis = 1)
   
    A = np.concatenate((A_temp,A_end),axis=1)
    #A_aux[i] = np.concatenate((np.tile(C_aux[i],(N,1)),v[i]*np.identity(6)),axis=0)
    #b = A*x
    return A

def GetB():
   
    b_temp = np.zeros((9,1))
    b_aux = np.zeros((6,1))
   
    b3 = np.reshape(Vtot*F0,(9,1))
    b_end = np.array([b3[0],b3[4],b3[8],b3[5],b3[2],b3[1]])
   
    for i in range(N):
       
        b_temp = np.reshape(CalSig(Fgb[i],C[i]),(9,1))
        b_aux = np.array([b_temp[0],b_temp[4],b_temp[8],b_temp[5],b_temp[2],b_temp[1]])
        if(i == 0):
            b = b_aux
        else:
            b = np.vstack((b,b_aux))
   
    b = np.vstack((b,b_end))
    return b
def UnpackX(x):
   
    F = np.zeros((N,3,3))
    lamb = np.zeros((3,3))
    lamb[0,0] = x[N*6]
    lamb[1,1] = x[N*6 + 1]
    lamb[2,2] = x[N*6 + 2]
    lamb[1,2] = x[N*6 + 3]
    lamb[2,0] = x[N*6 + 4]
    lamb[0,1] = x[N*6 + 5]
   
    lamb[1,0] = lamb[0,1]
    lamb[2,1] = lamb[1,2]
    lamb[2,0] = lamb[0,2]
    for i in range(N):
        F[i,0,0] = x[i*6]
        F[i,1,1] = x[i*6 + 1]
        F[i,2,2] = x[i*6 + 2]
        F[i,1,2] = x[i*6 + 3]
        F[i,2,0] = x[i*6 + 4]
        F[i,0,1] = x[i*6 + 5]
       
        #F[i,1,0] = F[i,0,1]
        #F[i,2,1] = F[i,1,2]
        #F[i,2,0] = F[i,0,2]
       
    return F,lamb

def CalcAstar(grain,Fstar):
   
    A = 0
    if(grain == 1):
        Ftemp = Fstar[0] - Fgb[0]
        P = CalSig(Ftemp,C[0])  
        for i in range(0,3):
            for j in range(0,3):
                A = A + 0.5*Ftemp[i,j]*P[i,j]
    else:
        Ftemp = Fstar[1] - Fgb[1]
        P = CalSig(Ftemp,C[1])  
        for i in range(0,3):
            for j in range(0,3):
                A = A + 0.5*Ftemp[i,j]*P[i,j]
               
    return A

angle = np.array([np.radians(30/2),np.radians(-30/2),.1])
Vtot = 1.0

C,Fgb,v = Init(angle,Vtot)

V1[0] = v[0]
V2[0] = v[1]

# set up rotation matrix
theta = np.radians(0)
c, s = np.cos(theta), np.sin(theta)
rot2 = np.array(((c,0, s),(0,1,0), (-s,0, c)))
Fgb12 = np.matmul(np.matmul(rot2,Fgb12),np.transpose(rot2))

for t in range(tt):
   
    F0[0,1] = .5*np.sin(20*m.pi/(tt)*t) # applied strain

   
    #if(t == 381):
    #    aaa = 1
        
    A = GetA(1)
    b = GetB()
    x = linalg.solve(A,b)
    F,lamb = UnpackX(x)
   
    A1 = CalcAstar(1,F)
    A2 = CalcAstar(2,F)
   
    dh12 = -(A1 - A2 + phi0)/phi1
   
    dh21 = -(-A1 + A2 + phi0)/phi1
    gib1[t] = A1 - A2
   
    if (h12[t] < -1/2/a12 or h12[t] > 1/2/a12):
        dh12 = 0
        dh21 = 0
    if(dh12 < 0):
        dh12 = 0
    if(dh21 < 0):
        dh21 = 0
    if(dh21 > 0):
        dh12 = -dh21
        
    if(h12[t] > 0):
        # V1 gets bigger
        Fgb1 = Fgb12
        Fgb2 = np.identity(3)
        
    if(h12[t] < 0):
        # V1 gets smaller 
        Fgb1 = np.identity(3)
        Fgb2 = Fgb12
       
    h12[t+1] = h12[t] + dt*dh12 # euler integration of interface position and volume
    v1dot = a12*dh12
    v2dot = -a12*dh12
   
    V1[t+1] = V1[t] + v1dot*dt
    V2[t+1] = V2[t] + v2dot*dt
   
    time[t+1] = t*dt
   
    # stress strain curve
    temp1 = F[0] - Fgb[0]
    P1 = CalSig(temp1,C[0])
    temp2 = F[1] - Fgb[1]
    P2 = CalSig(temp2,C[1])
   
    stress[t+1] = lamb[0,1]#(P1[0,1]*V1[t] + P2[0,1]*V2[t])/Vtot
    strain[t+1] = F0[0,1]
   
    v[0] = V1[t]
    v[1] = V2[t]

   
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
ax2.plot(stress, strain)
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

# %% 

num = 1000
angle = np.linspace(0,2*m.pi,num)
helm1 = np.zeros(num)
helm2 = np.zeros(num)
F0 = np.zeros((3,3))
F0[0,1] = 1
F0[1,1] = 1

for i in range(num):
    
    ang = np.array([angle[i],-angle[i]])
    
    C,Fgb,v = Init(ang,Vtot)
    A = GetA(1)
    b = GetB()
    x = linalg.solve(A,b)
    F,lamb = UnpackX(x)
   
    A1 = CalcAstar(1,F)
    A2 = CalcAstar(2,F)
    
    helm1[i] = A1-A2
    #helm2[i] = A2
    
fig1, ax1 = plt.subplots()
ax1.plot(angle,helm1)
ax1.plot(angle,helm2)
ax1.set_xlabel("Angle (rad)")
ax1.set_ylabel("Forces")
ax1.set_title('Thermodynamic Forces with differing angle')
plt.grid(alpha=.7,linestyle='-.')
