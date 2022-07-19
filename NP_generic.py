# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
import math as m

np.set_printoptions(edgeitems=30, linewidth=100000, 
    formatter=dict(float=lambda x: "%.5g" % x))
# input values
s = np.array([1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 10, 10, 11, 11, 12, 13, 13, 14, 14, 15, 16, 17, 17, 17, 17, 18, 19, 19, 20])-1
t = np.array([2, 3, 4, 5, 3, 18, 4, 17, 18, 5, 8, 11, 6, 8, 7, 8, 10, 10, 11, 7, 12, 15, 16, 9, 10, 12, 13, 10, 14, 10, 15, 16, 9, 4, 11, 19, 20, 17, 18, 20, 11])-1

N = 20 + len(s)
edges = len(s)
Fpq = np.array([[1,0.5,0],[0,1,0],[0,0,1]])
vol = np.ones((20,1))*1/20
Vtot = 1
a = np.linspace(1, edges,edges)
# initalize values
F0 = np.array([[1,0.25,0],[0,1,0],[0,0,1]])
F0 = F0.transpose().reshape((9,1))
alpha = np.zeros([N-edges,9,1])
xi_inv = np.zeros([N-edges,9,9])
Fgb = np.zeros([N,9,1])

v_map = np.arange(0,N-edges)
e_map = np.arange(N-edges,N)

def adjacency(s,t,weight):
    A = np.zeros([N-edges, N-edges])
    for i in range(0, edges):
        for j in range(0, edges):
            A[s[i],t[i]] = weight[i]
    return A
def incidence(s,t,weight):
    I = np.zeros([N-edges,edges])
    for i in range(0, edges):
        temp = I[:,i]
        temp[s[i]] = -weight[i]
        temp[t[i]] = weight[i]
        I[:,i] = temp
    return I      
def init(p1, p, p2, vol,Fpq):
    CC = np.zeros([N,9,9])
    CCinv = np.zeros([N,9,9])
    Fgb = np.zeros([N,3,3])
    v = np.zeros([N,1])
    for i in range(0, N-edges):
        # assign index value from v map
        index = v_map[i]
        CC[index], CCinv[index] = SetC(p1[i], p[i], p2[i])
        v[index] = vol[i]
        Fgb[index] = np.eye(3)
    for i in range(0, edges):
        # assign index value from e map
        index = e_map[i]
        Fgb[index] = Fpq
        
    return CC, CCinv, v, Fgb
def delta(i,j):
    if(i == j):
        val = 1
    else:
        val = 0
    return val
def SetC(p1,p,p2):
    C11 = 169.3097
    C12 = 122.5
    C44 = 76.0
    
    S11 = (C11 + C12)/( (C11 - C12)*(C11 + 2*C12) )
    S12 = -C12/( (C11 - C12)*(C11 + 2*C12) )
    S44 = 1/C44
    C0 = C11 - C12 - 2*C44
    S0 = S11 - S12 -1/2*S44
    
    C = np.zeros([9,9])
    Cinv = np.zeros([9,9])
    
    r11 = np.cos(p1)*np.cos(p2) - np.cos(p)*np.sin(p1)*np.sin(p2)
    r12 = -np.cos(p1)*np.sin(p2) - np.cos(p)*np.sin(p1)*np.cos(p2)
    r13 = np.sin(p)*np.sin(p1)

    r21 = np.cos(p2)*np.sin(p1) + np.cos(p)*np.cos(p1)*np.sin(p2)
    r22 = np.cos(p)*np.cos(p1)*np.cos(p2) -  np.sin(p1)*np.sin(p2)
    r23 = -np.sin(p)*np.cos(p1)

    r31 = np.sin(p)*np.sin(p2)
    r32 = np.sin(p)*np.cos(p2)
    r33 = np.cos(p)
    
    e1 = np.array([r11,r12,r13])
    e2 = np.array([r21,r22,r23])
    e3 = np.array([r31,r32,r33])
    
    d_ijkl = np.kron(np.tensordot(e1,e1,axes=0),np.tensordot(e1,e1,axes=0)) + np.kron(np.tensordot(e2,e2,axes=0),np.tensordot(e2,e2,axes=0)) + np.kron(np.tensordot(e3,e3,axes=0),np.tensordot(e3,e3,axes=0))
    counter = 1
    I = 0
    J = 0
    for i in range(0, 3):
        for j in range(0, 3):
            for k in range(0, 3):
                for l in range(0, 3):
                    d_ij = delta(i, j)
                    d_kl = delta(k, l)
                    d_ik = delta(i, k)
                    d_jl = delta(j, l)
                    #d_jk = delta(j, k)
                    #d_il = delta(i, l)
                    
                    C[I,J] = C12*d_ij*d_kl + 2*C44*d_ik*d_jl + C0*d_ijkl[I,J]
                    Cinv[I,J] = S12*d_ij*d_kl + 1/2*S44*d_ik*d_jl + S0*d_ijkl[I,J]
                    if(counter%9 == 0):
                        I += 1
                        counter = 0
                        
                    J = counter
                    counter += 1
                    
    return C, Cinv
def FstarAnalytic():
    Fstr = np.zeros([N-edges,9,1])
    for i in range(0, N-edges):
        # get (i) index from v map
        index = v_map[i]
        xi = v[index]*np.eye(9)
        alpha_temp = Vtot*F0
        for j in range(0, N):
            if(i == j):
                continue
            xi += abs(v[j])*np.matmul(Cinv[j],C[index])
            alpha_temp -= v[j]*Fgb[j]
        temp = (xi - v[index]*np.eye(9))
        alpha[i] = alpha_temp + np.matmul(temp,Fgb[index])
        if(v[i] == 0):
            Fstr[i] = Fgb[index]
            xi_inv[i] = np.eye(9)
        else:
            xi_inv[i] = np.linalg.inv(xi)
            Fstr[i] = np.matmul(xi_inv[i],alpha[i])
    return Fstr
def GetP0(Fstr,vol):
    for i in range(0, N-edges):
        index = v_map[i]
        if(vol[i] <= 0):
            continue
        else:
            P = np.matmul(C[index],(Fstr[index]-Fgb[index]))
            break
    return P
def Gethdot(Fstr,P0,edge,vol):
    
    Fstardh_ana_pq = np.zeros([N-edges,3,3])
    Fstardh_ana_qp = np.zeros([N-edges,3,3])
    dhpq = 0
    dhqp = 0
    II = np.array([[1],[0],[0],[0],[1],[0],[0],[0],[1]])
    
    pq = e_map[edge]
    for pp in range(0, N-edges):
        p = v_map[s[pp]]
        q = v_map[t[pp]]
        
        C_diff = np.eye(9) - np.matmul(Cinv[p],C[q])
        xiinv2 = np.matmul(xi_inv[q],xi_inv[q])
        alpha_dh = II - Fgb[pq] + np.matmul(C_diff,Fgb[q])
        dfstar_dhpq = -np.matmul(np.matmul(xiinv2,C_diff),alpha[q]) + np.matmul(xi_inv[q],alpha_dh)
        alpha_dh = -II + Fgb[pq] + np.matmul(C_diff,Fgb[q])
        dfstar_dhpq = np.matmul(np.matmul(xiinv2,C_diff),alpha[q]) + np.matmul(xi_inv[q],alpha_dh)
        
        Fstardh_ana_pq[pp] = dfstar_dhpq.reshape((3,3)).transpose()
        Fstardh_ana_qp[pp] = dfstar_dhpq.reshape((3,3)).transpose()
        
    Fstardh_pq = np.zeros([3,3])
    Fstardh_qp = np.zeros([3,3])
    
    for j in range(0, N-edges):
        Fstardh_pq += Fstardh_ana_pq[j]*vol[j]
        Fstardh_qp += Fstardh_ana_qp[j]*vol[j]
        
    q = s[edge]
    p = t[edge]
    
    F_diff = Fstr[q] - Fstr[p]
    F_diff = F_diff.reshape((3,3)).transpose()
    
    P0 = P0.reshape((3,3)).transpose()
    for i in range(0, 3):
        for j in range(0, 3):
            dhpq += (1/2*F_diff[i,j] + Fstardh_pq[i,j])*P0[i,j]     
    for i in range(0, 3):
        for j in range(0, 3):
            dhqp += (-1/2*F_diff[i,j] + Fstardh_qp[i,j])*P0[i,j]
            
    return dhpq, dhqp
def UpdateVol(dh, dt):
    h = dh*dt
    ha = abs(dh)
    
    Ah = adjacency(s, t, h)
    Aha = adjacency(s, t, ha)
    
    for i in range(0, edges):
        pq = e_map[i]
        v[pq] = v[pq] + a[i]*h[i]
        
    gimal = -Aa.transpose() + Aa
    gimal_p = Aa.transpose() + Aa
    gimal_h = -Ah.transpose() + Ah
    gimal_ha = -Aha.transpose() + Aha
    
    dV = -1/2* (np.matmul(gimal_ha,gimal_p) + np.matmul(gimal_h,gimal))
    dV = np.diagonal(dV)
    
    for i in range(0, N-edges):
        p = v_map[i]
        v[p] = v[p] + dV[p]
        
    vol = -np.matmul(I,dh)
    return vol, h

A = adjacency(s,t,a)
I = incidence(s,t,a)
p1 = np.random.standard_normal(N-edges)
p = np.random.standard_normal(N-edges)
p2 = np.random.standard_normal(N-edges)
C,Cinv, v, Fgbt = init(p1,p,p2,vol,Fpq)
for i in range(0, N):
    Fgb[i] = Fgbt[i].transpose().reshape((9,1))
    
    
#Fstar = FstarAnalytic()

#Fstar[0]*0.5 + Fstar[1]*0.5 - Vtot*F0
#np.matmul(C[0],(Fstar[0]-Fgb[0])) - np.matmul(C[1],(Fstar[1]-Fgb[1]))

dt = 2e-5
tf = np.floor(1/dt)
endt = 4005
time = np.zeros([endt,1])
Aa = adjacency(s, t, a)
I = incidence(s, t, a)

dhpq = 0
dhqp = 0;
doth = np.zeros([edges,1])
h = np.zeros([endt,edges])
V = np.zeros([endt,N-edges])
for i in range(0, N-edges):
    p = v_map[i]
    V[0,p]
#vol = Grain_vol

stress = np.zeros([endt,1])
strain = np.zeros([endt,1])

n = 0.3
mm = 20
phi_h1 = 0.2
phi_h2 = 2.0
kappa = np.zeros([edges,1])
for i in range(0, edges):
    p = v_map[s[i]]
    q = v_map[t[i]]
    kappa[i] = 1/np.random.uniform(0.3,0.5,1)*np.min([v[p],v[q]])
    
aeff = a
phi1 = np.random.uniform(0.01,0.08,edges)
phi0 = np.random.uniform(2.5,2.6,edges)

dv = np.zeros([edges,1])
dh = np.zeros([edges,1])
for tt in range(0, endt):
    F0[3] = 1*tt/tf
    
    Fstar = FstarAnalytic()
    P0 = GetP0(Fstar, vol)
    for i in range(0, edges):
        
        q = v_map[s[i]]
        p = v_map[t[i]]
        
        dFpq, dFqp = Gethdot(Fstar,P0,i,vol)
        
        phi_h = phi_h1*abs(h[tt,i]*kappa[i]*aeff[i])**n + phi_h2*abs(h[tt,i]*kappa[i]*aeff[i])**mm
        dhpq = -1/phi1[i]*(dFpq + phi0[i] + phi_h)
        dhqp = -1/phi1[i]*(dFpq - phi0[i] - phi_h)
        
        if(vol[q] <= 0 or vol[p] <= 0):
            dhpq = 0
            dhqp = 0
        if(dhpq < 0):
            dhpq = 0
        if(dhqp > 0):
            dhqp = 0
        elif(dhqp < 0):
            dhpq = dhqp
        pq = e_map[i]
        if(v[pq] < 0):
            C[pq] = C[p]
            Cinv[pq] = Cinv[p]
        elif(v[pq] > 0):
             C[pq] = C[q]
             Cinv[pq] = Cinv[q]  
        else:
             C[pq] = np.zeros([9,9])
             Cinv[pq] = np.zeros([9,9])
        doth[i] = dhpq
        bea = -2*abs(h[tt,i])/0.03
        
        if(bea == 0):
            continue
        else:
            aeff[i] = 0.5*(a[i]*m.sqrt(bea**2*a[i]**2 + 1) + m.asinh(bea*a[i])/bea)
                           
    dv, dh = UpdateVol(doth, dt)
    vol = vol + dv
    for j in range(0, N-edges):
        V[tt,j] = vol[j]
    for j in range(0, edges):
        h[tt,j] = h[tt,j] + dh[j]
    time[tt] = tt*dt
    
    stress[tt] = P0[3]
    strain[tt] = F0[3]
    
    if(tt%100 == 0):
        print(tt)

plt.plot(strain, stress)
plt.show()

   