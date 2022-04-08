clear all
global N Vtot C_v Cinv_v v_v xi_inv_v alpha_v Fgb_v C_e Cinv_e v_e Fgb_e F0 edges D Aa I ss tt

load C:\Users\icrma\Documents\Matlab\mtex-5.8.0\userScripts\GrainData2.mat
N = 61; edges = 41;
Vtot = sum(Grain_vol);

F0 = [1,.1 ,0;
      0, 1 ,0;
      0, 0 ,1];
Fgb_v = zeros(9,1,N-edges);Fgb_e = zeros(9,1,edges);
xi_inv_v = zeros(9,9,N-edges); 
alpha_v = zeros(9,1,N-edges); 
[C_v,Cinv_v,v_v] = init(phi1,Phi,phi2,Grain_vol);
F0 = reshape(F0, [9,1]);
th_gb = 0;Fgbt = zeros(3,3,edges);
for i = 1:edges
    rot1 = [1 0 0;
        0 cos(th_gb) sin(th_gb);
        0 -sin(th_gb) cos(th_gb)];
    Fgbt(:,:,i) = rot1'*Fpq(:,:,i)*rot1;
    Fgb_e(:,:,i) = reshape(Fgbt(:,:,i), [9,1]);
end
for i = 1:N-edges
    Fgb_v(:,:,i) = reshape(eye(3), [9,1]);
end

%% Tricrystal

dt = 1e-4;
tf = floor(1/dt);

ss = s;
tt = t;
D = digraph(ss,tt);

dhpq = 0; dhqp = 0; doth = zeros(edges,1);
h = zeros(tf,edges);
V = zeros(tf,N-edges);
counter = 1;
for i = 1:2:N-edges
    V(1,counter) = v(i);
    counter = counter + 1;
end
vol = [v(1); v(3); v(5)];
V0 = vol; V(1,:) = V0;
phi1 = 0.006;
phi0 = 0.2;

stress = zeros(tf,1);
strain = zeros(tf,1);
time = zeros(tf,1);

AA = zeros(tf,edges);

n = .3;
m = 20;
phi_h1 = .5;
phi_h2 = 10;
kappa = [1.6/0.3750; 2.36/0.3750; 1.1/0.3750];

phi_hs = zeros(tf,3); 
for t = 1:tf
    
    F0(4) = .2*t/tf;
    
    Fstar = FstarAnalytic();
    P0 = GetP0(Fstar);
    for edge = 1:edges
    
        q = ss(edge)*2-1; p = tt(edge)*2-1;
        [dFpq, dFqp] = Gethdot(Fstar,P0,edge,vol);

        phi_h = phi_h1*abs(h(t,edge))^n + phi_h2*abs(h(t,edge)*kappa(edge))^m;
        dhpq = -1/phi1*(dFpq + phi0 + phi_h);
        dhqp = -1/phi1*(dFqp - phi0 - phi_h);

        AA(t+1,edge) = dFqp;
        phi_hs(t+1,edge) = phi_h;

        if(h(t,edge) < -V0(ss(edge))/a(edge) || h(t,edge) > V0(tt(edge))/a(edge))
            dhpq = 0;
            dhqp = 0;
        end
        if(dhpq < 0); 
            dhpq = 0;
        end

        if(dhqp > 0); 
            dhqp = 0;
        elseif(dhqp < 0); 
            dhpq = dhqp;
        end
        pq = edge*2;
        if(v(pq) < 0)
            C(:,:,pq) = C(:,:,p);
            Cinv(:,:,pq) = Cinv(:,:,p);
        elseif(v(pq) > 0)
            C(:,:,pq) = C(:,:,q);
            Cinv(:,:,pq) = Cinv(:,:,q);
        else
            C(:,:,pq) = zeros(9,9);
            Cinv(:,:,pq) = zeros(9,9);
        end

        doth(edge) = dhpq;
    end

    [dv,dh] = UpdateVol(doth,dt);
    vol = vol + dv;
    for i = 1:N-edges
        V(t+1,i) = vol(i);
        h(t+1,i) = h(t,i) + dh(i);
    end

    time(t+1) = t*dt;

    stress(t+1) = P0(4);
    strain(t+1) = F0(4);
end

figure(2)

hold on
plot(time,h(:,1),'--','LineWidth',1.3)
plot(time,h(:,2),'LineWidth',1.3)
plot(time,h(:,3),'LineWidth',1.3)
xlabel('Time')
ylabel('Boundary displacment')
legend('h_{13}','h_{35}','h_{51}','Location','NorthWest')
grid on
hold off

figure(3)
hold on
plot(time,V(:,1),'LineWidth',1.3)
plot(time,V(:,2),'LineWidth',1.3)
plot(time,V(:,3),'LineWidth',1.3)
xlabel('Time')
ylabel('Volume')
legend('V_1','V_3','V_5','Location','SouthEast')
grid on
hold off
figure(4)
plot(strain,stress,'LineWidth',1.3)
xlabel('Strain')
ylabel('Stress')
grid on
figure(5)
hold on
plot(AA(:,1))
plot(AA(:,2))
plot(AA(:,3))
plot([1 tf],[1 1]*phi0)
plot([1 tf],-[1 1]*phi0)
title('dA^*/dh_{pq}')
hold off

figure(6)
hold on
plot(phi_hs(:,1))
plot(phi_hs(:,2))
plot(phi_hs(:,3))
title('dA^*/dh_{pq}')
hold off
        
function [vol, h] = UpdateVol(dh,dt)
global N Vtot C Cinv v Fgb F0  a D Aa I

    h = dh.*dt; ha = abs(h);
    Ah = full(adjacency(D,h));
    Aha = full(adjacency(D,ha));
    counter = 1;
    for i = 2:2:N
       v(i) = v(i) + a(counter)*h(counter); 
       counter = counter + 1;
    end
    
    gimal = -Aa' + Aa;
    gimal_p = Aa' + Aa;
    gimal_h = Ah + Ah';
    gimal_ha = Aha + Aha';
    
    dV = -1/2*(gimal_ha*gimal_p + gimal_h*gimal);
    dV = diag(dV);
    counter = 1;
    for i = 1:2:N
       v(i) = v(i) + dV(counter);
       counter = counter + 1;
    end
    
    vol = I*h;
  
  for i = 1:2:N
      

      if(v(i) < 0)
         v(i) = 0;
      end
  end
end

function [dhpq, dhqp] = Gethdot(Fstr,P0,i,vol)
    global N Vtot C Cinv  Fgb F0 a12 xi_inv edges alpha ss tt
   
    Fstardh_ana = zeros(3,3,N-edges);
    dhpq = 0; dhqp = 0; I = [1;0;0;0;1;0;0;0;1];
    % ****** vm < 0 ******
    pq = i*2;
    counter = 1;
    qs = ss*2 - 1; ps = tt*2 - 1; 
    for pp = 1:edges
        
        p = ps(pp); q = qs(pp);  
        C_diff = eye(9) - Cinv(:,:,p)*C(:,:,q);
        dfstar_dhpq = -xi_inv(:,:,q)^2*(C_diff)*alpha(:,:,q) + xi_inv(:,:,q)*(I-Fgb(:,:,pq) + C_diff*Fgb(:,:,q));

        Fstardh_ana(:,:,counter) = reshape(dfstar_dhpq,[3,3]);
        counter = counter + 1;
            
    end

    if(mod(N,2) == 1)
        
        C_diff = eye(9) - Cinv(:,:,q)*C(:,:,p);
        dfstar_dhpq = -xi_inv(:,:,p)^2*(C_diff)*alpha(:,:,p) + xi_inv(:,:,p)*(I - Fgb(:,:,pq) + C_diff*Fgb(:,:,p));

        Fstardh_ana(:,:,end) = reshape(dfstar_dhpq,[3,3]);

    end


    counter = 1; Fstardh = zeros(3,3);
    for j = 1:2:(N-edges)*2
        
        Fstardh = Fstardh + Fstardh_ana(:,:,counter)*vol(counter);
        counter = counter + 1;
    end
   

    p = ps(i); q = qs(i); 
    Fdiff = Fstr(:,:,q)-Fstr(:,:,p);

    Fdiff = reshape(Fdiff,[3,3]);
    P0 = reshape(P0,[3,3]);
   
    for i = 1:3
        for j = 1:3
            dhpq = dhpq + (1/2*Fdiff(i,j) + Fstardh(i,j))*P0(i,j);
%             dhpq = dhpq + (1/2*Fdiff(i,j) + vp/Vtot*Fstardh_p_ana(i,j) + vq/Vtot*Fstardh_q_ana(i,j))*P(i,j);
%             dhpq = dhpq +  Fstardh_q(i,j)*P(i,j);
        end
    end
   
    for i = 1:3
        for j = 1:3
          dhqp = dhqp + (-1/2*Fdiff(i,j) + Fstardh(i,j))*P0(i,j);
%             dhqp = dhqp + Fstardh_p(i,j)*P(i,j);
        end
    end
   
end


function [Fstr] = FstarAnalytic()
global N Vtot C_v Cinv_v v_v xi_inv_v alpha_v Fgb_v C_e Cinv_e v_e  Fgb_e F0 edges 
   
    Fstr = zeros(9,1,N);
%     for i = 1:edges
%         
%         xi_e = 
%     end
    for i = 1:N-edges
        xi = v(i)*eye(9);
        alpha_temp = Vtot*F0;
        for j = 1:N-edges            
            if(i == j)
                continue
            end

            xi = xi + abs(v_v(j))*Cinv_v(:,:,j)*C_v(:,:,i);
            alpha_temp = alpha_temp - v(j)*Fgb_v(:,:,j);
        end
        alpha_v(:,:,i) = alpha_temp + (xi-v(i)*eye(9))*Fgb_v(:,:,i);
        if(v(i) == 0)
            Fstr(:,:,i) = Fgb_v(:,:,i);
            xi_inv_v(:,:,i) = eye(9);
        else
            xi_inv_v(:,:,i) = xi\eye(9);
            Fstr(:,:,i) = xi_inv_v(:,:,i)*alpha_v(:,:,i);
        end
       
    end
   
end
function [C,Cinv] = SetC(p1,p,p2)

    C11 = 169.3097; C12 = 122.5; C44 = 76;
%     C11 = 169.3097; C12 = 87.2201 ; C44 = 41.0448;

    S11 = (C11+C12)/((C11-C12)*(C11+2*C12));
    S12 = -C12/((C11-C12)*(C11+2*C12));
    S44 = 1/C44;
    C0 = C11 - C12 -2*C44;
    S0 = S11-S12-1/2*S44;
   
    C = zeros(9,9);
    Cinv = zeros(9,9);
   
    r11 = cos(p1)*cos(p2) - cos(p)*sin(p1)*sin(p2);
    r12 = -cos(p1)*sin(p2) - cos(p)*sin(p1)*cos(p2);
    r13 = sin(p)*sin(p1);

    r21 = cos(p2)*sin(p1) + cos(p)*cos(p1)*sin(p2);
    r22 = cos(p)*cos(p1)*cos(p2) -  sin(p1)*sin(p2);
    r23 = -sin(p)*cos(p1);

    r31 = sin(p)*sin(p2);
    r32 = sin(p)*cos(p2);
    r33 = cos(p);
    

    e1 = [r11;r12;r13];
    e2 = [r21;r22;r23];
    e3 = [r31;r32;r33];
    d_ijkl = kron(e1*e1',e1*e1') + kron(e2*e2',e2*e2') + kron(e3*e3',e3*e3');
    counter = 1;
    I = 1; J = 1;
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    d_ij = del(i,j);
                    d_kl = del(k,l);
                    d_ik = del(i,k);
                    d_jl = del(j,l);
                    d_jk = del(j,k);
                    d_il = del(i,l);
                   
                    C(I,J) = C12*d_ij*d_kl + C44*(d_ik*d_jl*2 + d_il*d_jk*0) + C0*d_ijkl(I,J);
                    Cinv(I,J) = S12*d_ij*d_kl + 1/4*S44*(d_ik*d_jl*2 + d_il*d_jk*0) + S0*d_ijkl(I,J);
                    if(mod(counter,9) == 0)
                        I = I + 1;
                        counter = 0;
                    end
                    counter = counter + 1;
                    J = counter;
                end
            end
        end
    end
   
end

function val = del(i,j)
   
    if(i == j)
        val = 1;
    else
        val = 0;
    end
end
function P = GetP0(Fstr)
    global N Vtot C Fgb v
    
    P = zeros(9,1);
    for i = 1:N
        
        P = P + abs(v(i))/Vtot*C(:,:,i)*(Fstr(:,:,i)-Fgb(:,:,i));
    end
end
function [CC,CCinv,v,Fgb] = init(p1,p,p2,vol)
    global N Vtot edges
   
    CC = zeros(9,9,N-edges);
    CCinv = CC; v = zeros(N-edges,1); Fgb = zeros(3,3,N-edges);
    for i = 1:N-edges
       
      [CC(:,:,i),CCinv(:,:,i)] = SetC(p1(i),p(i),p2(i));
        
      v(i) = vol(i);
      Fgb(:,:,i) = eye(3);

    end

end