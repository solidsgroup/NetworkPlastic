clear all
global N Vtot C Cinv v xi_inv alpha Fgb F0 edges D Aa I grain_map edge_map ss tt a

% load C:\Users\icrman\Documents\Matlab\mtex-5.8.0\userScripts\GrainData2.mat
load C:\Users\icrman\Documents\MATLAB\NetworkPlasticity\GrainData2.mat
N = 62; edges = 42;
Vtot = sum(Grain_vol);
% s = [1 1 1 1 2 2 3 3 3 4 4 4 5 5 6 7 7 8 8 9 10 10 10 10 11 11 12 13 13 14 14 15 16 17 17 17 17 18 19 19 20]';
% t = [2 3 4 5 3 18 4 17 18 5 8 11 6 8 7 8 10 10 11 7 12 15 16 9 10 12 13 10 14 10 15 16 9 4 11 19 20 17 18 20 11]';

 c = unique(s); % the unique values in s  
 grain_map = zeros(N-edges,3);
 edge_map = zeros(edges,1);
 for i = 1:length(c)
     grain_map(i,1) = c(i); % vertex number
     grain_map(i,2) = sum(s==c(i)); % number of times each unique value is repeated
     grain_map(i,3) = sum(grain_map(1:i-1,2)) + grain_map(i,1); % coresponding vertex index in N-arrays
 end
 counter = 2; temp = cumsum(grain_map(:,2)) + grain_map(:,1);
 counter2 = 1;
 for i = 1:edges
     edge_map(i) = counter; % coresponding edge index in N-arrays
     if(counter == temp(counter2))
         counter = counter + 1;
         counter2 = counter2 + 1;
     end
     counter = counter + 1;
 end

F0 = [1,.1 ,0;
      0, 1 ,0;
      0, 0 ,1];
Fgb = zeros(9,1,N);
xi_inv = zeros(9,9,N-edges); 
alpha= zeros(9,1,N-edges); 
[C,Cinv,v,Fgbt] = init(phi1,Phi,phi2,Grain_vol,Fpq,grain_map, edge_map);
F0 = reshape(F0, [9,1]);
th_gb = 0;
for i = 1:edges
    rot1 = [1 0 0;
        0 cos(th_gb) sin(th_gb);
        0 -sin(th_gb) cos(th_gb)];
    index = edge_map(i);
    Fgbt(:,:,index) = rot1'*Fgbt(:,:,index)*rot1;
end
for i = 1:N
    Fgb(:,:,i) = reshape(Fgbt(:,:,i), [9,1]);
end

%% 20-grains

dt = 2e-5;
tf = floor(1/dt);

ss = s;
tt = t;
D = digraph(ss,tt);
I = full(incidence(D));

dhpq = 0; dhqp = 0; doth = zeros(edges,1);
h = zeros(tf,edges);
V = zeros(tf,N-edges);
for i = 1:N-edges
    p = grain_map(i,3);
    V(1,i) = v(p);
end
vol = Grain_vol';
V0 = vol; V(1,:) = V0;

% stress = zeros(tf,1);
% strain = zeros(tf,1);
time = zeros(tf,1);

AA = zeros(tf,edges);

n = .3;
m = 20;
phi_h1 = .5;
phi_h2 = 50;
for i = 1:edges
    p = grain_map(ss(i),1); q = grain_map(tt(i),1);
    vv = [Grain_vol(p);Grain_vol(q)];
    kappa = (rand(edges,1)*2+.11)/min(vv);
end

phi_hs = zeros(tf,3); 
phi11 = phi11*5; 
phi0 = rand(edges,1)*.1+.7;
for t = 1:tf
    
    F0(4) = .1*t/tf;
    
    Fstar = FstarAnalytic();
    P0 = GetP0(Fstar);
    for edge = 1:edges
    
        p = grain_map(ss(edge),3); q = grain_map(tt(edge),3);  
        [dFpq, dFqp] = Gethdot(Fstar,P0,edge,vol);

        phi_h = phi_h1*abs(h(t,edge)/kappa(edge))^n + phi_h2*abs(h(t,edge)*kappa(edge))^m;
        dhpq = -1/phi11(edge)*(dFpq + phi0(edge) + phi_h);
        dhqp = -1/phi11(edge)*(dFqp - phi0(edge) - phi_h);

        AA(t+1,edge) = dFqp;
        phi_hs(t+1,edge) = phi_h;
        
        veff = h(t,edge)*a(edge);
        if(veff < -vol(ss(edge)) || veff > vol(tt(edge)))
            dhpq = 0;
            dhqp = 0;
        end
        if(veff > vol(ss(edge)) || veff < -vol(tt(edge)))
            dhpq = 0;
            dhqp = 0;
        end
        if(dhpq < 0); dhpq = 0; end
        if(dhqp > 0); dhqp = 0;
        elseif(dhqp < 0); dhpq = dhqp;end

        pq = edge_map(edge);
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
    end
    for i = 1:edges
        h(t+1,i) = h(t,i) + dh(i);
    end

    time(t+1) = t*dt;

    stress(t+1) = P0(4);
    strain(t+1) = F0(4);
end
%%
figure(2)

hold on
plot(time,h,'--','LineWidth',1.3)
xlabel('Time')
ylabel('Boundary displacment')
% legend('h_{13}','h_{35}','h_{51}','Location','NorthWest')
grid on
hold off

figure(3)
hold on
plot(time,V,'LineWidth',1.3)
xlabel('Time')
ylabel('Volume')
% legend('V_1','V_3','V_5','Location','SouthEast')
grid on
hold off
figure(4)
plot(strain,stress,'LineWidth',1.3)
xlabel('Strain')
ylabel('Stress')
grid on
figure(5)
hold on
plot(AA)
% plot([1 tf],[1 1]*phi0)
% plot([1 tf],-[1 1]*phi0)
title('dA^*/dh_{pq}')
hold off

figure(6)
hold on
plot(phi_hs(:,:))
title('dA^*/dh_{pq}')
hold off
        
function [vol, h] = UpdateVol(dh,dt)
global N edges Vtot C Cinv v Fgb F0  a D Aa I grain_map edge_map

    h = dh.*dt; ha = abs(h);
    Ah = full(adjacency(D,h));
    Aha = full(adjacency(D,ha));
    
    for i = 1:edges % calculate Vpq
       pq = edge_map(i);
       v(pq) = v(pq) + a(i)*h(i);
    end
    
    gimal = -Aa' + Aa;
    gimal_p = Aa' + Aa;
    gimal_h = Ah + Ah';
    gimal_ha = Aha + Aha';
    
    dV = -1/2*(gimal_ha*gimal_p + gimal*gimal_h);
    dV = diag(dV);
    
    for i = 1:N-edges % add dv to vertices
       p = grain_map(i,3);
%        if(v(p) <= 0)
%            dV(i) = 0;
%        end
       v(p) = v(p) + dV(i);
    end
    
    vol = -I*h; % volume of grains + quasi-grains
  
%   for i = 1:2:N 
%       if(v(i) < 0)
%          v(i) = 0;
%       end
%   end
end

function [dhpq, dhqp] = Gethdot(Fstr,P0,i,vol)
    global N Vtot C Cinv  Fgb F0 a12 xi_inv edges alpha grain_map edge_map ss tt
   
    Fstardh_ana = zeros(3,3,N-edges);
    dhpq = 0; dhqp = 0; I = [1;0;0;0;1;0;0;0;1];
    % ****** vm < 0 ******
    pq = edge_map(i);
    
    counter = 1;
    for pp = 1:N-edges
        p = grain_map(ss(pp),3); q = grain_map(tt(pp),3);  

        C_diff = eye(9) - Cinv(:,:,p)*C(:,:,q);
        dfstar_dhpq = -xi_inv(:,:,q)^2*(C_diff)*alpha(:,:,q) + xi_inv(:,:,q)*(I-Fgb(:,:,pq) + C_diff*Fgb(:,:,q));

        Fstardh_ana(:,:,counter) = reshape(dfstar_dhpq,[3,3]);
        counter = counter + 1;
            
    end

    Fstardh = zeros(3,3);
    for j = 1:N-edges
        Fstardh = Fstardh + Fstardh_ana(:,:,j)*vol(j);
    end

    p = grain_map(ss(i),3); q = grain_map(tt(i),3);  
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
    global N edges Vtot C Cinv v Fgb F0 xi_inv alpha grain_map
   
    Fstr = zeros(9,1,N-edges);
    for i = 1:N-edges
        index = grain_map(i,3);
       
        xi = v(index)*eye(9);
        alpha_temp = Vtot*F0;
        for j = 1:N            
            if(i == j);continue; end        
            xi = xi + abs(v(j))*Cinv(:,:,j)*C(:,:,index);
            alpha_temp = alpha_temp - v(j)*Fgb(:,:,j);
        end
        alpha(:,:,i) = alpha_temp + (xi-v(index)*eye(9))*Fgb(:,:,index);
        if(v(i) == 0)
            Fstr(:,:,i) = Fgb(:,:,index);
            xi_inv(:,:,i) = eye(9);
        else
            xi_inv(:,:,i) = xi\eye(9);
            Fstr(:,:,i) = xi_inv(:,:,i)*alpha(:,:,i);
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

    if(i == j); val = 1;
    else; val = 0;end
end
function P = GetP0(Fstr)
    global N edges Vtot C Fgb v grain_map
    
    P = zeros(9,1);
    for i = 1:N-edges
        index = grain_map(i,3);
        if(v(index) <= 0)
            continue
        else
            P = P + C(:,:,index)*(Fstr(:,:,i)-Fgb(:,:,index));
            break
        end
    end
end
function [CC,CCinv,v,Fgb] = init(p1,p,p2,vol,Fpq,v_map,e_map)
    global N Vtot edges
   
    CC = zeros(9,9,N);
    CCinv = CC; v = zeros(N,1); Fgb = zeros(3,3,N);
    for i = 1:N-edges
      index = v_map(i,3);
      [CC(:,:,index),CCinv(:,:,index)] = SetC(p1(i),p(i),p2(i));
        
      v(index) = vol(i);
      Fgb(:,:,index) = eye(3);
    end
    for i = 1:edges
        index = e_map(i);
        Fgb(:,:,index) = Fpq(:,:,i);
%         Fgb(:,:,index) = eye(3);
%         Fgb(1,2,index) = rand(1)*2-1;
    end
end
