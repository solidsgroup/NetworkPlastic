clear all
global N Vtot C Cinv v Fgb F0  a13 a35 a51 xi_inv alpha edges a D Aa I

N = 6; edges = 3;
Vtot = 1;

F0 = [1,.1 ,0;
      0, 1 ,0;
      0, 0 ,1];
Fgb = zeros(9,1,N);
xi_inv = zeros(9,9,N);
alpha = zeros(9,1,N);
angle = [0.2 0.2 0.3 0.3 .4 .4];
% angle = [0.2, 0.2, .3];
[C,Cinv,v,Fgbt] = init(angle);
F0 = reshape(F0, [9,1]);
th_gb = atan(.5); 
rot1 = [1 0 0;
    0 cos(th_gb) sin(th_gb);
    0 -sin(th_gb) cos(th_gb)];
Fgbt(:,:,6) = rot1'*Fgbt(:,:,6)*rot1;
Fgbt(:,:,4) = rot1*Fgbt(:,:,4)*rot1';
for i = 1:N
    Fgb(:,:,i) = reshape(Fgbt(:,:,i), [9,1]);
end

% delu = (F0-Fgb(:,:,1))*(C(:,:,1)-C(:,:,3))*(F0-Fgb(:,:,1))
% vm1 = 0.1; vm2 = .2;
% v(1) = 1/3; v(2) = vm1; v(3) = (1/3 - vm1); v(4) = vm2; v(5) = 1/3-vm2; 
% Fstr = FstarAnalytic();
% 
% p = zeros(9,1);
% vv = zeros(9,1);
% for i = 1:N-1
%     p = p + C(:,:,i)*(Fstr(:,:,i) - Fgb(:,:,i));
%     vv = vv + v(i)*Fstr(:,:,i);
%     C_avg = C_avg + v(i)*C(:,:,i);
% end
% P_test = p - (N-1)*C(:,:,end)*(Fstr(:,:,end) - Fgb(:,:,end))
% V_test = vv + v(end)*Fstr(:,:,end) - Vtot*F0

%% Tricrystal

dt = 1e-4;
tf = floor(1/dt);
a13 = 1/2; a35 = sqrt(5)/4; a51 = sqrt(5)/4;
v(1) = 3/16*2; v(3) = 3/16*2; v(5) = 1/8*2; % rect. base = 1 height = .5

s = [1 2 3];
t = [2 3 1];
a = [a13 a35 a51];
D = digraph(s,t);
Aa = full(adjacency(D,a));
I = -full(incidence(D))*diag(a);

dh13 = 0; dh31 = 0; % first boundary
dh35 = 0; dh53 = 0; % second boundary
dh51 = 0; dh15 = 0; % third boundary

phi1 = .006;
phi0 = 0.2;

h13 = zeros(tf,1); h35 = zeros(tf,1); h51 = zeros(tf,1);
V1 = zeros(tf,1); V1(1) = v(1);
V3 = zeros(tf,1); V3(1) = v(3);
V5 = zeros(tf,1); V5(1) = v(5);
stress = zeros(tf,1);
strain = zeros(tf,1);
time = zeros(tf,1);

AA1 = zeros(tf,1); AA2 = zeros(tf,1);AA3 = zeros(tf,1);
n = .3;
m = 20;
phi_h1 = .5;
phi_h2 = 10;
vol = [V1(1); V3(1); V5(1)];

phi_hs = zeros(tf,3); 

for t = 1:tf
    
    F0(4) = .2*t/tf;
    
    Fstar = FstarAnalytic();
    P0 = GetP0(Fstar);
        
    [dhpq, dhqp] = Gethdot(Fstar,P0,1,vol);
    
    phi_h = phi_h1*abs(h13(t))^n + phi_h2*abs(h13(t)*1.6/0.3750)^m;
    dh13 = -1/phi1*(dhpq + phi0 + phi_h);
    dh31 = -1/phi1*(dhqp - phi0 - phi_h);

    phi_hs(t+1,1) = phi_h;
    AA1(t+1) = dhqp;

    [dhpq, dhqp] = Gethdot(Fstar,P0,2,vol);
    
    phi_h = phi_h1*abs(h35(t))^n + phi_h2*abs(h35(t)*2.36/0.3750)^m;
    dh35 = -1/phi1*(dhpq + phi0 + phi_h);
    dh53 = -1/phi1*(dhqp - phi0 - phi_h);
    
    AA2(t+1) = dhqp;
    phi_hs(t+1,2) = phi_h;

    [dhpq, dhqp] = Gethdot(Fstar,P0,3,vol);
    
    phi_h = phi_h1*abs(h51(t))^n + phi_h2*abs(h51(t)*1.1/0.3750)^m;
    dh51 = -1/phi1*(dhpq + phi0 + phi_h);
    dh15 = -1/phi1*(dhqp - phi0 - phi_h);
    
    AA3(t+1) = dhqp;
    phi_hs(t+1,3) = phi_h;

    if(h13(t) < -3/8/a13 || h13(t) > 3/8/a13)
        dh13 = 0;
        dh31 = 0;
    end
    if(h35(t) < -3/8/a35 || h35(t) > 1/4/a35)
        dh35 = 0;
        dh53 = 0;
    end
    if(h51(t) < -1/4/a51 || h51(t) > 3/8/a51)
        dh51 = 0;
        dh15 = 0;
    end
    if(dh13 < 0)
        dh13 = 0;
    end
    if(dh31 > 0)
        dh31 = 0;
    elseif(dh31 < 0)
        dh13 = dh31;
    end
    if(dh35 < 0)
        dh35 = 0;
    end
    if(dh53 > 0)
        dh53 = 0;
    elseif(dh53 < 0)
        dh35 = dh53;
    end
    if(dh51 < 0)
        dh51 = 0;
    end
    if(dh15 > 0)
        dh15 = 0;
    elseif(dh15 < 0)
        dh51 = dh15;
    end
    if(v(2) < 0)
        C(:,:,2) = C(:,:,3);
        Cinv(:,:,2) = Cinv(:,:,3);
    elseif(v(2) > 0)
        C(:,:,2) = C(:,:,1);
        Cinv(:,:,2) = Cinv(:,:,1);
    else
        C(:,:,2) = zeros(9,9);
        Cinv(:,:,2) = zeros(9,9);
    end
    if(v(4) < 0)
        C(:,:,4) = C(:,:,5);
        Cinv(:,:,4) = Cinv(:,:,5);
    elseif(v(4) > 0)
        C(:,:,4) = C(:,:,3);
        Cinv(:,:,4) = Cinv(:,:,3);
    else
        C(:,:,4) = zeros(9,9);
        Cinv(:,:,4) = zeros(9,9);
    end
    if(v(6) < 0)
        C(:,:,6) = C(:,:,1);
        Cinv(:,:,6) = Cinv(:,:,1);
    elseif(v(6) > 0)
        C(:,:,6) = C(:,:,5);
        Cinv(:,:,6) = Cinv(:,:,5);
    else
        C(:,:,6) = zeros(9,9);
        Cinv(:,:,6) = zeros(9,9);
    end 
    
%     h13(t+1) = h13(t) + dt*dh13;
%     h35(t+1) = h35(t) + dt*dh35;
%     h51(t+1) = h51(t) + dt*dh51;
%     v1dot = a13*dh13 - a51*dh51;
%     v3dot = -a13*dh13 + a35*dh35;
%     v5dot = -a35*dh35 + a51*dh51;
%     
%     vdot = [v1dot; v3dot; v5dot]*dt;
%     V1(t+1) = V1(t) + v1dot*dt;
%     V3(t+1) = V3(t) + v3dot*dt;
%     V5(t+1) = V5(t) + v5dot*dt;
%     vol = [V1(t); V3(t); V5(t)];
    
    [dv,dh] = UpdateVol(dh13,dh35,dh51,dt);
    vol = vol + dv;
    V1(t+1) = vol(1);
    V3(t+1) = vol(2);
    V5(t+1) = vol(3);
    h13(t+1) = h13(t) + dh(1);
    h35(t+1) = h35(t) + dh(2);
    h51(t+1) = h51(t) + dh(3);

    time(t+1) = t*dt;

    stress(t+1) = P0(4);
    strain(t+1) = F0(4);
end

figure(2)

hold on
plot(time,h13,'--','LineWidth',1.3)
plot(time,h35,'LineWidth',1.3)
plot(time,h51,'LineWidth',1.3)
xlabel('Time')
ylabel('Boundary displacment')
legend('h_{13}','h_{35}','h_{51}','Location','NorthWest')
grid on
hold off

figure(3)
hold on
plot(time,V1,'LineWidth',1.3)
plot(time,V3,'LineWidth',1.3)
plot(time,V5,'LineWidth',1.3)
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
plot(AA1)
plot(AA2)
plot(AA3)
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
        
function [vol, h] = UpdateVol(dh13,dh35,dh51,dt)
global N Vtot C Cinv v Fgb F0  a D Aa I

    h13 = dh13*dt;
    h35 = dh35*dt;
    h51 = dh51*dt;
    
    h = [h13; h35; h51]; ha = abs(h);
    Ah = full(adjacency(D,h));
    Aha = full(adjacency(D,ha));
    vtest = v; 
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
%     dV(2,1) = 1/2*(h13*a(1) + h51*a(3) + h13*(a(1)) + h35*(-a(2))) = h13*a13
%     dV(1,3) = 1/2*(h13*a(2) + h51*a(3) + h13*(a(2)) + h51*(-a(3))) = h13*a35
%     dV(3,2) = 1/2*(h51*a(1) + h35*a(3) + h51*(a(1)) + h35*(-a(3))) = h51*a13
    dV = diag(dV);
    counter = 1;
    for i = 1:2:N
       v(i) = v(i) + dV(counter);
       counter = counter + 1;
    end
    
    vol = I*h;
%     vtest(2) = vtest(2) + a(1)*h13;
%     vtest(4) = vtest(4) + a(2)*h35;
%     vtest(6) = vtest(6) + a(3)*h51;
% 
%     vtest(1) = vtest(1) - a(1)*(abs(h(1)) - h(1))/2 - a(3)*(abs(h(3)) + h(3))/2;
%     vtest(3) = vtest(3) - a(1)*(abs(h(1)) + h(1))/2 - a(2)*(abs(h(2)) - h(2))/2;
%     vtest(5) = vtest(5) - a(2)*(abs(h(2)) + h(2))/2 - a(3)*(abs(h(3)) - h(3))/2;
%     dv_test = [vtest(1);vtest(3);vtest(5)];
%   if(v(2) > 0)
%       v(3) = v(3) - a13*dh13*dt;
%   elseif(v(2) < 0)
%       v(1) = v(1) + a13*dh13*dt;
%   end
%   if(v(4) > 0)
%       v(5) = v(5) - a35*dh35*dt;
%   elseif(v(4) < 0)
%       v(3) = v(3) + a35*dh35*dt;
%   end
%   if(v(6) > 0)
%       v(1) = v(1) - a51*dh51*dt;
%   elseif(v(6) < 0)
%       v(5) = v(5) + a51*dh51*dt;
%   end
%       
  for i = 1:2:N
      
%       if(abs(v(i)) > 1/3)
%           v(i) = 1/3*sign(v(i));
%       end
      if(v(i) < 0)
         v(i) = 0;
      end
  end
end

function [dhpq, dhqp] = Gethdot(Fstr,P0,i,vol)
    global N Vtot C Cinv  Fgb F0 a12 xi_inv edges alpha
   
    Fstardh_ana = zeros(3,3,N-edges);
    dhpq = 0; dhqp = 0; I = [1;0;0;0;1;0;0;0;1];
    % ****** vm < 0 ******
    pq = i*2;
    counter = 1;
    
    for pp = 1:edges
            
        m = pp*2; q = m - 1; p = m + 1;
        
        if(p > N)
            p = 1;
        end
        
        C_diff = eye(9) - Cinv(:,:,p)*C(:,:,q);
%         alpha = Vtot*F0 - v(pq)*Fgb(:,:,pq) - v(p)*Fgb(:,:,p) + (v(p)*Cinv(:,:,p)*C(:,:,q))*Fgb(:,:,q);
        dfstar_dhpq = -xi_inv(:,:,q)^2*(C_diff)*alpha(:,:,q) + xi_inv(:,:,q)*(I-Fgb(:,:,pq) + C_diff*Fgb(:,:,q));

        Fstardh_ana(:,:,counter) = reshape(dfstar_dhpq,[3,3]);
        counter = counter + 1;
            
    end

    if(mod(N,2) == 1)
        
        C_diff = eye(9) - Cinv(:,:,q)*C(:,:,p);
%         alpha = Vtot*F0 - v(pq)*Fgb(:,:,pq) - v(q)*Fgb(:,:,q) + (v(q)*Cinv(:,:,q)*C(:,:,p))*Fgb(:,:,p);
        dfstar_dhpq = -xi_inv(:,:,p)^2*(C_diff)*alpha(:,:,p) + xi_inv(:,:,p)*(I - Fgb(:,:,pq) + C_diff*Fgb(:,:,p));

        Fstardh_ana(:,:,end) = reshape(dfstar_dhpq,[3,3]);

    end

%     if(vm < 0)
%         
%         for pp = 1:edges
%             
%             m = pp*2; q = m - 1; p = m + 1;
%             
%             alpha = Vtot*F0 - v(pq)*Fgb(:,:,pq) - v(p)*Fgb(:,:,p) + (abs(v(pq))*Cinv(:,:,pq)*C(:,:,q)+v(p)*Cinv(:,:,p)*C(:,:,q))*Fgb(:,:,q);
%             dfstar_dhpq = -xi_inv(:,:,q)^2*(eye(9) - Cinv(:,:,pq)*C(:,:,q))*alpha + xi_inv(:,:,q)*(-Fgb(:,:,pq) + Cinv(:,:,pq)*C(:,:,q)*Fgb(:,:,q));
%             
%             Fstardh_ana(:,:,counter) = reshape(dfstar_dhpq,[3,3]);
%             counter = counter + 1;
%             
%         end
%         
%         if(mod(N,2) == 1)
%             
%             alpha = Vtot*F0 - v(pq)*Fgb(:,:,pq) - v(q)*Fgb(:,:,q) + (abs(v(pq))*Cinv(:,:,pq)*C(:,:,p) + v(q)*Cinv(:,:,q)*C(:,:,p))*Fgb(:,:,p);
%             dfstar_dhpq = -xi_inv(:,:,p)^2*(eye(9) - Cinv(:,:,q)*C(:,:,p))*alpha + xi_inv(:,:,p)*(-Fgb(:,:,pq) + Cinv(:,:,pq)*C(:,:,p)*Fgb(:,:,p));
%             
%             Fstardh_ana(:,:,end) = reshape(dfstar_dhpq,[3,3]);
%             
%         end
%         % for dfdh fstar2/q -- v1 const --
%        
% %         alpha = Vtot*F0 - v(m)*Fgb(:,:,m)-v(p)*Fgb(:,:,p) + (v(m)*Cinv(:,:,m)*C(:,:,q)+v(p)*Cinv(:,:,p)*C(:,:,q))*Fgb(:,:,q);
% % 
% %         dfstar_dhpq = -xi_inv(:,:,q)^2*(eye(9) - Cinv(:,:,m)*C(:,:,q))*alpha + xi_inv(:,:,q)*(-Fgb(:,:,m) + Cinv(:,:,m)*C(:,:,q)*Fgb(:,:,q));
% % 
% %         Fstardh_q_ana = reshape(dfstar_dhpq,[3,3]);
% %        
% %         % for dfdh fstar2/pq
% % %         xi_inv = (v(1)*Cinv(:,:,1)*C(:,:,2) + v(2)*eye(6) + v(3)*Cinv(:,:,3)*C(:,:,2))\eye(6);
% % %         alpha = Vtot*F0 - v(3)*Fgb(:,:,3)-v(1)*Fgb(:,:,1) + (v(3)*Cinv(:,:,3)*C(:,:,2)+v(1)*Cinv(:,:,1)*C(:,:,2))*Fgb(:,:,2);
% % %
% % %         dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,3)*C(:,:,2))*alpha + xi_inv*(Fgb(:,:,3) - Cinv(:,:,3)*C(:,:,2)*Fgb(:,:,2));
% % %
% % %         Fstardh_pq_ana = unpackF(dfstar_dhpq);
% % 
% %         % dfdh fstar 3/p
% %    
% %         alpha = Vtot*F0 - v(m)*Fgb(:,:,m) - v(q)*Fgb(:,:,q) + (v(m)*Cinv(:,:,m)*C(:,:,p) + v(q)*Cinv(:,:,q)*C(:,:,p))*Fgb(:,:,p);
% % 
% %         dfstar_dhpq = xi_inv(:,:,p)^2*(-eye(9) + Cinv(:,:,q)*C(:,:,p))*alpha + xi_inv(:,:,p)*(-Fgb(:,:,m) + Cinv(:,:,m)*C(:,:,p)*Fgb(:,:,p));
% % 
% %         Fstardh_p_ana = reshape(dfstar_dhpq,[3,3]);
%    
%    
%     % ***** vm > 0 *****
%     elseif(vm >= 0)
%        
%         
%         for pp = 1:edges
%             
%             m = pp*2; q = m - 1; p = m + 1;
%             
%             alpha = Vtot*F0 - v(pq)*Fgb(:,:,pq) - v(p)*Fgb(:,:,p) + (abs(v(pq))*Cinv(:,:,pq)*C(:,:,q) + v(p)*Cinv(:,:,p)*C(:,:,q))*Fgb(:,:,q);
%             dfstar_dhpq = -xi_inv(:,:,q)^2*(eye(9) - Cinv(:,:,pq)*C(:,:,q))*alpha + xi_inv(:,:,q)*(-Fgb(:,:,pq) + Cinv(:,:,pq)*C(:,:,q)*Fgb(:,:,q));
%             
%             Fstardh_ana(:,:,counter) = reshape(dfstar_dhpq,[3,3]);
%             counter = counter + 1;
%             
%         end
%         
%         if(mod(N,2) == 1)
%             
%             alpha = Vtot*F0 - v(pq)*Fgb(:,:,pq) - v(q)*Fgb(:,:,q) + (abs(v(pq))*Cinv(:,:,pq)*C(:,:,p) + v(q)*Cinv(:,:,q)*C(:,:,p))*Fgb(:,:,p);
%             dfstar_dhpq = -xi_inv(:,:,p)^2*(eye(9) - Cinv(:,:,q)*C(:,:,p))*alpha + xi_inv(:,:,p)*(-Fgb(:,:,pq) + Cinv(:,:,pq)*C(:,:,p)*Fgb(:,:,p));
%             
%             Fstardh_ana(:,:,end) = reshape(dfstar_dhpq,[3,3]);
%             
%         end
%         
%         % dfdh fstar 1/q -- v3 const --
%         
% %         alpha = Vtot*F0 - v(m)*Fgb(:,:,m) - v(p)*Fgb(:,:,p) + (v(m)*Cinv(:,:,m)*C(:,:,q) + v(p)*Cinv(:,:,p)*C(:,:,q))*Fgb(:,:,q);
% % 
% %         dfstar_dhpq = -xi_inv(:,:,q)^2*(eye(9) - Cinv(:,:,m)*C(:,:,q))*alpha + xi_inv(:,:,q)*(-Fgb(:,:,m) + Cinv(:,:,m)*C(:,:,q)*Fgb(:,:,q));
% % 
% %         Fstardh_q_ana = reshape(dfstar_dhpq,[3,3]);
% % 
% % %         % for dfdh fstar2/pq
% % %         xi_inv = (v(1)*Cinv(:,:,1)*C(:,:,2) + v(2)*eye(6) + v(3)*Cinv(:,:,3)*C(:,:,2))\eye(6);
% % %         alpha = Vtot*F0 - v(3)*Fgb(:,:,3)-v(1)*Fgb(:,:,1) + (v(3)*Cinv(:,:,3)*C(:,:,2)+v(1)*Cinv(:,:,1)*C(:,:,2))*Fgb(:,:,2);
% % %
% % %         dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,1)*C(:,:,2))*alpha + xi_inv*(Fgb(:,:,1) - Cinv(:,:,1)*C(:,:,2)*Fgb(:,:,2));
% % %
% % %         Fstardh_pq_ana = unpackF(dfstar_dhpq);
% %        
% %         % dfdh fstar 3/p
% %         
% %         alpha = Vtot*F0 - v(m)*Fgb(:,:,m) - v(q)*Fgb(:,:,q) + (v(m)*Cinv(:,:,m)*C(:,:,p) + v(q)*Cinv(:,:,q)*C(:,:,p))*Fgb(:,:,p);
% % 
% %         dfstar_dhpq = -xi_inv(:,:,p)^2*(eye(9) - Cinv(:,:,q)*C(:,:,p))*alpha + xi_inv(:,:,p)*(-Fgb(:,:,m) + Cinv(:,:,m)*C(:,:,p)*Fgb(:,:,p));
% % 
% %         Fstardh_p_ana = reshape(dfstar_dhpq,[3,3]);
%        
%     else % vm = 0
%        
%     end    
%     V_e = zeros(1,N);
%     for j = 1:edges
%         m = j*2; q = m - 1; p = m + 1;
%         if(v(m) < 0)
%             V_e(p) = V_e(p) + abs(v(m));
%         else%(v(m) > 0)
%             V_e(q) = V_e(q) + abs(v(m));
%         end
%     end
%     for j = 1:N
%         
%         if(mod(j,2) == 0)
%             V_e(j) = [];
%         end
%     end
    counter = 1; Fstardh = zeros(3,3);
    for j = 1:2:(N-edges)*2
        
        Fstardh = Fstardh + Fstardh_ana(:,:,counter)*vol(counter);
        counter = counter + 1;
    end
   
%    [Fstardh_p,Fstardh_q,Fstardh_pq] = dFstr_num(Fstr,vm);
   
%    t1 = Fstardh_q_ana - Fstardh_q
%    t2 = Fstardh_pq_ana - Fstardh_pq
%    t3 = Fstardh_p_ana - Fstardh_p
%     if(vm < 0)
%         vp = v(p); vq = v(q) + abs(v(m)); % df_q = df_pq
%     elseif(vm > 0)
%         vp = v(p)+ abs(v(m)); vq = v(q); % df_p = df_pq
%     else % vm = 0
%         vp = v(p); vq = v(q);
%     end

    m = i*2; q = m - 1; p = m + 1;
    if(p > N)
        p = 1;
    end
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
    global N Vtot C Cinv v Fgb F0 xi_inv alpha
   
    Fstr = zeros(9,1,N);
    
    for i = 1:N
        xi = v(i)*eye(9);
        alpha_temp = Vtot*F0;
        for j = 1:N            
            if(i == j)
                continue
            end        
            xi = xi + abs(v(j))*Cinv(:,:,j)*C(:,:,i);
            alpha_temp = alpha_temp - v(j)*Fgb(:,:,j);
        end
        alpha(:,:,i) = alpha_temp + (xi-v(i)*eye(9))*Fgb(:,:,i);
        if(v(i) == 0)
            Fstr(:,:,i) = Fgb(:,:,i);
            xi_inv(:,:,i) = eye(9);
        else
            xi_inv(:,:,i) = xi\eye(9);
            Fstr(:,:,i) = xi_inv(:,:,i)*alpha(:,:,i);
        end
       
    end
   
end
function [C,Cinv] = SetC(th)

    C11 = 169.3097; C12 = 122.5; C44 = 76;
%     C11 = 169.3097; C12 = 87.2201 ; C44 = 41.0448;

    S11 = (C11+C12)/((C11-C12)*(C11+2*C12));
    S12 = -C12/((C11-C12)*(C11+2*C12));
    S44 = 1/C44;
    C0 = C11 - C12 -2*C44;
    S0 = S11-S12-1/2*S44;
   
    C = zeros(9,9);
    Cinv = zeros(9,9);
   
    e1 = [cos(th);sin(th);0];
    e2 = [-sin(th);cos(th);0];
    e3 = [0;0;1];
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
function [CC,CCinv,v,Fgb] = init(angle)
    global N Vtot
   
    CC = zeros(9,9,N);
    CCinv = CC; v = zeros(N,1); Fgb = zeros(3,3,N);
     
    for i = 1:N
       
      [CC(:,:,i),CCinv(:,:,i)] = SetC(angle(i));
        
      v(i) = Vtot/(N-2);
      Fgb(:,:,i) = eye(3);
     
      if( mod(i,2) == 0 )
          CC(:,:,i) = zeros(9,9);
          CCinv(:,:,i) = zeros(9,9);
          Fgb(1,2,i) = 2*tan(angle(i-1)/2);
          v(i) = 0;
      end
    end

end
