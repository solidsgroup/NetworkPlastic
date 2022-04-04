clear all
global N Vtot C Cinv v Fgb F0  a13 a35 xi_inv edges alpha

N = 5; edges = 2;
Vtot = 1;

F0 = [1,.1 ,0;
      0, 1 ,0;
      0, 0 ,1];
Fgb = zeros(9,1,N);
xi_inv = zeros(9,9,N);
alpha = zeros(9,1,N);
C_avg = zeros(9,9);
angle = [0.2 0.2 0.3 0.3 .4];
% angle = [0.2, 0.2, .3];
[C,Cinv,v,Fgbt] = init(angle);
F0 = reshape(F0, [9,1]);
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

dt = 1e-5;
tf = floor(1/dt);
a13 = 1/3; a35 = 1/3;
dh13 = 0; dh31 = 0; % first boundary
dh35 = 0; dh53 = 0; % second boundary
phi1 = .006;
phi0 = .2;

h13 = zeros(tf,1); h35 = zeros(tf,1);
V1 = zeros(tf,1); V1(1) = v(1);
V3 = zeros(tf,1); V3(1) = v(3);
V5 = zeros(tf,1); V5(1) = v(5);
stress = zeros(tf,1);
strain = zeros(tf,1);
time = zeros(tf,1);

vm13 = 0;
vm35 = 0;

AA1 = zeros(tf,1);
AA2 = zeros(tf,1);
vol = [V1(1); V3(1); V5(1)];
for t = 1:tf
    
    F0(4) = .25*t/tf;
    
    Fstar = FstarAnalytic();
    P0 = GetP0(Fstar);
    
    [dhpq, dhqp] = Gethdot(Fstar,P0,1,vol);
    
    dh13 = -1/phi1*(dhpq + phi0);
    dh31 = -1/phi1*(dhqp - phi0);
    
    AA1(t+1) = dhqp;
    [dhpq, dhqp] = Gethdot(Fstar,P0,2,vol);
    
    dh35 = -1/phi1*(dhpq + phi0);
    dh53 = -1/phi1*(dhqp - phi0);
    
    AA2(t+1) = dhqp;
    if(abs(v(2)) >= 1/3)
        dh13 = 0;
        dh31 = 0;
    end
    if(abs(v(4)) >= 1/3)
        dh35 = 0;
        dh53 = 0;
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
    if(v(2) < 0)
        C(:,:,2) = C(:,:,1);
        Cinv(:,:,2) = Cinv(:,:,1);
    elseif(v(2) > 0)
        C(:,:,2) = C(:,:,3);
        Cinv(:,:,2) = Cinv(:,:,3);
    else
        C(:,:,2) = eye(9);
        Cinv(:,:,2) = eye(9);
    end
    if(v(4) < 0)
        C(:,:,4) = C(:,:,3);
        Cinv(:,:,4) = Cinv(:,:,3);
    elseif(v(4) > 0)
        C(:,:,4) = C(:,:,5);
        Cinv(:,:,4) = Cinv(:,:,5);
    else
        C(:,:,4) = eye(9);
        Cinv(:,:,4) = eye(9);
    end
    
    %%%%%% calc Astar
%     GetAstar(dh13,dh35,dt);
    

    h13(t+1) = h13(t) + dt*dh13;
    h35(t+1) = h35(t) + dt*dh35;
    v1dot = a13*dh13;
    v3dot = -a13*dh13 + a35*dh35;
    v5dot = -a35*dh35;
    
    V1(t+1) = V1(t) + v1dot*dt;
    V3(t+1) = V3(t) + v3dot*dt;
    V5(t+1) = V5(t) + v5dot*dt;
    vol = [V1(t); V3(t); V5(t)];
    vm13 = a13*h13(t);
    vm35 = a35*h35(t);
    
%     v(1) = V1(t); v(2) = vm13; v(3) = 1/3 + vm35; v(4) = vm35; v(5) = V5(1);
    
    UpdateVol(dh13,dh35,dt);
      
%     if(vm13 < 0)
%         v(1) = V1(t); v(2) = vm13;
%     elseif(vm13 > 0)
%         v(1) = V1(1); v(2) = vm13;
%     end
%     if(vm35 < 0)
%         v(5) = V5(1); v(4) = vm35;
%     elseif(vm35 > 0)
%         v(5) = V5(t); v(4) = vm35;
%     end
%     if(dh13 > 0 && dh35 > 0)
%         v(3) = v(3) - a13*dh13*dt;
%     elseif(dh13 < 0 && dh35 > 0)
%         v(3) = v(3);
%     elseif(dh13 > 0 && dh35 < 0)
%         v(3) = v(3) - (a13*dh13 - a35*dh35)*dt;
%     elseif(dh13 < 0 && dh35 < 0)
%         v(3) = v(3) + a35*dh35;
%     end

%     if(vm13 > 0 && vm35 > 0)
%         v(3) = V3(t) - vm13;
%     elseif(vm13 < 0 && vm35 < 0)
%         v(3) = V3(t) - vm35;
%     end
    time(t+1) = t*dt;

    temp = C(:,:,1)*(Fstar(:,:,1) - Fgb(:,:,1));
    stress(t+1) = temp(4);
    strain(t+1) = F0(4);
end

figure(2)

hold on
plot(time,h13,'--','LineWidth',1.3)
plot(time,h35,'LineWidth',1.3)
xlabel('Time')
ylabel('Boundary displacment')
legend('h_{13}','h_{35}')
grid on
hold off

figure(3)
hold on
plot(time,V1,'LineWidth',1.3)
plot(time,V3,'LineWidth',1.3)
plot(time,V5,'LineWidth',1.3)
xlabel('Time')
ylabel('Volume')
legend('V_1','V_3','V_5','Location','SouthWest')
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
plot([1 tf],[1 1]*phi0)
plot([1 tf],-[1 1]*phi0)
title('dA^*/dh_{pq}')
hold off

function GetAstar(dh13,dh35,dt)
global N Vtot C Cinv v Fgb F0  a13 a35
    
    v_save = v; A = zeros(2,1);
    
    UpdateVol(dh13,0,dt);
    Fstr = FstarAnalytic();
    
    for i = 1:N
        P = C(:,:,i)*(Fstr(:,:,i)-Fgb(:,:,i));
        P = reshape(P,[3,3]);
        Fstar = reshape(Fstr(:,:,i),[3,3]);
        FGB = reshape(Fgb(:,:,i),[3,3]);
        for j = 1:3
            for k = 1:3
                A(1) = A(1) + 1/2*(Fstar(j,k)-FGB(j,k))*P(j,k);
            end
        end
    end
    v = v_save;
    UpdateVol(0,dh35,dt);
    Fstr = FstarAnalytic();
    for i = 1:N
        P = C(:,:,i)*(Fstr(:,:,i)-Fgb(:,:,i));
        P = reshape(P,[3,3]);
        Fstar = reshape(Fstr(:,:,i),[3,3]);
        FGB = reshape(Fgb(:,:,i),[3,3]);
        for j = 1:3
            for k = 1:3
                A(2) = A(2) + 1/2*(Fstar(j,k)-FGB(j,k))*P(j,k);
            end
        end
    end
    v = v_save;
    
end
        
   
function UpdateVol(dh13,dh35,dt)
global N Vtot C Cinv v Fgb F0  a13 a35

    h13 = dh13*dt;
    h35 = dh35*dt;
    v(2) = v(2) + a13*h13;
    v(4) = v(4) + a35*h35;
    
    v(3) = v(3) - a13*(abs(h13) + h13)/2 - a35*(abs(h35) - h35)/2;
    v(1) = v(1) - a13*(abs(h13) - h13)/2;
    v(5) = v(5) - a35*(abs(h35) + h35)/2;
    
    
%   v(2) = v(2) + a13*dh13*dt;
%   v(4) = v(4) + a35*dh35*dt;
%   
%   if(v(2) > 0)
%       v(1) = 1/3;
%       v(3) = v(3) - a13*dh13*dt;
%   elseif(v(2) < 0)
%       v(1) = v(1) + a13*dh13*dt;
%       v(3) = v(3);
%   end
%   if(v(4) > 0)
%       v(3) = v(3);
%       v(5) = v(5) - a35*dh35*dt;
%   elseif(v(4) < 0)
%       v(3) = v(3) + a35*dh35*dt;
%       v(5) = 1/3;
%   end
%   for i = 1:N
%       
%       if(abs(v(i)) > 1/3)
%           v(i) = 1/3*sign(v(i));
%       end
%       if(mod(i,2) == 1 && v(i) < 0)
%          v(i) = 0;
%       end
%   end
end

function [dhpq, dhqp] = Gethdot(Fstr,P0,i,vol)
    global N Vtot C Cinv v Fgb F0 a12 xi_inv edges alpha
   
    Fstardh_ana = zeros(3,3,N-edges);
    dhpq = 0; dhqp = 0; I = [1;0;0;0;1;0;0;0;1];
    % ****** vm < 0 ******
    pq = i*2;
    counter = 1;
    
    for pp = 1:edges
            
        m = pp*2; q = m - 1; p = m + 1;
        
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
    Fstr_test = zeros(9,1,N);
   
     if(N == 2)
     % N = 2
         temp = (v(1)*eye(6) + v(3)*Cinv(:,:,3)*C(:,:,1))\(Vtot*F0 + v(3)*(-Fgb(:,:,3) + Cinv(:,:,3)*C(:,:,1)*Fgb(:,:,1)) );
         Fstr_test(:,1,1) = temp;

         temp = (v(3)*eye(6) + v(1)*Cinv(:,:,1)*C(:,:,3))\(Vtot*F0 + v(1)*(-Fgb(:,:,1) + Cinv(:,:,1)*C(:,:,3)*Fgb(:,:,3)) );
         Fstr_test(:,1,3) = temp;
     elseif(N == 3)
     % N = 3
         temp = (v(1)*eye(9) + abs(v(2))*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))\(Vtot*F0 - v(2)*Fgb(:,:,2) - v(3)*Fgb(:,:,3) + ...
            (abs(v(2))*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))*Fgb(:,:,1) );
        Fstr_test(:,:,1) = temp;
       
        temp = (abs(v(2))*eye(9) + v(1)*Cinv(:,:,1)*C(:,:,2) + v(3)*Cinv(:,:,3)*C(:,:,2))\(Vtot*F0 - v(1)*Fgb(:,:,1) - v(3)*Fgb(:,:,3) +...
                (v(1)*Cinv(:,:,1)*C(:,:,2) + v(3)*Cinv(:,:,3)*C(:,:,2))*Fgb(:,:,2) );
        Fstr_test(:,:,2) = temp;

        temp = (v(3)*eye(9) + v(1)*Cinv(:,:,1)*C(:,:,3) + abs(v(2))*Cinv(:,:,2)*C(:,:,3))\(Vtot*F0 - v(1)*Fgb(:,:,1) - v(2)*Fgb(:,:,2) +...
            (v(1)*Cinv(:,:,1)*C(:,:,3) + abs(v(2))*Cinv(:,:,2)*C(:,:,3))*Fgb(:,:,3) );
        Fstr_test(:,:,3) = temp;
     end
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
          CC(:,:,i) = eye(9);
          CCinv(:,:,i) = eye(9);
          Fgb(1,2,i) = 2*tan(angle(i-1)/2);
          v(i) = 0;
      end
    end

end