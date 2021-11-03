clear all
global N Vtot C Cinv v Fgb F0 a12

N = 3;
Vtot = 1;

F0 = [1,.1 ,0;
      0, 1 ,0;
      0, 0 ,1];
Fgb = zeros(6,1,N);
% angle = [-0.4118977,-0.4118977,0.4118977];
angle = [-pi/6,-pi/6,pi/6];
[C,Cinv,v,Fgbt] = init(angle);
% v(1) = 0.25; v(2) = 1-v(1);
F0 = pack(F0);
for i = 1:N
   Fgb(:,1,i) = pack(Fgbt(:,:,i)); % put into viogt notation
end
% dv = 0.001;
% vm = -.5:dv:.5;
% % vm = .49*sin(.1:dv:2*pi)+.5;
vm = .2;
v(1) = .5; v(2) = vm; v(3) = Vtot-(vm+v(3));
%v(2) = vm; v(1) = Vtot-(vm+0.5);
Fstr = FstarAnalytic(vm);
for i = 1:N
   F = unpackF(Fstr(:,:,i)); 
end

% for i = 1:lengt;h(vm)
%     A1 = 0;
%     A2 = 0;
%     A3 = 0;
%     if(vm(i) < 0)
%              F0(6) = (.25-.25/floor(length(vm)/2)*i)*(-1)^j;
%         v(1) = .6; v(2) = abs(vm(i)); v(3) = Vtot-(abs(vm(i))+v(1));
%         C(:,:,2) = C(:,:,1);
%         Cinv(:,:,2) = Cinv(:,:,1);
%     elseif(vm(i) > 0)
%              F0(6) = (.25/floor(length(vm)/2)*(i-floor(length(vm)/2)))*(-1)^j;
%         v(3) = .4; v(2) = abs(vm(i)); v(1) = Vtot-(abs(vm(i))+v(3));
%         C(:,:,2) = C(:,:,3);
%         Cinv(:,:,2) = Cinv(:,:,3);
%     else
%         C(:,:,2) = zeros(6,6);
%         Cinv(:,:,2) = zeros(6,6);
%     end
%     Fstr = FstarAnalytic(vm(i));
%     [da1,da2] = Gethdot(Fstr,vm(i));
% 
%     P1 = v(1)*C(:,:,1)*(Fstr(:,:,1)-Fgb(:,:,1));
%     P2 = v(2)*C(:,:,2)*(Fstr(:,:,2)-Fgb(:,:,2));
%     P3 = v(3)*C(:,:,3)*(Fstr(:,:,3)-Fgb(:,:,3));
%     P1 = unpack(P1); P2 = unpack(P2); P3 = unpack(P3);
%     F1 = Fstr(:,:,1)-Fgb(:,:,1);  F2 = Fstr(:,:,2)-Fgb(:,:,2);  F3 = Fstr(:,:,3)-Fgb(:,:,3); 
%     F1 = unpackF(F1); F2 = unpackF(F2); F3 = unpackF(F3);
% 
%     for m = 1:3
%         for n = 1:3
%             A1 = A1 + 1/2*P1(m,n)*F1(m,n);
%             A2 = A2 + 1/2*P2(m,n)*F2(m,n);
%             A3 = A3 + 1/2*P3(m,n)*F3(m,n);
%         end
%     end
%     
%     AA1(i) = A1;
%     AA2(i) = A2;
%     AA3(i) = A3;
%     
%     da1_ana(i) = da1;
%     da2_ana(i) = da2;
% %         F(i) = F0(6);
% end
% 
% dA1dh_calc = diff(AA1)./diff(vm);
% dA2dh_calc = diff(AA2)./diff(vm);
% dA3dh_calc = diff(AA3)./diff(vm);
% figure(1)
% hold on 
% plot(dA1dh_calc,'k-')
% plot(dA2dh_calc,'r--')
% plot(dA3dh_calc,'r-')
% plot(da1_ana,'k-.')
% plot(da2_ana,'r-.')
% ylim([-20 20])
% hold off
% A1_ana = cumsum(da1_ana.*dv);
% A2_ana = cumsum(-da2_ana.*dv);
% figure(2)
% hold on
% plot(AA1,'k-')
% plot(AA3,'r-')
% plot(AA2,'b-')
% plot(A1_ana,'k-.')
% plot(A2_ana,'r-.')
% hold off
% %
% Gethdot(Fstr,-vm)
%
% C(:,:,1)*(Fstr(:,:,1)-Fgb(:,:,1)) + C(:,:,2)*(Fstr(:,:,2)-Fgb(:,:,2)) - 2*C(:,:,3)*(Fstr(:,:,3)-Fgb(:,:,3))
% v(1)*Fstr(:,:,1) + v(2)*Fstr(:,:,2) + v(3)*Fstr(:,:,3) - Vtot*F0
%% bicrystal

v(1) = 0.5; v(2) = 0; v(3) = Vtot - v(1);

dt = 1e-5;
a12 = .5;
tf = floor(1/dt);
dh12 = 0;
dh21 = 0;
phi1 = .006;
phi0 = 0.2;
Fgb12 = [1, .4, 0;
         0, 1, 0;
         0, 0, 1];

I = [1;1;1;0;0;0];
h = zeros(tf,1);
V1 = zeros(tf,1); V1(1) = v(1);
V2 = zeros(tf,1); V2(1) = v(3);
stress = zeros(tf,1);
strain = zeros(tf,1);
time = zeros(tf,1);
AA1 = zeros(tf,1);
vm = 0;

for t = 1:tf

    F0(6) = -.2/tf*t;
    % F0(2) = 1-.1/tf*t;

    Fstar = FstarAnalytic(vm);
    
    [dhpq, dhqp] = Gethdot(Fstar,vm);
    
    dh12 = -1/phi1*(dhpq + phi0);
    dh21 = -1/phi1*(dhqp - phi0);

    AA1(t+1) = dhpq;

    if(h(t) <-1/2/a12 || h(t)> 1/2/a12)
        dh12 = 0;
        dh21 = 0;
    end
    if(dh12 < 0)
        dh12 = 0;
    end
    if(dh21 > 0)
        dh21 = 0;
    elseif(dh21 < 0)
        dh12 = dh21;
    end

    if(h(t) < 0)
        C(:,:,2) = C(:,:,1);
        Cinv(:,:,2) = Cinv(:,:,1);
    elseif(h(t) > 0)
        C(:,:,2) = C(:,:,3);
        Cinv(:,:,2) = Cinv(:,:,3);
    else
        C(:,:,2) = zeros(6,6);
        Cinv(:,:,2) = zeros(6,6);
    end

    h(t+1) = h(t) - dt*dh12; % negative due to graph orientation
    v1dot = a12*dh12;
    v2dot = -a12*dh12;
    vm = -V2(1)+V2(t);

    V1(t+1) = V1(t) + v1dot*dt;
    V2(t+1) = V2(t) + v2dot*dt;

    if(vm < 0)
        v(1) = V1(1); v(3) = V2(t); v(2) = -vm;
    elseif(vm > 0)
        v(1) = V1(t); v(3) = V2(1); v(2) = -vm;
    end
    time(t+1) = t*dt;

    temp = C(:,:,1)*(Fstar(:,:,1) - Fgb(:,:,1));
    stress(t+1) = temp(6);
    strain(t+1) = F0(6)*2;
end


figure(2)
yyaxis left
plot(time,h)
xlabel('Time')
ylabel('Boundary displacment')

yyaxis right
plot(time,strain,'--')
ylabel('F_{12}')

figure(3)
hold on
plot(time,V1)
plot(time,V2)
xlabel('Time')
ylabel('Volume')
legend('V_1','V_2')
hold off
figure(4)
plot(strain,stress)
xlabel('Strain')
ylabel('Stress')
figure(5)
hold on
plot(AA1)
plot([1 tf],[1 1]*phi0)
plot([1 tf],-[1 1]*phi0)
title('dA^*/dh_{pq}')
hold off

function Fstr = FstarAnalytic(vm)
    global N Vtot C Cinv v Fgb F0
   
    Fstr = zeros(6,1,N);
   
%     if(vm == 0)
%         % N = 2
%         temp = (v(1)*eye(6) + v(3)*Cinv(:,:,3)*C(:,:,1))\(Vtot*F0 + v(3)*(-Fgb(:,:,3) + Cinv(:,:,3)*C(:,:,1)*Fgb(:,:,1)) );
%         Fstr(:,1,1) = temp;
%        
%         Fstr(:,1,2) = [1;1;1;0;0;0];
%        
%         temp = (v(3)*eye(6) + v(1)*Cinv(:,:,1)*C(:,:,3))\(Vtot*F0 + v(1)*(-Fgb(:,:,1) + Cinv(:,:,1)*C(:,:,3)*Fgb(:,:,3)) );
%         Fstr(:,1,3) = temp;
    if(N == 3)
        temp = (v(1)*eye(6) + abs(v(2))*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))\(Vtot*F0 - v(2)*Fgb(:,:,2) - v(3)*Fgb(:,:,3) +...
            (v(2)*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))*Fgb(:,:,1) );
        Fstr(:,1,1) = temp;
        
        if(vm ~= 0)
            temp = (abs(v(2))*eye(6) + v(1)*Cinv(:,:,1)*C(:,:,2) + v(3)*Cinv(:,:,3)*C(:,:,2))\(Vtot*F0 - v(1)*Fgb(:,:,1) - v(3)*Fgb(:,:,3) +...
                (v(1)*Cinv(:,:,1)*C(:,:,2) + v(3)*Cinv(:,:,3)*C(:,:,2))*Fgb(:,:,2) );
            Fstr(:,1,2) = temp;
        else 
            Fstr(:,1,2) = [1;1;1;0;0;0];
        end
       
        temp = (v(3)*eye(6) + v(1)*Cinv(:,:,1)*C(:,:,3) + abs(v(2))*Cinv(:,:,2)*C(:,:,3))\(Vtot*F0 - v(1)*Fgb(:,:,1) - v(2)*Fgb(:,:,2) +...
            (v(1)*Cinv(:,:,1)*C(:,:,3) + v(2)*Cinv(:,:,2)*C(:,:,3))*Fgb(:,:,3) );
        Fstr(:,1,3) = temp;
    end
end
   
function ff = pack(F)

    ff = zeros(6,1);
    ff(1) = F(1,1); ff(2) = F(2,2); ff(3) = F(3,3);
    ff(4) = F(2,3)/2; ff(5) = F(1,3)/2; ff(6) = F(1,2)/2;
end
function mat = unpackF(F)

    mat = zeros(3,3);
    mat(1,1) = F(1); mat(2,2) = F(2); mat(3,3) = F(3);
    mat(2,3) = F(4)*2; mat(1,3) = F(5)*2; mat(1,2) = F(6)*2;
end
function mat = unpack(F)

    mat = zeros(3,3);
    mat(1,1) = F(1); mat(2,2) = F(2); mat(3,3) = F(3);
    mat(2,3) = F(4); mat(1,3) = F(5); mat(1,2) = F(6);
end
function [dhpq, dhqp] = Gethdot(Fstr,vm)
    global N Vtot C Cinv v Fgb F0 a12
   
    dhpq = 0; dhqp = 0;
   
    % ****** vm < 0 ******
    if(vm < 0)
       
        vp = v(3); vq = v(1) + v(2);
       
        % for dfdh fstar2/q -- v1 const --
       
        xi_inv = (v(1)*eye(6) + abs(v(2))*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))\eye(6);
        alpha = Vtot*F0 - v(2)*Fgb(:,:,2)-v(3)*Fgb(:,:,3) + (v(2)*Cinv(:,:,2)*C(:,:,1)+v(3)*Cinv(:,:,3)*C(:,:,1))*Fgb(:,:,1);

        dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,3)*C(:,:,1))*alpha + xi_inv*(-Fgb(:,:,2) + Cinv(:,:,3)*C(:,:,1)*Fgb(:,:,1));

        Fstardh_q_ana = unpackF(dfstar_dhpq);
       
        % for dfdh fstar2/pq
%         xi_inv = (v(1)*Cinv(:,:,1)*C(:,:,2) + v(2)*eye(6) + v(3)*Cinv(:,:,3)*C(:,:,2))\eye(6);
%         alpha = Vtot*F0 - v(3)*Fgb(:,:,3)-v(1)*Fgb(:,:,1) + (v(3)*Cinv(:,:,3)*C(:,:,2)+v(1)*Cinv(:,:,1)*C(:,:,2))*Fgb(:,:,2);
% 
%         dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,3)*C(:,:,2))*alpha + xi_inv*(Fgb(:,:,3) - Cinv(:,:,3)*C(:,:,2)*Fgb(:,:,2));
% 
%         Fstardh_pq_ana = unpackF(dfstar_dhpq);

        % dfdh fstar 3/p
        xi_inv = (v(1)*Cinv(:,:,1)*C(:,:,3)+ abs(v(2))*Cinv(:,:,2)*C(:,:,3) + v(3)*eye(6))\eye(6);
        alpha = Vtot*F0 - v(2)*Fgb(:,:,2) - v(1)*Fgb(:,:,1) + (v(2)*Cinv(:,:,2)*C(:,:,3) + v(1)*Cinv(:,:,1)*C(:,:,3))*Fgb(:,:,3);

        dfstar_dhpq = xi_inv^2*(-eye(6) + Cinv(:,:,2)*C(:,:,3))*alpha + xi_inv*(-Fgb(:,:,2) + Cinv(:,:,2)*C(:,:,3)*Fgb(:,:,3));

        Fstardh_p_ana = unpackF(dfstar_dhpq);
   
   
    % ***** vm > 0 *****
    elseif(vm >= 0)
       
        vp = v(3)+ v(2); vq = v(1);
       
        % dfdh fstar 1/q -- v3 const --
        xi_inv = (v(1)*eye(6) + abs(v(2))*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))\eye(6);
        alpha = Vtot*F0 - v(2)*Fgb(:,:,2) - v(3)*Fgb(:,:,3) + (v(2)*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))*Fgb(:,:,1);

        dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,2)*C(:,:,1))*alpha + xi_inv*(-Fgb(:,:,2) + Cinv(:,:,2)*C(:,:,1)*Fgb(:,:,1));

        Fstardh_q_ana = unpackF(dfstar_dhpq);

%         % for dfdh fstar2/pq
%         xi_inv = (v(1)*Cinv(:,:,1)*C(:,:,2) + v(2)*eye(6) + v(3)*Cinv(:,:,3)*C(:,:,2))\eye(6);
%         alpha = Vtot*F0 - v(3)*Fgb(:,:,3)-v(1)*Fgb(:,:,1) + (v(3)*Cinv(:,:,3)*C(:,:,2)+v(1)*Cinv(:,:,1)*C(:,:,2))*Fgb(:,:,2);
% 
%         dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,1)*C(:,:,2))*alpha + xi_inv*(Fgb(:,:,1) - Cinv(:,:,1)*C(:,:,2)*Fgb(:,:,2));
% 
%         Fstardh_pq_ana = unpackF(dfstar_dhpq);
       
        % dfdh fstar 3/p
        xi_inv = (v(1)*Cinv(:,:,1)*C(:,:,3)+ abs(v(2))*Cinv(:,:,2)*C(:,:,3) + v(3)*eye(6))\eye(6);
        alpha = Vtot*F0 - v(2)*Fgb(:,:,2) - v(1)*Fgb(:,:,1) + (v(2)*Cinv(:,:,2)*C(:,:,3) + v(1)*Cinv(:,:,1)*C(:,:,3))*Fgb(:,:,3);

        dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,1)*C(:,:,3))*alpha + xi_inv*(-Fgb(:,:,2) + Cinv(:,:,2)*C(:,:,3)*Fgb(:,:,3));

        Fstardh_p_ana = unpackF(dfstar_dhpq);
       
    else % vm = 0
       
%         vp = v(3); vq = v(1);
%         % for dfdh fstar1/p
%         xi_inv = (v(1)*eye(6) + v(3)*Cinv(:,:,3)*C(:,:,1))\eye(6);
%         alpha = Vtot*F0 + v(3)*(-Fgb(:,:,3) + Cinv(:,:,3)*C(:,:,1)*Fgb(:,:,1));
% 
%         dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,3)*C(:,:,1))*alpha + xi_inv*(Fgb(:,:,3) - Cinv(:,:,3)*C(:,:,1)*Fgb(:,:,1));
%        
%         Fstardh_p_ana = unpackF(dfstar_dhpq);
%        
%         % for dfdh fstar 3/q
%         xi_inv = (v(3)*eye(6) + v(1)*Cinv(:,:,1)*C(:,:,3))\eye(6);
%         alpha = Vtot*F0 + v(1)*(-Fgb(:,:,1) + Cinv(:,:,1)*C(:,:,3)*Fgb(:,:,3));
% 
%         dfstar_dhpq = -xi_inv^2*(-eye(6) + Cinv(:,:,1)*C(:,:,3))*alpha + xi_inv*(-Fgb(:,:,1) + Cinv(:,:,1)*C(:,:,3)*Fgb(:,:,3));
%        
%         Fstardh_q_ana = unpackF(dfstar_dhpq);
%         Fstardh_pq_ana = Fstardh_q_ana;
       
    end    
    
%    [Fstardh_p,Fstardh_q,Fstardh_pq] = dFstr_num(Fstr,vm);
    
%    t1 = Fstardh_q_ana - Fstardh_q
%    t2 = Fstardh_pq_ana - Fstardh_pq
%    t3 = Fstardh_p_ana - Fstardh_p
    if(vm < 0)
        vp = v(3); vq = v(1) + abs(v(2)); % df_q = df_pq
    elseif(vm > 0)
        vp = v(3)+ abs(v(2)); vq = v(1); % df_p = df_pq
    else % vm = 0
        vp = v(3); vq = v(1);
    end

    Fdiff = Fstr(:,:,1)-Fstr(:,:,3);
    P = C(:,:,1)*(Fstr(:,:,1)-Fgb(:,:,1));
    Fdiff = unpackF(Fdiff);
    P = unpack(P);
   
    for i = 1:3
        for j = 1:3
            dhpq = dhpq + (1/2*Fdiff(i,j) + vp/Vtot*Fstardh_p_ana(i,j) + vq/Vtot*Fstardh_q_ana(i,j))*P(i,j);
%             dhpq = dhpq +  Fstardh_q(i,j)*P(i,j);
        end
    end
   
    Fdiff = Fstr(:,:,3)-Fstr(:,:,1);
    Fdiff = unpackF(Fdiff);
    for i = 1:3
        for j = 1:3
            dhqp = dhqp + (1/2*Fdiff(i,j) + vp/Vtot*Fstardh_p_ana(i,j) + vq/Vtot*Fstardh_q_ana(i,j))*P(i,j);
%             dhqp = dhqp + Fstardh_p(i,j)*P(i,j);
        end
    end
    
end

function [df_p, df_q, df_pq] = dFstr_num(Fstr,vm)
    global N Vtot C Cinv v Fgb F0 a12
    
    F1_cur = unpackF(Fstr(:,:,1));
    F2_cur = unpackF(Fstr(:,:,2));
    F3_cur = unpackF(Fstr(:,:,3));
    dv = 1e-5;
    V_save = v;
   
    if(vm > 0)

        v(1) = v(1) - dv; v(2) = v(2) + dv; 
        Fstar = FstarAnalytic(vm);
        F1_nxt = unpackF(Fstar(:,:,1));
        F2_nxt = unpackF(Fstar(:,:,2));
        F3_nxt = unpackF(Fstar(:,:,3));

        df_q = (F1_nxt - F1_cur)/dv;
        df_pq = (F2_nxt - F2_cur)/dv;
        df_p = (F3_nxt - F3_cur)/dv;
        
    elseif(vm < 0)
        
        v(3) = v(3) - dv; v(2) = v(2) + dv; 
        Fstar = FstarAnalytic(vm);
        F1_nxt = unpackF(Fstar(:,:,1));
        F2_nxt = unpackF(Fstar(:,:,2));
        F3_nxt = unpackF(Fstar(:,:,3));

        df_q = (F1_nxt - F1_cur)/dv;
        df_pq = (F2_nxt - F2_cur)/dv;
        df_p = (F3_nxt - F3_cur)/dv;
        
    else % vm == 0
        xi_inv = (v(1)*eye(6) + v(3)*Cinv(:,:,3)*C(:,:,1))\eye(6);
        alpha = Vtot*F0 + v(3)*(-Fgb(:,:,3) + Cinv(:,:,3)*C(:,:,1)*Fgb(:,:,1));

        dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,3)*C(:,:,1))*alpha + xi_inv*(Fgb(:,:,3) - Cinv(:,:,3)*C(:,:,1)*Fgb(:,:,1));
       
        df_p = unpackF(dfstar_dhpq);
       
        % for dfdh fstar 3/q
        xi_inv = (v(3)*eye(6) + v(1)*Cinv(:,:,1)*C(:,:,3))\eye(6);
        alpha = Vtot*F0 + v(1)*(-Fgb(:,:,1) + Cinv(:,:,1)*C(:,:,3)*Fgb(:,:,3));

        dfstar_dhpq = -xi_inv^2*(-eye(6) + Cinv(:,:,1)*C(:,:,3))*alpha + xi_inv*(-Fgb(:,:,1) + Cinv(:,:,1)*C(:,:,3)*Fgb(:,:,3));
       
        df_q = unpackF(dfstar_dhpq);
        df_pq = df_q;
        
    end

    v = V_save;

end

function [CC,CCinv,v,Fgb] = init(angle)
    global N Vtot
   
    CC = zeros(6,6,N);
    CCinv = CC; v = zeros(N,1); Fgb = zeros(3,3,N);
   
    C11 = 169.3097; C12 = 87.2201 ; C44 = 41.0448;
%     C11 = 169.3097; C12 = 122.5; C44 = 76;
   
    C =       [C11, C12, C12, 0, 0, 0;
          C12, C11, C12, 0, 0, 0;
          C12, C12, C11, 0, 0, 0;
          0, 0, 0, C44, 0, 0;
          0, 0, 0, 0, C44, 0;
          0, 0, 0, 0, 0, C44];
     
    for i = 1:N
       
        R =       [cos(angle(i))^2, sin(angle(i))^2, 0, 0, 0, -2*sin(angle(i))*cos(angle(i));
          sin(angle(i))^2, cos(angle(i))^2 ,0 ,0 ,0 ,2*sin(angle(i))*cos(angle(i));
          0, 0, 1, 0, 0, 0;
          0, 0, 0, cos(angle(i)),sin(angle(i)), 0;
          0 ,0 ,0 ,-sin(angle(i)), cos(angle(i)) ,0;
          sin(angle(i))*cos(angle(i)) ,-sin(angle(i))*cos(angle(i)) ,0 ,0 ,0 ,2*cos(angle(i))^2-1];
     
      CC(:,:,i) = R'*C*R;
      CCinv(:,:,i) = CC(:,:,i)\eye(6);
     
      v(i) = Vtot/N;
      Fgb(:,:,i) = eye(3);
    end
    Fgb(1,2,2) = .42;
    %CC(:,:,2) = zeros(6,6);
    %CCinv(:,:,2) = zeros(6,6);
    v(1) = 0.5; v(2) = 0; v(3) = Vtot - v(1);
end