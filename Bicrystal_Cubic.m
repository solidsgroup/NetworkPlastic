clear all
global N Vtot C Cinv Fgb v F0

N = 3;
Vtot = 1;
th = [-pi/12, -pi/12, pi/12];

[C, Cinv,v, Fgbt] = init(th);

F = [1,0.6,0;
    0,1,0;
    0,0,1];
F0 = reshape(F, [9,1]);
for i = 1:N
    Fgb(:,:,i) = reshape(Fgbt(:,:,i), [9,1]);
end

vm = .4;
v(1) = .5; v(2) = vm; v(3) = Vtot-(vm+v(3));
Fstr = FstarAnalytic(vm);

v(1)*Fstr(:,:,1) + v(2)*Fstr(:,:,2) + v(3)*Fstr(:,:,3) - Vtot*F0
C(:,:,1)*(Fstr(:,:,1)-Fgb(:,:,1)) + C(:,:,2)*(Fstr(:,:,2)-Fgb(:,:,2)) - 2*C(:,:,3)*(Fstr(:,:,3)-Fgb(:,:,3))

Fstarq = reshape(Fstr(:,:,1),[3,3]);
Fstarpq = reshape(Fstr(:,:,2),[3,3]);
Fstarp = reshape(Fstr(:,:,3),[3,3]);

% P1 = C(:,:,1)*(Fstr(:,:,1)-Fgb(:,:,1));
% P1 = reshape(P1,[3,3])
% P2 = mul(C(:,:,1),Fstarq-eye(3))

% C_test = [C(1,1) C(1,5) C(1,9) C(1,2) C(1,6) C(1,7);
%          C(5,1) C(5,5), C(5,9) C(5,2) C(5,6) C(5,7);
%          C(9,1) C(9,5) C(9,9) C(9,2) C(9,6) C(9,7)];
% P_test = C_test*[F0(1,1);F0(2,2);F0(3,3);0;0;0]    
%
% P = mul(C1,F0)
%
% F = mul(Cinv1,P)



%% Bicrystal

v(1) = 0.5; v(2) = 0; v(3) = Vtot - v(1);

dt = 5e-6;
a12 = .5;
tf = floor(1/dt);
dh12 = 0;
dh21 = 0;
phi1 = .006;
phi0 = 0.2;

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

    F0(4) = .2/tf*t;
%     F0(1) = 1-.1/tf*t;

    Fstar = FstarAnalytic(vm);
    
%     Fstarq = reshape(Fstar(:,:,1),[3,3]);
%     Fstarpq = reshape(Fstar(:,:,2),[3,3]);
%     Fstarp = reshape(Fstar(:,:,3),[3,3]);
    
    [dhpq, dhqp] = Gethdot(Fstar,vm);
    
    dh12 = -1/phi1*(dhpq + phi0);
    dh21 = -1/phi1*(dhqp - phi0);

    AA1(t+1) = dhqp;

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
        C(:,:,2) = zeros(9,9);
        Cinv(:,:,2) = zeros(9,9);
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
    stress(t+1) = temp(4);
    strain(t+1) = F0(4);
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

function [dhpq, dhqp] = Gethdot(Fstr,vm)
    global N Vtot C Cinv v Fgb F0 a12
   
    dhpq = 0; dhqp = 0;
   
    % ****** vm < 0 ******
    if(vm < 0)
       
        vp = v(3); vq = v(1) + v(2);
       
        % for dfdh fstar2/q -- v1 const --
       
        xi_inv = (v(1)*eye(9) + abs(v(2))*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))\eye(9);
        alpha = Vtot*F0 - v(2)*Fgb(:,:,2)-v(3)*Fgb(:,:,3) + (v(2)*Cinv(:,:,2)*C(:,:,1)+v(3)*Cinv(:,:,3)*C(:,:,1))*Fgb(:,:,1);

        dfstar_dhpq = -xi_inv^2*(eye(9) - Cinv(:,:,3)*C(:,:,1))*alpha + xi_inv*(-Fgb(:,:,2) + Cinv(:,:,3)*C(:,:,1)*Fgb(:,:,1));

        Fstardh_q_ana = reshape(dfstar_dhpq,[3,3]);
       
        % for dfdh fstar2/pq
%         xi_inv = (v(1)*Cinv(:,:,1)*C(:,:,2) + v(2)*eye(6) + v(3)*Cinv(:,:,3)*C(:,:,2))\eye(6);
%         alpha = Vtot*F0 - v(3)*Fgb(:,:,3)-v(1)*Fgb(:,:,1) + (v(3)*Cinv(:,:,3)*C(:,:,2)+v(1)*Cinv(:,:,1)*C(:,:,2))*Fgb(:,:,2);
% 
%         dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,3)*C(:,:,2))*alpha + xi_inv*(Fgb(:,:,3) - Cinv(:,:,3)*C(:,:,2)*Fgb(:,:,2));
% 
%         Fstardh_pq_ana = unpackF(dfstar_dhpq);

        % dfdh fstar 3/p
        xi_inv = (v(1)*Cinv(:,:,1)*C(:,:,3)+ abs(v(2))*Cinv(:,:,2)*C(:,:,3) + v(3)*eye(9))\eye(9);
        alpha = Vtot*F0 - v(2)*Fgb(:,:,2) - v(1)*Fgb(:,:,1) + (v(2)*Cinv(:,:,2)*C(:,:,3) + v(1)*Cinv(:,:,1)*C(:,:,3))*Fgb(:,:,3);

        dfstar_dhpq = xi_inv^2*(-eye(9) + Cinv(:,:,2)*C(:,:,3))*alpha + xi_inv*(-Fgb(:,:,2) + Cinv(:,:,2)*C(:,:,3)*Fgb(:,:,3));

        Fstardh_p_ana = reshape(dfstar_dhpq,[3,3]);
   
   
    % ***** vm > 0 *****
    elseif(vm >= 0)
       
        vp = v(3)+ v(2); vq = v(1);
       
        % dfdh fstar 1/q -- v3 const --
        xi_inv = (v(1)*eye(9) + abs(v(2))*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))\eye(9);
        alpha = Vtot*F0 - v(2)*Fgb(:,:,2) - v(3)*Fgb(:,:,3) + (v(2)*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))*Fgb(:,:,1);

        dfstar_dhpq = -xi_inv^2*(eye(9) - Cinv(:,:,2)*C(:,:,1))*alpha + xi_inv*(-Fgb(:,:,2) + Cinv(:,:,2)*C(:,:,1)*Fgb(:,:,1));

        Fstardh_q_ana = reshape(dfstar_dhpq,[3,3]);

%         % for dfdh fstar2/pq
%         xi_inv = (v(1)*Cinv(:,:,1)*C(:,:,2) + v(2)*eye(6) + v(3)*Cinv(:,:,3)*C(:,:,2))\eye(6);
%         alpha = Vtot*F0 - v(3)*Fgb(:,:,3)-v(1)*Fgb(:,:,1) + (v(3)*Cinv(:,:,3)*C(:,:,2)+v(1)*Cinv(:,:,1)*C(:,:,2))*Fgb(:,:,2);
% 
%         dfstar_dhpq = -xi_inv^2*(eye(6) - Cinv(:,:,1)*C(:,:,2))*alpha + xi_inv*(Fgb(:,:,1) - Cinv(:,:,1)*C(:,:,2)*Fgb(:,:,2));
% 
%         Fstardh_pq_ana = unpackF(dfstar_dhpq);
       
        % dfdh fstar 3/p
        xi_inv = (v(1)*Cinv(:,:,1)*C(:,:,3)+ abs(v(2))*Cinv(:,:,2)*C(:,:,3) + v(3)*eye(9))\eye(9);
        alpha = Vtot*F0 - v(2)*Fgb(:,:,2) - v(1)*Fgb(:,:,1) + (v(2)*Cinv(:,:,2)*C(:,:,3) + v(1)*Cinv(:,:,1)*C(:,:,3))*Fgb(:,:,3);

        dfstar_dhpq = -xi_inv^2*(eye(9) - Cinv(:,:,1)*C(:,:,3))*alpha + xi_inv*(-Fgb(:,:,2) + Cinv(:,:,2)*C(:,:,3)*Fgb(:,:,3));

        Fstardh_p_ana = reshape(dfstar_dhpq,[3,3]);
       
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
    Fdiff = reshape(Fdiff,[3,3]);
    P = reshape(P,[3,3]);
   
    for i = 1:3
        for j = 1:3
            dhpq = dhpq + (1/2*Fdiff(i,j) + vp/Vtot*Fstardh_p_ana(i,j) + vq/Vtot*Fstardh_q_ana(i,j))*P(i,j);
%             dhpq = dhpq +  Fstardh_q(i,j)*P(i,j);
        end
    end
   
    for i = 1:3
        for j = 1:3
            dhqp = dhqp + (-1/2*Fdiff(i,j) + vp/Vtot*Fstardh_p_ana(i,j) + vq/Vtot*Fstardh_q_ana(i,j))*P(i,j);
%             dhqp = dhqp + Fstardh_p(i,j)*P(i,j);
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
function P = mul(C,F)
   
    P = zeros(3,3);
    counter = 1;
    I = 1; J = 1;
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                   
                    P(i,j) = P(i,j) + C(I,J)*F(k,l);
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

function [CC,CCinv,v,Fgb] = init(angle)
    global N Vtot
   
    CC = zeros(9,9,N);
    CCinv = CC; v = zeros(N,1); Fgb = zeros(3,3,N);
     
    for i = 1:N
     
     [CC(:,:,i),CCinv(:,:,i)] = SetC(angle(i));
     
      v(i) = Vtot/N;
      Fgb(:,:,i) = eye(3);
    end
    Fgb(1,2,2) = .2;
%     Fgb(3,1,2) = -.2;
    %CC(:,:,2) = zeros(6,6);
    %CCinv(:,:,2) = zeros(6,6);
    v(1) = 0.5; v(2) = 0; v(3) = Vtot - v(1);
end

function Fstr = FstarAnalytic(vm)
    global N Vtot C Cinv v Fgb F0
   
    Fstr = zeros(9,1,N);
   
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
        temp = (v(1)*eye(9) + abs(v(2))*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))\(Vtot*F0 - v(2)*Fgb(:,:,2) - v(3)*Fgb(:,:,3) + ...
            (abs(v(2))*Cinv(:,:,2)*C(:,:,1) + v(3)*Cinv(:,:,3)*C(:,:,1))*Fgb(:,:,1) );
        Fstr(:,:,1) = temp;
       
        if(vm ~= 0)
            temp = (abs(v(2))*eye(9) + v(1)*Cinv(:,:,1)*C(:,:,2) + v(3)*Cinv(:,:,3)*C(:,:,2))\(Vtot*F0 - v(1)*Fgb(:,:,1) - v(3)*Fgb(:,:,3) +...
                (v(1)*Cinv(:,:,1)*C(:,:,2) + v(3)*Cinv(:,:,3)*C(:,:,2))*Fgb(:,:,2) );
            Fstr(:,:,2) = temp;
        else
            Fstr(:,:,2) = temp;
           
        end
       
        temp = (v(3)*eye(9) + v(1)*Cinv(:,:,1)*C(:,:,3) + abs(v(2))*Cinv(:,:,2)*C(:,:,3))\(Vtot*F0 - v(1)*Fgb(:,:,1) - v(2)*Fgb(:,:,2) +...
            (v(1)*Cinv(:,:,1)*C(:,:,3) + abs(v(2))*Cinv(:,:,2)*C(:,:,3))*Fgb(:,:,3) );
        Fstr(:,:,3) = temp;
    end
end