clear all
global N Vtot C Cinv v xi_inv alpha Fgb F0 edges D Aa I grain_map edge_map ss tt a Neighbor_v Neighbor_e
notes = '';
% load C:\Users\icrman\Documents\Matlab\mtex-5.8.0\userScripts\GrainData2.mat
% load C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\GrainData_GBeng_v2.mat; fname = 'test_case20.mat';
% load C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\GrainData_GB_NONeng_v2.mat; fname = 'GBresults_GB_NONeng_v2.mat';
load microstructureData_Bimodal_N3.mat

% ******** Bicrystal init
% need to change last val in init function to use "th" and e1 e2 e3 in setC
% s = 1;t = 2; notes = 'Bicrystal';
% Grain_vol = [0.5 0.5]; a = 1;
% ph_i1 = [pi/12 -pi/12].*rand(2,1); Ph_i = [pi/8 -pi/8].*rand(2,1); ph_i2 = [pi/10 -pi/10].*rand(2,1);
% th = [-53.1301/2 53.1301/2]*pi/180; % sigma 5 GB
% phi1 = .005; phi0 = 0.6;
% Fpq = [1 -0.8 0;
%        0 1 0;
%        0 0 1];
% th_gb = 0*pi/180;
% rot1 = [cos(th_gb) 0 -sin(th_gb);
%         0  1  0;
%      sin(th_gb) 0 cos(th_gb)];
% Fpq = rot1'*Fpq*rot1;
% ******** Tricrystal vert
% s = [1 2]; t = [2 3]; notes = 'TricrystalVert';
% Grain_vol =[1/3; 1/3;1/3]; a = [0.5 0.5];
% ph_i1 = [0 -pi/12 pi/5]; Ph_i = [0 -pi/8 pi/4]; ph_i2 = [0 -pi/10 pi/7];
% th = [-pi/6 pi/6 pi/6];
% phi1 = [.005 .005]; phi0 = [0.2 0.2]; Fpq = zeros(3,3,2);
% Fpq(:,:,1) = [1 0.8 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,2) = [1 -0.8 0;
%        0 1 0;
%        0 0 1];
% % ******** 5 stack crystal vert
% s = [1 2 3 4]; t = [2 3 4 5]; notes = '5Stack';
% Grain_vol =[1/5; 1/5; 1/5; 1/5; 1/5]; a = [0.5 0.5 0.5 0.5];
% ph_i1 = [pi/12 -pi/12 pi/12 -pi/12 pi/12].*rand(1,5); Ph_i = [pi/8 -pi/8 pi/8 -pi/8 pi/8].*rand(1,5); ph_i2 = [pi/10 -pi/10 pi/10 -pi/10 pi/10].*rand(1,5);
% th = [-pi/6 pi/6 pi/6];
% phi1 = [.005 .005 .005 .005]; phi0 = [0.2 0.2 0.2 0.2]; Fpq = zeros(3,3,4);
% Fpq(:,:,1) = [1 0.8 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,2) = [1 -0.8 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,3) = [1 0.8 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,4) = [1 -0.8 0;
%        0 1 0;
%        0 0 1];
% ******** Tricrystal triple point
% s = [1 2 3]; t = [2 3 1]; notes = 'TriplePoint';
% ph_i1 = [pi/12 -pi/12 pi/8]; Ph_i = [pi/8 -pi/8 pi/12]; ph_i2 = [pi/10 -pi/10 pi/6];
% a = [1/2 sqrt(5)/4 sqrt(5)/4];
% Grain_vol = [3/16*2 ; 3/16*2; 1/8*2]; % rect. base = 1 height = .5
% phi1 = [.005 .005 .005]; phi0 = [0.2 0.35 0.2]; Fpq = zeros(3,3,3);
% Fpq(:,:,1) = [1 0.8 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,2) = [1 -0.8 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,3) = [1 0.4 0;
%        0 1 0;
%        0 0 1];
% th_gb = atan(.5); 
% rot1 = [1 0 0;
%     0 cos(th_gb) sin(th_gb);
%     0 -sin(th_gb) cos(th_gb)];
% Fpq(:,:,2) = rot1'*Fpq(:,:,2)*rot1;
% Fpq(:,:,3) = rot1*Fpq(:,:,3)*rot1';

N = length(s)+length(Grain_vol); edges = length(s); 
Grain_vol = Grain_vol./sum(Grain_vol); a = a./sum(a);
Vtot = sum(Grain_vol);

grain_map = 1:N-edges;
edge_map = N-edges+1:N;

F0 = [1,0 ,0;
      0, 1 ,0;
      0, 0 ,1];
Fgb = zeros(9,1,N);
xi_inv = zeros(9,9,N-edges); 
alpha= zeros(9,1,N-edges); 
[C,Cinv,v,Fgbt] = init(ph_i1,Ph_i,ph_i2,Grain_vol,Fpq,grain_map, edge_map,zeros(N,1));
F0 = reshape(F0, [9,1]);
for i = 1:N
    Fgb(:,:,i) = reshape(Fgbt(:,:,i), [9,1]);
end
%% 20-grains

dt = 2e-6;
tf = floor(1/dt);
endt = 5000;
ss = s;
tt = t;
D = digraph(s,t,a);
Aa = sparse(adjacency(D,a)); 
I = sparse(incidence(D)*diag(a));
% I = incidence(D);

stp_e = zeros(edges,1);
stp_v = zeros(N-edges,1);
stress_shearMap = zeros(endt,N-edges);
A = zeros(endt,N-edges);
dhpq = 0; dhqp = 0; doth = zeros(edges,1);
h = zeros(endt,edges);
V = zeros(endt,N-edges);
for i = 1:N-edges
    p = grain_map(i);
    V(1,i) = v(p);
end
% get neighbors
Neighbor_v =  cell(N-edges,1); % neighboring grains to grain p
Neighbor_e = cell(N-edges,1); % neighboring gbs to grain p
Aa_temp = Aa + Aa';
for vertex = 1:N-edges 
    Neighbor_v{vertex} = find(Aa_temp(vertex,:));
    Neighbor_e{vertex} = find(I(vertex,:));
end
% GetGammaS(Grain_vol)

vol = Grain_vol;
V0 = vol; V(1,:) = V0;

stress = zeros(endt,1);
strain = zeros(endt,1);
time = zeros(endt,1);

n = .3;
m = 10;
phi_h1 = .003*0;
phi_h2 = 2e-5;
for i = 1:edges
    p = grain_map(ss(i)); q = grain_map(tt(i));
    vv = [v(p);v(q)];
    if(rand(1,1) < 0)
        kappa(i) = 0;
    else
        kappa(i) = 1/(a(i)*min(vv))*0.005;
    end
end
aeff = a;
for i = 1:edges
    if(phi0(i) < 0.1)
        phi0(i) = phi0(i)*10;
    end
    while(phi1(i) > 0.005)
        phi1(i) = phi1(i)*0.6;
    end
end
% phi1 = phi0*.1;
%  phi1 = rand(edges,1)*.01 + .001;
% phi0 = rand(edges,1)*1.5 + .5;
for t = 1:endt

%     F0(1) = 1-0.25*t/endt;
%     F0(5) = 1/sqrt(1-0.25*t/endt);
%     F0(9) = 1/sqrt(1-0.25*t/endt);
    F0(4) = 0.25*t/endt;
%     F0(4) =  F0(4) + getSrain(t,100,endt);
    
    Fstar = FstarAnalytic();
    P0 = GetP0(Fstar);
    for edge = 1:edges
        
        if(stp_e(edge)) % stop gb motion to conserve volume
            doth(edge) = 0;
            continue
        end
        q = grain_map(ss(edge)); p = grain_map(tt(edge));  
%         [dFpq, dFqp] = Gethdot(Fstar,P0,edge,vol);
        [dFpq, dFqp] = Gethdot_Neighbors(Fstar,P0,edge,vol);

        phi_h = phi_h1*abs(h(t,edge)*kappa(edge)*aeff(edge))^n + phi_h2*abs(h(t,edge)*kappa(edge)*aeff(edge))^m;
        dhpq = -1/phi1(edge)*(dFpq + phi0(edge) + phi_h);
        dhqp = -1/phi1(edge)*(dFqp - phi0(edge) - phi_h);
        

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

        doth(edge) = dhpq; bea = -2*abs(h(t,edge))/.3;
        if(bea == 0) 
            continue; 
        else 
            aeff(edge) = .5*(a(edge)*sqrt(bea^2*a(edge)^2 + 1) + asinh(bea*a(edge))/bea ); 
        end

    end
    
    [dv,dh,stp_e,stp_v] = UpdateVol(doth,dt);
    vol = vol + dv;
    for i = 1:N-edges
        V(t+1,i) = vol(i);
    end
    for i = 1:edges
        h(t+1,i) = h(t,i) + dh(i);
    end

    time(t+1) = t*dt;

    stress(t) = P0(4);
    stress_shearMap(t,:) = Fstar(4,:);
    strain(t) = F0(4); 
    if(mod(t,100) == 0)
        fprintf('t = %f\n',t)
    end
    
    [stp_e,stp_v] = stop(vol);
%     if(P0(4) < 0)
%         fprintf('STOP! t = %f\n',t)
%         break
%     end
    for i = 1:N-edges
        for j = 1:9
            A(t,i) = A(t,i) + Fstar(j,:,i)*P0(j);
        end
    end
end
%%
figure(2)

hold on
plot(time,h,'LineWidth',1.3)
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
hold on
plot(strain,stress,'LineWidth',1.3)
% xline([strain(floor(1500+endt/4)) strain(floor(1500+endt*1.3/4)) strain(floor(1500+endt*3/4))])
xlabel('Strain')
ylabel('Stress')
grid on
hold off 

%%
np = pwd;
cd 'C:\NP_results'
fname = append('Results',num2str(randi([1 90000],1,1)));
fname = append(fname,append('_g',num2str(length(Grain_vol))));
fname = append(fname,notes);
np_data_dir = append('C:\NP_results\',fname);
mkdir(fname)
cd(fname)
save GBresults.mat stress strain V time h stress_shearMap A
save input.mat a Fpq ph_i1 Ph_i ph_i2 phi0 phi1 s t
polymap = 'poly'; polymap = append(polymap,append('_',num2str(length(Grain_vol))));
polymap = append(polymap,'.mat');
cd(np)
if(isfile(polymap))
    copyfile(polymap,np_data_dir)
else
    polymap = 'Poly_N*';
    polymap = append(polymap,append('_g',num2str(length(Grain_vol))));
    polymap = append(polymap,'.mat');
    list = dir(polymap);
    if(~isempty(list))
        copyfile(list.name,np_data_dir)
    else
        disp('no poly.mat file found')
    end
end
% save GBresults_bicrystal_90.mat stress strain V time h stress_shearMap A

% save GBresults_500_full.mat stress strain V time h stress_shearMap

function [stp_e, stp_v]  = stop(vol)
global ss tt N edges Aa I
    stp_e = zeros(edges,1);
    stp_v = zeros(N-edges,1);
    for i = 1:length(vol)
        if(vol(i) <= 0)
            neigh = abs(I(i,:));
            stp_v(i) = true;
            for j = 1:edges
                if(neigh(j) > 0)
                    stp_e(j) = true;
                end 
            end
        end
    end
end
function [vol, h,stp_e,stp_v] = UpdateVol(dh,dt)
global N edges Vtot C Cinv v Fgb F0  a D Aa I grain_map edge_map

    stp_e = zeros(edges,1);
    stp_v = zeros(N-edges,1);

    h = dh.*dt; ha = abs(h);
    Ah = sparse(adjacency(D,h));
    Aha = sparse(adjacency(D,ha));
    
    for i = 1:edges % calculate Vpq
       pq = edge_map(i);
       v(pq) = v(pq) + a(i)*h(i);
    end
    
    gimal = (Aa' + Aa);
    gimal_h = (Ah - Ah');
    gimal_ha = (Aha + Aha');
    
    dV = -1/2*(gimal_ha + gimal_h')*gimal;
    dV = diag(dV);
    
    for i = 1:N-edges % add dv to vertices
       p = grain_map(i);
       if(v(p) <= 0)
           dV(i) = 0;
           neigh = abs(I(i,:));
           stp_v(i) = true;
           for j = 1:edges
                if(neigh(j) > 0)
                    stp_e(j) = true;
                end 
           end
       end
       v(p) = v(p) + dV(i);
    end
    
    vol = -I*h; % volume of grains + quasi-grains
end

function [dhpq, dhqp] = Gethdot_Neighbors(Fstr,P0,i,vol)
    global N Vtot C Cinv  Fgb F0 a12 xi_inv edges alpha grain_map edge_map ss tt Neighbor_v Neighbor_e I 
    
    V_e = Neighbor_e{ss(i)};
    V_v = Neighbor_v{ss(i)};
    Fstardh_ana_pq = zeros(3,3); Fstardh_ana_qp = zeros(3,3); Fstardh_num_pq = zeros(3,3);
    dhpq = 0; dhqp = 0; Ii = [1;0;0;0;1;0;0;0;1];
    % ****** vm < 0 ******
    pq = edge_map(i);
    counter = 1;

    pi = tt(i); qi = ss(i);
    C_diff = (Cinv(:,:,pi)-Cinv(:,:,qi)); 

    XI = (Cinv(:,:,pi)-Cinv(:,:,qi))*(P0 + C(:,:,1)*Fgb(:,:,1)); 
    Fstardh_ana_pq = (Ii - Fgb(:,:,pq));
    Fstardh_ana_pq = reshape(Fstardh_ana_pq,[3,3]);
    XI = reshape(XI,[3,3]);

    Fdiff = Fstr(:,:,qi)-Fstr(:,:,pi);

    Fdiff = reshape(Fdiff,[3,3]);
    P0 = reshape(P0,[3,3]);
   
    for k = 1:3
        for j = 1:3
            dhpq = dhpq + (1/2*Fdiff(k,j) + (Fstardh_ana_pq(k,j)+XI(k,j)))*P0(k,j);
        end
    end
   
    for k = 1:3
        for j = 1:3
          dhqp = dhqp + (-1/2*Fdiff(k,j) + (Fstardh_ana_pq(k,j)+XI(k,j)))*P0(k,j);
        end
    end
   
end

function [dhpq, dhqp] = Gethdot(Fstr,P0,i,vol)
    global N Vtot C Cinv  Fgb F0 a12 xi_inv edges alpha grain_map edge_map ss tt
   
    Fstardh_ana_pq = zeros(3,3,N-edges); Fstardh_ana_qp = zeros(3,3,N-edges);
    dhpq = 0; dhqp = 0; I = [1;0;0;0;1;0;0;0;1];
    % ****** vm < 0 ******
    pq = edge_map(i);
    
    counter = 1;
    for pp = 1:N-edges
        p = grain_map(ss(pp)); q = grain_map(tt(pp));  

        C_diff = Cinv(:,:,q)*C(:,:,p);
        dfstar_dhpq = xi_inv(:,:,tt(pp))*(C_diff*Fgb(:,:,p)-Fgb(:,:,pq)-I) - xi_inv(:,:,tt(pp))^2*(C_diff-eye(9))*alpha(:,:,tt(pp));
        Fstardh_ana_pq(:,:,counter) = reshape(dfstar_dhpq,[3,3]);
%         Fstardh_ana_qp(:,:,counter) = reshape(dfstar_dhqp,[3,3]);
        counter = counter + 1;
            
    end

    Fstardh_pq = zeros(3,3); 
%     Fstardh_qp = zeros(3,3);
    for j = 1:N-edges
        Fstardh_pq = Fstardh_pq + Fstardh_ana_pq(:,:,j)*vol(j);
%         Fstardh_qp = Fstardh_qp + Fstardh_ana_qp(:,:,j)*vol(j);
    end

    q = ss(i); p = tt(i);  
    Fdiff = Fstr(:,:,q)-Fstr(:,:,p);

    Fdiff = reshape(Fdiff,[3,3]);
    P0 = reshape(P0,[3,3]);
   
    for i = 1:3
        for j = 1:3
            dhpq = dhpq + (1/2*Fdiff(i,j) + Fstardh_pq(i,j))*P0(i,j);
        end
    end
   
    for i = 1:3
        for j = 1:3
          dhqp = dhqp + (-1/2*Fdiff(i,j) + Fstardh_pq(i,j))*P0(i,j);
        end
    end
   
end

% function [Fstr] = FstarAnalytic()
%     global N edges Vtot C Cinv v Fgb F0 xi_inv alpha grain_map
%    
%     Fstr = zeros(9,1,N-edges);
%     for i = 1:N-edges
%         index = grain_map(i,3);
%        
%         xi = v(index)*eye(9);
%         alpha_temp = Vtot*F0;
%         for j = 1:N            
%             if(i == j);continue; end        
%             xi = xi + abs(v(j))*Cinv(:,:,j)*C(:,:,index);
%             alpha_temp = alpha_temp - v(j)*Fgb(:,:,j);
%         end
%         alpha(:,:,i) = alpha_temp + (xi-v(index)*eye(9))*Fgb(:,:,index);
%         if(v(i) == 0)
%             Fstr(:,:,i) = Fgb(:,:,index);
%             xi_inv(:,:,i) = eye(9);
%         else
%             xi_inv(:,:,i) = xi\eye(9);
%             Fstr(:,:,i) = xi_inv(:,:,i)*alpha(:,:,i);
%         end
%     end
%    
% end
function [Fstr] = FstarAnalytic()
    global N edges Vtot C Cinv v Fgb F0 xi_inv alpha grain_map
   
    Fstr = zeros(9,1,N-edges);
    xi = zeros(9,9);
    alpha_temp = Vtot*F0;
    for i = 1:N
        xi = xi + abs(v(i))*Cinv(:,:,i);
        alpha_temp = alpha_temp - v(i)*Fgb(:,:,i);
    end
    for i = 1:N-edges
        index = grain_map(i);
        xi_p = (xi - v(index)*Cinv(:,:,index))*C(:,:,index) + v(index)*eye(9);

        alpha(:,:,i) = alpha_temp + v(index)*Fgb(:,:,index) + (xi_p - v(index)*eye(9))*Fgb(:,:,index);
        xi_inv(:,:,i) = xi_p\eye(9);
        Fstr(:,:,i) = xi_inv(:,:,i)*alpha(:,:,i);
    end
   
end
function [C,Cinv] = SetC(p1,p,p2,th)

    % Copper elastic constants
%     C11 = 169.3097; C12 = 122.5; C44 = 76; % cubic GPa
%     C11 = 169.3097; C12 = 87.2201 ; C44 = 41.0448; % Isotropic GPa

    % Aluminum elastic constants
    C11 = 116.3; C12 = 64.8; C44 = 30.9; % cubic GPa

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

%     e1 = [cos(th);sin(th);0];
%     e2 = [-sin(th);cos(th);0];
%     e3 = [0;0;1];
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
        index = grain_map(i);
        if(v(index) <= 0)
            continue
        else
            P = P + C(:,:,index)*(Fstr(:,:,i)-Fgb(:,:,index));
            break
        end
    end
end
function [CC,CCinv,v,Fgb] = init(p1,p,p2,vol,Fpq,v_map,e_map,th)
    global N Vtot edges
   
    CC = zeros(9,9,N);
    CCinv = CC; v = zeros(N,1); Fgb = zeros(3,3,N);
    for i = 1:N-edges
      index = v_map(i);
      [CC(:,:,index),CCinv(:,:,index)] = SetC(p1(i),p(i),p2(i),th(i));
        
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
function e = getSrain(t,ints,intf)

%     interval = floor(abs(intf - ints)/4);
%     interval = floor(abs(intf - ints)/5);
    f = .03;
    e = 1e-3*f;
%     if(t < ints + interval)
%         e = 1e-3*f;
%     else
%         e = 1e-4*f;
%     end

%     if(t < ints + interval)
%         e = 1e-4*f;
%     elseif(t >= ints + interval && t < ints + interval*1.3 )
%         e = 1e-3*f;
%     elseif(t >= ints + interval*1.3 && t < ints + interval*3)
%         e = 1e-5*f;
%     elseif(t >= ints + interval*3 && t < ints + interval*4)
%         e = 1e-4*f;
%     else
%         e = 1e-4*f;
%     end
end
function GetGammaS(Grain_vol)
global N edges Neighbor_v Neighbor_e

    omega = 6/((N-edges));
  
    for i = 1:N-edges
        
        N_p = Neighbor_v{i};
        N_pnext = N_p;
        N_pq = Neighbor_e{i};
        N_pqnext = N_pq;
        vol = Grain_vol(N_p);
        vtf = vol < omega/10;
        while sum(vtf) > 0
            for j = 1:length(vtf)
                if(vtf(j))
                    VV = Neighbor_v{N_pnext(j)};
                    PQ = Neighbor_e{N_pnext(j)};
                    N_pnext = [N_pnext VV];
                    N_pqnext = [N_pqnext PQ];
                end
            end

            N_pnext = unique(N_pnext);
            N_pnext = setdiff(N_pnext,N_p);

            N_pqnext = unique(N_pqnext);
            N_pqnext = setdiff(N_pqnext,N_pq);

            vol = Grain_vol(N_pnext);
            vtf = vol < omega/10;
            N_p = [N_p N_pnext];
            N_pq = [N_pq N_pqnext];
            if(sum(Grain_vol(N_p)) >= omega)
                break
            end
        end
        N_p = unique(N_p);
        N_pq = unique(N_pq);
        for j = 1:length(N_p)
            if(N_p(j) == i)
                N_p(j) = [];
                break
            end
        end
        Neighbor_v{i} = N_p;
        Neighbor_e{i} = N_pq;
        vtest(i) = sum(Grain_vol(N_p));
    end

end
function [dFdh] = Numbericaldfdh(pq,h)
 global v N edges 

 dFdh = zeros(9,1,N-edges);
 v_save = v;
 h_save = h;
 dh = 1e-6;

 h(pq) = h(pq) + dh;
 UpdateVol(h,1);
 fstr_pdh = FstarAnalytic();

 v = v_save;
 h = h_save;
 h(pq) = h(pq) - dh;
 UpdateVol(h,1);
 fstr_ndh = FstarAnalytic();

 for i = 1:N-edges
    dFdh(:,:,i) = (fstr_pdh(:,:,i) - fstr_ndh(:,:,i))/(2*dh);
 end
 v = v_save;
end
