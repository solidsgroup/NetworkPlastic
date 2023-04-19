clear all
global N Vtot C Cinv v xi_inv alpha Fgb F0 edges D Aa I grain_map edge_map ss tt a Neighbor_V2e
notes = '';
% load C:\Users\icrman\Documents\Matlab\mtex-5.8.0\userScripts\GrainData2.mat
% load C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\GrainData_GBeng_v2.mat; notes = 'GBeng_comp';
% load C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\GrainData_GB_NONeng_v2.mat; notes = 'NONGBeng_comp';
% load C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\GrainData_GB_NONeng_twin.mat; notes = 'EBSDTwin';
% load C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\NP_stat_varOri\microstructureData_GBengSubset_1.mat; notes = 'FullOutputengSubSet';
% load C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\NP_stat_varOri\microstructureData_testnon_1.mat; notes = 'FullOutputNon';
load 'GrainData_GB_NONeng_twin.mat'; notes = 'twin_Negshear';
% load 'GrainData_GB_NONeng.mat'; notes = 'nonGBeng_comp';
% ******** Bicrystal init

% % % need to change last val in init function to use "th" and e1 e2 e3 in setC
% s = 1;t = 2; notes = 'Bicrystal_tens_stress';
% Grain_vol = [0.5 0.5]; a = 1;
% ph_i1 = [pi/12 -pi/12]; Ph_i = [pi/8 -pi/8]; ph_i2 = [pi/10 -pi/10];
% % th = [-53.1301/2 53.1301/2]*pi/180; % sigma 5 GB
% % th = [-pi/12 pi/12];
% phi1 = .005; phi0 = 0.6;  
% Fpq = [1 -1 0;
%        0 1 0;
%        0 0 1];
% th_gb = 63*pi/180;
% rot1 = [cos(th_gb) -sin(th_gb) 0;
%      sin(th_gb) cos(th_gb) 0;
%      0 0 1];
% Fpq = rot1'*Fpq*rot1;

% ******** Tricrystal vert
% s = [1 2]; t = [2 3]; notes = 'TricrystalVert';
% Grain_vol =[1/3; 1/3;1/3]; a = [0.5 0.5];
% ph_i1 = [0 -pi/12 pi/5]; Ph_i = [0 -pi/8 pi/4]; ph_i2 = [0 -pi/10 pi/7];
% th = [-pi/6 pi/6 pi/6];
% phi1 = [.005 .005]; phi0 = [0.2 0.2]; Fpq = zeros(3,3,2);
% Fpq(:,:,1) = [1 -0.8 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,2) = [1 0.8 0;
%        0 1 0;
%        0 0 1];
% % ******** 5 stack crystal vert
% s = [1 2 3 4]; t = [2 3 4 5]; notes = '5Stack';
% Grain_vol =[1/5; 1/5; 1/5; 1/5; 1/5]; a = [0.5 0.5 0.5 0.5];
% ph_i1 = [pi/12 -pi/12 pi/12 -pi/12 pi/12].*rand(1,5); Ph_i = [pi/8 -pi/8 pi/8 -pi/8 pi/8].*rand(1,5); ph_i2 = [pi/10 -pi/10 pi/10 -pi/10 pi/10].*rand(1,5);
% th = [-pi/6 pi/6 pi/6];
% phi1 = [.005 .005 .005 .005]; phi0 = [0.2 0.2 0.2 0.2]; Fpq = zeros(3,3,4);
% Fpq(:,:,1) = [1 -0.8 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,2) = [1 0.8 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,3) = [1 -0.8 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,4) = [1 0.8 0;
%        0 1 0;
%        0 0 1];
% ******** Tricrystal triple point
% s = [1 2 3]; t = [2 3 1]; notes = 'TriplePoint';
% ph_i1 = [pi/12 -pi/12 pi/8]; Ph_i = [pi/8 -pi/8 pi/12]; ph_i2 = [pi/10 -pi/10 pi/6];
% a = [1/2 sqrt(5)/4 sqrt(5)/4];
% Grain_vol = [3/16*2 3/16*2 1/8*2]; % rect. base = 1 height = .5
% phi1 = [.0005 .0005 .0005]; phi0 = [0.2 0.2 0.2]; Fpq = zeros(3,3,3);
% Fpq(:,:,1) = [1 0.5 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,2) = [1 -0.5 0;
%        0 1 0;
%        0 0 1];
% Fpq(:,:,3) = [1 0.4 0;
%        0 1 0;
%        0 0 1];
% th_gb = atan(.5); 
% rot1 = [cos(th_gb) -sin(th_gb) 0;
%      sin(th_gb) cos(th_gb) 0;
%      0 0 1];
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
[C,Cinv,v,Fgbt] = init(ph_i1,Ph_i,ph_i2,Grain_vol,Fpq,grain_map, edge_map);
F0 = reshape(F0, [9,1]);
for i = 1:N
    Fgb(:,:,i) = reshape(Fgbt(:,:,i), [9,1]);
end

dt = 5e-7;
tf = floor(1/dt);
endt = 2500;
ss = s;
tt = t;
D = digraph(s,t,a);
Aa = sparse(adjacency(D,a)); 
I = sparse(incidence(D)*diag(a));

stp_e = zeros(edges,1); flag_edge = zeros(edges,1);
stp_v = zeros(N-edges,1);
dhpq = 0; dhqp = 0; doth = zeros(edges,1);
h = zeros(endt,edges);
V = zeros(endt,N-edges);
for i = 1:N-edges
    p = grain_map(i);
    V(1,i) = v(p);
end

% get neighbors
Neighbor_V2e = cell(N-edges,1); % neighboring gbs to grain p
for vertex = 1:N-edges
    Neighbor_V2e{vertex} = find(I(vertex,:));
end

vol = Grain_vol;
V0 = vol; V(1,:) = V0;

stress = zeros(endt,9,1);
strain = zeros(endt,9,1);
time = zeros(endt,1);
aa1 = zeros(endt,1);

n = .3;
m = 10;
phi_h1 = .003*0;
phi_h2 = 2e-5;
for i = 1:edges
    p = grain_map(ss(i)); q = grain_map(tt(i));
    vv = [v(p);v(q)];
%     gg = 5e-3;  
    gg = .00005;
    kappa(i) = gg/(a(i)*min(vv));
end
aeff = a;
phi1 = phi0*0.0005;

Ii = [1;0;0;0;1;0;0;0;1];
C_ref = zeros(9,9); Cinv_ref = zeros(9,9); mask = [0;0;0;0;0;0;0;0;0]; 
P0 = zeros(9,1);

load_type = 'shear negative'; strain_rate = 0.05;

if(strcmpi(load_type,'compression'))
    mask(1) = true;
    F0 = [-abs(strain_rate)/endt;0;0;0;0;0;0;0;0];
elseif(strcmpi(load_type,'tension'))
    mask(1) = true;
    F0 = [abs(strain_rate)/endt;0;0;0;0;0;0;0;0];
elseif(strcmpi(load_type,'shear positive'))
    mask(4) = true;
    F0 = [0;0;0;abs(strain_rate)/endt;0;0;0;0;0];
elseif(strcmpi(load_type,'shear negative'))
    mask(4) = true;
    F0 = [0;0;0;-abs(strain_rate)/endt;0;0;0;0;0];
else
    error('bad loading condition')
end

for t = 1:endt
    
    [C_ref, Cinv_ref, Fgb_avg] = GetAVG(vol);
    % _____Stress driven motion_____
%     P0 = [10*t/endt;0;0;10*t/endt*0;0;0;0;0;0]; 
%     for i = 1:9
%         if(mask(i))
%             if(sign(Fgb_avg(i)) + sign(F0(i)) ~= 0)
%                 Fgb_avg(i) = -Fgb_avg(i);
%             end
%         end
%     end
%     F0 = Cinv_ref*P0 + Ii - Fgb_avg;
    % _____Strain driven motion_____
    T = (Cinv_ref)\eye(9);
    F0_t = F0*t;
    P = -T*Fgb_avg;
    F0_t = CG(P, F0_t, T, mask);
    F0_t = F0_t + Ii;
    
%     Fstar = FstarAnalytic(Fgb_avg,Cinv_ref);
%     P0 = GetP0(Fstar);
    P0 = T*(Vtot*F0_t + Fgb_avg - Ii);
    for edge = 1:edges
        
        if(stp_e(edge) || edge == -409) % stop gb motion to conserve volume
            doth(edge) = 0;
            continue
        end
        q = grain_map(ss(edge)); p = grain_map(tt(edge));  pq = edge_map(edge);

        [dFpq, dFqp] = Gethdot(P0,edge);

        phi_h = phi_h1*abs(h(t,edge)*kappa(edge)*aeff(edge))^n + phi_h2*abs(h(t,edge)*kappa(edge)*aeff(edge))^m;
        dhpq = -1/phi1(edge)*(dFpq + phi0(edge) + phi_h);
        dhqp = -1/phi1(edge)*(dFqp + phi0(edge) + phi_h);
        
        aa1(t) = dFqp;
        if(dhpq < 0); dhpq = 0; end
        if(dhqp > 0); dhpq = -dhqp;
        elseif(dhqp < 0); dhqp = 0; end
        
        doth(edge) = dhpq; 
        bea = -2*abs(h(t,edge))/.3;
        if(bea == 0) 
            continue; 
        else 
            aeff(edge) = .5*(a(edge)*sqrt(bea^2*a(edge)^2 + 1) + asinh(bea*a(edge))/bea ); 
        end
    end
    
    [dv,dh,stp_e,stp_v] = UpdateVol(doth,dt);
    for i = 1:N-edges
        vol(i) = vol(i) + dv(i);
        V(t+1,i) = vol(i);
    end
    for i = 1:edges
        h(t+1,i) = h(t,i) + dh(i);
    end

    time(t+1) = t*dt;

    stress(t,:) = P0;
%     stress_shearMap(t,:) = Fstar(4,:);
    strain(t,:) = F0_t; 
    if(mod(t,100) == 0)
        fprintf('t = %f\n',t)
    end
end
%%

% figure(1)
% hold on
% plot(aa1,'LineWidth',1.3)
% plot([0 endt],[phi0(end) phi0(end)])
% plot([0 endt],[-phi0(end) -phi0(end)])
% xlabel('Iteration')
% ylabel('dAdh')
% legend('dadh','\phi_0^*','\phi_0^*')
% hold off 
% 
% figure(2)
% 
% hold on
% plot(time,h,'LineWidth',1.3)
% xlabel('Time')
% ylabel('Boundary displacment')
% grid on
% hold off
% 
% figure(3)
% hold on
% plot(time,V,'LineWidth',1.3)
% xlabel('Time')
% ylabel('Volume')
% 
% grid on
% hold off
figure(4)
hold on
% plot(abs(strain(:,1)-1),abs(stress(:,1)),'LineWidth',1.3)
% plot(strain(:,2),stress(:,2),'LineStyle','-.','LineWidth',1.3)
% plot(strain(:,3),stress(:,3),'LineStyle','-.','LineWidth',1.3)
plot(strain(:,4),stress(:,4),'LineStyle','-.','LineWidth',1.3)
% plot(strain(:,5)-1,stress(:,5),'LineWidth',1.3)
% plot(strain(:,6),stress(:,6),'LineStyle','-.','LineWidth',1.3)
% plot(strain(:,7),stress(:,7),'LineStyle','-.','LineWidth',1.3)
% plot(strain(:,8),stress(:,8),'LineStyle','-.','LineWidth',1.3)
% plot(strain(:,9)-1,stress(:,9),'LineWidth',1.3)
xlabel('Strain')
ylabel('Stress')
% legend('x_{11}','x_{21}','x_{31}','x_{12}','x_{22}','x_{23}','x_{31}','x_{32}','x_{33}','Location','northeastoutside')
grid on
hold off 

%%
np = pwd;
cd 'C:\NP_results'
fname = append('Results_test',num2str(randi([1 90000],1,1)));
fname = append(fname,append('_g',num2str(length(Grain_vol))));
fname = append(fname,notes);
np_data_dir = append('C:\NP_results\',fname);
mkdir(fname)
cd(fname)
save GBresults.mat stress strain V time h
save input.mat a Fpq ph_i1 Ph_i ph_i2 phi0 phi1 ss tt
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

function [vol, h,stp_e,stp_v] = UpdateVol(dh,dt)
global N edges Vtot C Cinv v Fgb F0  a D Aa I grain_map edge_map Neighbor_V2e

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
           neigh = Neighbor_V2e{p};
           stp_v(i) = true;
           for j = 1:length(neigh)
                stp_e(neigh(j)) = true;
           end
       end
       v(p) = v(p) + dV(i);
    end
    
    vol = -I*h; % volume of grains + quasi-grains
end

function [dhpq, dhqp] = Gethdot(P0,i)
    global N Vtot C Cinv  Fgb F0 a12 xi_inv edges alpha grain_map edge_map ss tt Neighbor_v Neighbor_e I 
    
    dhpq = 0; dhqp = 0; Ii = [1;0;0;0;1;0;0;0;1];
    pq = edge_map(i); pi = tt(i); qi = ss(i);

    XI = (Cinv(:,:,pi)-Cinv(:,:,qi))*P0; 
    Fstardh_ana_pq = (Ii - Fgb(:,:,pq));
    Fstardh_ana_pq = reshape(Fstardh_ana_pq,[3,3]);  
    XI = reshape(XI,[3,3]);

    P0 = reshape(P0,[3,3]);
   
    for k = 1:3
        for j = 1:3
            dhpq = dhpq + (1/2*XI(k,j) + Fstardh_ana_pq(k,j))*P0(k,j);
        end
    end
  dhqp = -dhpq;
end

function [Fstr] = FstarAnalytic(Fgb_avg,Cinv_ref)
    global N edges Vtot C Cinv v Fgb F0 xi_inv alpha grain_map
   
    Fstr = zeros(9,1,N-edges);
    alpha_temp = Vtot*F0 + Fgb_avg;

    for i = 1:N-edges
        index = grain_map(i);
        xi_p = Cinv_ref*C(:,:,index);

        alpha(:,:,i) = alpha_temp;
        xi_inv(:,:,i) = xi_p\eye(9);
        Fstr(:,:,i) = xi_inv(:,:,i)*alpha(:,:,i);
    end
   
end
function [C,Cinv] = SetC(p1,p,p2)

    % Copper elastic constants
    C11 = 169.3097; C12 = 122.5; C44 = 76; % cubic GPa
%     C11 = 169.3097; C12 = 87.2201 ; C44 = 41.0448; % Isotropic GPa

    % Aluminum elastic constants
%     C11 = 116.3; C12 = 64.8; C44 = 30.9; % cubic GPa

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
                   
                    C(I,J) = C12*d_ij*d_kl + 2*C44*(d_ik*d_jl) + C0*d_ijkl(I,J);
                    Cinv(I,J) = S12*d_ij*d_kl + 1/2*S44*(d_ik*d_jl) + S0*d_ijkl(I,J);

%                     C(I,J) = C12*d_ij*d_kl + C44*(d_ik*d_jl + d_il*d_jk) + C0*d_ijkl(I,J);
%                     Cinv(I,J) = S12*d_ij*d_kl + 1/4*S44*(d_ik*d_jl + d_il*d_jk) + S0*d_ijkl(I,J);

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
function [CC,CCinv,v,Fgb] = init(p1,p,p2,vol,Fpq,v_map,e_map)
    global N Vtot edges
   
    CC = zeros(9,9,N-edges);
    CCinv = CC; v = zeros(N,1); Fgb = zeros(3,3,N);
    for i = 1:N-edges
      index = v_map(i);
      [CC(:,:,index),CCinv(:,:,index)] = SetC(p1(i),p(i),p2(i));
        
      v(index) = vol(i);
      Fgb(:,:,index) = eye(3);
    end
    for i = 1:edges
        index = e_map(i);
        Fgb(:,:,index) = Fpq(:,:,i);
    end
end

function [C_ref, Cinv_ref, Fgb_avg] = GetAVG(vol)
    global N edges Vtot C Cinv v Fgb F0 xi_inv alpha grain_map edge_map
    
    Fgb_avg = zeros(9,1);
    C_ref = zeros(9,9); Cinv_ref = zeros(9,9);

    for p = 1:N-edges
%         C_ref = C_ref +  vol(p)*C(:,:,p); 
        Cinv_ref = Cinv_ref +  vol(p)*Cinv(:,:,p); 
    end
    for i = 1:edges
        p = edge_map(i);
        Fgb_avg = Fgb_avg - v(p)*(Fgb(:,:,p) - [1;0;0;0;1;0;0;0;1]);
    end
end
function F = CG(P, F0, C, mask)

    F = F0;
    r = P - C*F;
    for i = 1:9 
        if(mask(i)); r(i) = 0; end
    end
    if(sqrt(norm(r,2)) < 1e-7); return; end

    p = r;
    for i = 1:1000

        alpha = r'*r/(p'*C*p);
        F = F + alpha*p;
        rnew = r - alpha*C*p;

        for j = 1:9 
            if(mask(j)); rnew(j) = 0; end
        end
        if(sqrt(norm(rnew,2)) < 1e-7); return; end
        
        beta = rnew'*rnew/(r'*r);
        p = rnew + beta*p;
        r = rnew;
    end

end
