clear all
global N Vtot C Cinv v Fgb F0 edges D Aa I grain_map edge_map ss tt a Neighbor_V2e

loadPath = 'C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\NP_stat_varOri'; notes = 'Var_ori';
% loadPath = 'C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\EngEBSD'; notes = 'EBSD_eng'; % 1-9, 10-18, 19-27
% loadPath = 'C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\NonEngEBSD'; notes = 'EBSD_nonEng'; % 1-8, 9-16, 17-24
% loadPath = 'C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\NP_stat_varMico'; notes = 'Var_micro';
cd(loadPath)
% 1-12,13-24,25-36,37-48,49-61,62-74,75-87,88-100 
for iter = 88:100
    clear N Vtot C Cinv v Fgb F0 edges D Aa I grain_map edge_map ss tt a Neighbor_V2e
    global N Vtot C Cinv v Fgb F0 edges D Aa I grain_map edge_map ss tt a Neighbor_V2e
    fname = 'microstructureData_GBengSubset_';
    fname = append(fname,num2str(iter));
    load(fname) 
    
    N = length(s)+length(Grain_vol); edges = length(s); 
    Grain_vol = Grain_vol./sum(Grain_vol); a = a./sum(a);
    Vtot = sum(Grain_vol);
    
    grain_map = 1:N-edges;
    edge_map = N-edges+1:N; 
    
    F0 = [1,0 ,0;
          0, 1 ,0;
          0, 0 ,1];
    Fgb = zeros(9,1,N);
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
    I = incidence(D)*spdiags(a,0,length(a),length(a));
    
    stp_e = zeros(edges,1); flag_edge = zeros(edges,1);
    stp_v = zeros(N-edges,1);
    stress_shearMap = zeros(endt,N-edges);
    dhpq = 0; dhqp = 0; doth = zeros(edges,1);
    h = zeros(1,edges);
    V = zeros(1,N-edges);
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
%         if(Flag(tt(i)))
%             gg = .05;
%         elseif(Flag(ss(i)))
%             gg = .05;
%         else
%             gg = .00005; 
%         end
%         gg = .005; 
%         gg = 5e-3;  
        gg = .00005;
        kappa(i) = gg/(a(i)*min(vv));
    end
    aeff = a; 
    phi1 = phi0*0.001;
    
    Ii = [1;0;0;0;1;0;0;0;1];dfdhpq = 0; dfdhqp = 0;
    C_ref = zeros(9,9); Cinv_ref = zeros(9,9); mask = [1;0;0;0;0;0;0;0;0];
    P0 = zeros(9,1); F0 = zeros(9,1); Fgb_avg = zeros(9,1);
    for t = 1:endt
        
        [C_ref, Cinv_ref, Fgb_avg] = GetAVG(vol);

        % _____Stress driven motion_____
    %     P0 = [10*t/endt;0;0;10*t/endt*0;0;0;0;0;0]; 
    %     F0 = Cinv_ref*P0 + Ii - Fgb_avg;

        % _____Strain driven motion_____
        T = (Cinv_ref)\eye(9);
        F0 = [-.05*t/endt;0;0;0;0;0;0;0;0];
        P = -T*Fgb_avg;
        F0 = CG(P, F0, T, mask);
        F0 = F0 + Ii;
        
        P0 = T*(Vtot*F0 + Fgb_avg - Ii);
        for edge = 1:edges

            if(stp_e(edge)) % stop gb motion to conserve volume
                doth(edge) = 0;
                continue
            end
            q = grain_map(ss(edge)); p = grain_map(tt(edge));  pq = edge_map(edge);
    
            [dFpq, dFqp] = Gethdot(P0,edge);
    
            phi_h = phi_h1*abs(h(1,edge)*kappa(edge)*aeff(edge))^n + phi_h2*abs(h(1,edge)*kappa(edge)*aeff(edge))^m;
            dhpq = -1/phi1(edge)*(dFpq + phi0(edge) + phi_h);
            dhqp = -1/phi1(edge)*(dFqp + phi0(edge) + phi_h);
            
            aa1(t) = dFqp;
            if(dhpq < 0); dhpq = 0; end
            if(dhqp > 0); dhpq = -dhqp;
            elseif(dhqp < 0); dhqp = 0; end
    
            doth(edge) = dhpq; 
            bea = -2*abs(h(1,edge))/.3;
            if(bea == 0) 
                continue; 
            else 
                aeff(edge) = .5*(a(edge)*sqrt(bea^2*a(edge)^2 + 1) + asinh(bea*a(edge))/bea ); 
            end
    
        end
        
        [dv,dh,stp_e,stp_v] = UpdateVol(doth,dt);
        for i = 1:N-edges
            vol(i) = vol(i) + dv(i);
            V(1,i) = vol(i);
        end
        for i = 1:edges
            h(1,i) = h(1,i) + dh(i);
        end
    
        time(t+1) = t*dt;
     
        stress(t,:) = P0;
        strain(t,:) = F0; 
        if(mod(t,100) == 0)
            fprintf('t = %f, iter = %f, g = %f\n',t,iter, N-edges)
        end
    end
    
    np = pwd; file = append('C:\NP_results\NP_Stats\', notes);
    cd(file)
    fname = append('Results_engSubset_',num2str(iter));
%     fname = append(fname,append('_g',num2str(length(Grain_vol))));
    fname = append(fname,notes); fname = append(fname,'.mat');
    save(fname, 'stress', 'strain', 'time','phi0','phi1','ph_i1','Ph_i','ph_i2') 

    cd(loadPath)
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
    global N Vtot C Cinv  Fgb F0 edges grain_map edge_map ss tt I 
    
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
    global N edges Vtot C Cinv v Fgb F0 grain_map edge_map
    
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