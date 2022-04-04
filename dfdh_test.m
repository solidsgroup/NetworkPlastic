clear all
global N Vtot C Cinv v Fgb F0 xi_inv phi0 alpha
Vtot = 1;
F0 = [1,.5 ,0;
      0, 1 ,0;
      0, 0 ,1];
F0 = reshape(F0, [9,1]);

grains = 100; 
num = 3:2:grains; num = [num]; 
dfdhpq = zeros(1,length(num)); dfdhqp = zeros(1,length(num));
phi0 = 0.3;
yield = zeros(1,length(num));
guess = 1; p0_l = 0;
edges = (num(end)-1)/2;
vol = rand(num(end)-edges,1); 
angle = (2*5*pi/12*rand(1,num(end))-5*pi/12)*.1; % -70->70 degrees


% vol = Vtot/(num(end)-edges)*ones(num(end)-edges,1);
for i = 1:length(num)
    N = num(i);
    
    fprintf('number of grains = %f, N = %f', (N+1)/2,N)

    Fgb = zeros(9,1,N);
    xi_inv = zeros(9,9,N);
%     angle = linspace(0,pi/2.6,N); % 85 degrees

    [C,Cinv,v,Fgbt] = init(angle,vol);
    v = Vtot*v./sum(v);
    alpha = zeros(9,1,N);
    for j = 1:N
        Fgb(:,:,j) = reshape(Fgbt(:,:,j), [9,1]);
    end
    counter = 1;
    for k = 2:2:N
        beta(counter,1) = Fgb(4,1,k); % check if least edge is alway given by the smallest beta
        counter = counter + 1;
    end
    least = inf; 
    for edges = 1:(N-1)/2
        shear = GoldenSectionSearch(0, 2,1e-6,edges);
        F0(4) = shear;
        Fstr = FstarAnalytic();
        P0 = C(:,:,1)*(Fstr(:,:,1)-Fgb(:,:,1));
        [dhpq, dhqp] = Gethdot(Fstr,P0, 0,edges);
        
    %     yield(i) = P0(4);
    %     dfdhpq(i) = dhpq;
    %     dfdhqp(i) = dhqp;
         if(abs(dhpq) - phi0 <= 2*1e-4)

            if(least >= shear) 
                least_ed = edges*2;
                least = shear;
                yield(i) = P0(4);
                dfdhpq(i) = dhpq;
                dfdhqp(i) = dhqp;
            end
            
         end
        
    end
    fprintf(', yield = %f, shear = %f, least edge = %f, beta = %f \n',yield(i), least, least_ed,Fgb(4,:,least_ed))    
end
expected_Yield = phi0*(N-1)/(2*sum(Fgb(4,:,:)));
%%
figure(1)
hold on
plot(num-1,dfdhpq)
plot(num-1,dfdhqp)
xlabel('Number of grains')
ylabel('value of dfdh')
% grid on
hold off
vol = (2*Vtot./(num-1)).^(-1/2);
figure(2)
hold on
set(gca, 'YScale', 'log', 'XScale', 'log')
plot(vol,yield,'.')
plot([vol(1) vol(end)],expected_Yield*ones(2,1))
grid on
xlabel('Volume (log(1/v^{1/2}))')
ylabel('Yield')
hold off

function dfdh = grain(gam,edge)
    global N Vtot C Cinv v Fgb F0 a12 xi_inv phi0
        F0(4) = gam;
        Fstr = FstarAnalytic();

        P0 = C(:,:,1)*(Fstr(:,:,1)-Fgb(:,:,1));
        [dhpq, dhqp] = Gethdot(Fstr,P0, 0,edge);

        dfdh = abs(dhpq) - phi0;

end
function [range] = GoldenSectionSearch(a, b,tol,edge)

golden = (sqrt(5)+1)/2;
c = b - (b-a)/golden;
d = a + (b-a)/golden;
while abs(a-b) >= 2*tol
    g_c = -grain(c,edge);
    g_d = grain(d,edge);
    if g_c < g_d 
        b = d;
    else
        a = c;
    end
    c = b - (b-a)/golden;
    d = a + (b-a)/golden;
end
range = (a+b)/2;
end

function [dhpq, dhqp] = Gethdot(Fstr,P0, vm,i)
    global N Vtot C Cinv v Fgb F0 a12 xi_inv alpha
    
    edges = (N-1)/2;
    Fstardh_ana = zeros(3,3,N-edges);
    dhpq = 0; dhqp = 0; I = [1;0;0;0;1;0;0;0;1];
    vol = v(1:2:end);
    pq = i*2;
    counter = 1;
    
    for pp = 1:edges
            
        m = pp*2; q = m - 1; p = m + 1;
        
        C_diff = eye(9) - Cinv(:,:,p)*C(:,:,q); 
        dfstar_dhpq = -xi_inv(:,:,q)^2*(C_diff)*alpha(:,:,q) + xi_inv(:,:,q)*(I - Fgb(:,:,pq) - C_diff*Fgb(:,:,q));

        Fstardh_ana(:,:,counter) = reshape(dfstar_dhpq,[3,3]);
        counter = counter + 1;
            
    end

    if(mod(N,2) == 1)
        
        C_diff = eye(9) - Cinv(:,:,q)*C(:,:,p);
        dfstar_dhpq = -xi_inv(:,:,p)^2*(C_diff)*alpha(:,:,p) + xi_inv(:,:,p)*(I - Fgb(:,:,pq) - C_diff*Fgb(:,:,p));

        Fstardh_ana(:,:,end) = reshape(dfstar_dhpq,[3,3]);

    end
    
    counter = 1; Fstardh = zeros(3,3);
    for j = 1:2:(N-edges)*2
        
        Fstardh = Fstardh + Fstardh_ana(:,:,counter)*vol(counter);
        counter = counter + 1;
    end

    P0 = reshape(P0,[3,3]);
   
    for i = 1:3
        for j = 1:3
%             dhpq = dhpq + (vp/Vtot*Fstardh_p_ana(i,j) + vq/Vtot*Fstardh_q_ana(i,j))*P0(i,j);
            dhpq = dhpq + Fstardh(i,j)*P0(i,j);
%             dhpq = dhpq +  Fstardh_q(i,j)*P(i,j);
        end
    end
   
    for i = 1:3
        for j = 1:3
          dhqp = dhqp + Fstardh(i,j)*P0(i,j);
%             dhqp = dhqp + Fstardh_p(i,j)*P(i,j);
        end
    end
   
end

function [CC,CCinv,v,Fgb] = init(angle,vol)
    global N Vtot
   
    CC = zeros(9,9,N);
    CCinv = CC; v = zeros(N,1); Fgb = zeros(3,3,N);
%     
%     edges = (N-1)/2;
%     vol = Vtot/(N-edges);
    counter = 1;
    for i = 1:2:N
       
      [CC(:,:,i),CCinv(:,:,i)] = SetC(angle(i));
        
      v(i) = vol(counter);
      Fgb(:,:,i) = eye(3);
     counter = counter + 1;
    end
    for i = 2:2:N
        
      Fgb(:,:,i) = eye(3);     
      CC(:,:,i) = CC(:,:,i-1);
      CCinv(:,:,i) = CCinv(:,:,i-1);
      Fgb(1,2,i) = 2*tan(angle(i-1)/2);
      v(i) = 0;

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
function [Fstr] = FstarAnalytic()
    global N Vtot C Cinv v Fgb F0 xi_inv alpha
   
    Fstr = zeros(9,1,(N-1)/2);
    
    for i = 1:2:N
        xi = v(i)*eye(9);
        alpha_temp = Vtot*F0;
        for j = 1:2:N            
            if(i == j)
                continue
            end        
            xi = xi + v(j)*Cinv(:,:,j)*C(:,:,i);
            alpha_temp = alpha_temp - v(j)*Fgb(:,:,j);
            
        end
        alpha_temp = alpha_temp + (xi-v(i)*eye(9))*Fgb(:,:,i);
        alpha(:,:,i) = alpha_temp + (xi-v(i)*eye(9))*Fgb(:,:,i);
        xi_inv(:,:,i) = xi\eye(9);
    
        Fstr(:,:,i) = xi_inv(:,:,i)*alpha_temp;
  
    end
   
end