%% For single files
% clear all
% filename = 'C:\Users\icrma\source\repos\test\test\Microstructure.txt';
% delimiterIn = ' ';
% headerlinesIn = 1;
% A = importdata(filename,delimiterIn,headerlinesIn);
% V = A.data(1,1); % number of grains
% bimodal = A.data(1,2);
% numSec = A.data(1,3);
% s = A.data(2:end,1) + 1;
% t = A.data(2:end,2) + 1;
% a = A.data(2:end,3);
% th_gb = A.data(2:end,4);
% Grain_vol = A.data(2:V+1,5);
% 
% edges = length(s);
% 
% % rotation of grains
% ph_i1 = randn(V,1);
% Ph_i = randn(V,1);
% ph_i2 = randn(V,1);
% % calculate Fpq and phi0
% for edge = 1:length(s)
%     p = s(edge);
%     q = t(edge);
%     Rgb = [cos(th_gb(edge)) -sin(th_gb(edge)) 0;
%     sin(th_gb(edge)) cos(th_gb(edge)) 0;
%     0 0 1];
% 
%     g1 = [ph_i1(p) Ph_i(p) ph_i2(p)];
%     g2 = [ph_i1(q) Ph_i(q) ph_i2(q)];
% 
%     r11 = cos(g1(1))*cos(g1(3)) - cos(g1(2))*sin(g1(1))*sin(g1(3));
%     r12 = -cos(g1(1))*sin(g1(3)) - cos(g1(2))*sin(g1(1))*cos(g1(3));
%     r13 = sin(g1(2))*sin(g1(1));
%     r21 = cos(g1(3))*sin(g1(1)) + cos(g1(2))*cos(g1(1))*sin(g1(3));
%     r22 = cos(g1(2))*cos(g1(1))*cos(g1(3)) -  sin(g1(1))*sin(g1(3));
%     r23 = -sin(g1(2))*cos(g1(1));
%     r31 = sin(g1(2))*sin(g1(3));
%     r32 = sin(g1(2))*cos(g1(3));
%     r33 = cos(g1(2));
% 
%     rot1 = [r11 r12 r13;
%     r21 r22 r23;
%     r31 r32 r33];
%     r11 = cos(g2(1))*cos(g2(3)) - cos(g2(2))*sin(g2(1))*sin(g2(3));
%     r12 = -cos(g2(1))*sin(g2(3)) - cos(g2(2))*sin(g2(1))*cos(g2(3));
%     r13 = sin(g2(2))*sin(g2(1));
%     r21 = cos(g2(3))*sin(g2(1)) + cos(g2(2))*cos(g2(1))*sin(g2(3));
%     r22 = cos(g2(2))*cos(g2(1))*cos(g2(3)) -  sin(g2(1))*sin(g2(3));
%     r23 = -sin(g2(2))*cos(g2(1));
%     r31 = sin(g2(2))*sin(g2(3));
%     r32 = sin(g2(2))*cos(g2(3));
%     r33 = cos(g2(2));
%     rot2 = [r11 r12 r13;
%         r21 r22 r23;
%         r31 r32 r33];
%     tt = rot2*rot1';
%     axisangle = zeros(2,1);
%     [G,H] = eig(tt);
%     for j = 1:3
%         if(isreal(G(:,j)))
%             axisangle(1) = G(1,j);
%             axisangle(2) = G(2,j);
%             axisangle(3) = G(3,j);
%             break
%         end
%     end
%     a_dir = Rgb*[1;0;0];
%     phia = atan2(norm(cross(a_dir,axisangle)),dot(a_dir,axisangle));
% %     th = acos((trace(tt)-1)/2); % symmetric rotation angle
%     th = atan2(norm(axisangle),(trace(tt)-1)/2);
%    
%     if(phia < pi/18 && phia > -pi/18 )
%         beta = 0;
%         atgb(edge) = true;
%     else
%         beta = GetBeta(th);
%         atgb(edge) = false;
%     end
%     if(beta > 4)
%         beta = 0;
%         atgb(edge) = true;
%     end
%     shearCouling(edge) = beta;
%     Fpq_temp = eye(3); Fpq_temp(1,2) = beta;
%     Fpq(:,:,edge) = Fpq_temp;
%     Fpq(:,:,edge) = Rgb'*Fpq(:,:,edge)*Rgb;
% 
%     phi0(edge) = GetPhi0(phia);
% end
% mean(shearCouling)
% phi1 = phi0*0.1;
% fname = 'microstructureData';
% if(bimodal)
%     fname = append(fname,append('_Bimodal_N',num2str(numSec)));
% else
%     fname = append(fname,append('_',num2str(V)));
% end
% fname = append(fname,'.mat');
% save(fname, 'Fpq', 'a', 'Grain_vol', 's', 't', 'ph_i1', 'Ph_i', 'ph_i2', 'phi1', 'phi0')  
%% For grain orientation

clear all
filename = 'C:\Users\icrma\source\repos\test\test\NP_stat\microstructure_0';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);
V = A.data(1,1); % number of grains
bimodal = A.data(1,2);
numSec = A.data(1,3);
s = A.data(2:end,1) + 1;
t = A.data(2:end,2) + 1;
a = A.data(2:end,3);
th_gb = A.data(2:end,4);
Grain_vol = A.data(2:V+1,5);
Flag = A.data(2:V+1,6);

edges = length(s);
cd 'C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\NP_stat_varOri'

% load C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\optimalAngles.mat
% ph_i1 = opt_angles(:,1);
% Ph_i = opt_angles(:,2);
% ph_i2 = opt_angles(:,3);
for i = 1:100
% rotation of grains -- STD and Mean values from GB engineered EBSD data

% Ph_i = 0.1833*randn(15000,1) + 0.578;
% ph_i1 = 1.8178*randn(15000,1) + 3.1165;
% ph_i2 = mod(0.3181*randn(15000,1) + .8361,2*pi); %0.5681*randn(15000,1) + 0.8361;


% Ph_i = mod(0.1963*randn(15000,1) + 0.5615,2*pi);
% ph_i1 = mod(2.5784*randn(15000,1) + 3.4146,2*pi);
% ph_i2 = mod(0.3051*randn(15000,1) + 0.8092,2*pi);

% rotation of grains -- GB engineered subset

Ph_i = mod(0.1748*randn(15000,1) + 0.5901,2*pi);
ph_i1 = mod(1.9557*randn(15000,1) + 1.2191,2*pi);
ph_i2 = mod(0.3081*randn(15000,1) + 0.8226,2*pi);

% rotation of grains -- STD and Mean values from nonGB engineered EBSD data

%     Ph_i = 0.1912*randn(15000,1) + 0.5735;
%     ph_i1 = 1.8017*randn(15000,1) + 3.1006;
%     ph_i2 = mod(0.3551*randn(15000,1) + 0.8592,2*pi); % 0.7151*randn(15000,1) + 0.8592;

%     Ph_i = mod(0.2066*randn(15000,1) + 0.5573,2*pi);
%     ph_i1 = mod(2.5857*randn(15000,1) + 3.1085,2*pi);
%     ph_i2 = mod(0.3081*randn(15000,1) + 0.7858,2*pi);
    
    % calculate Fpq and phi0
    for edge = 1:length(s)
        p = s(edge);
        q = t(edge);
        Rgb = [cos(th_gb(edge)) -sin(th_gb(edge)) 0;
        sin(th_gb(edge)) cos(th_gb(edge)) 0;
        0 0 1];
    
        g1 = [ph_i1(p) Ph_i(p) ph_i2(p)];
        g2 = [ph_i1(q) Ph_i(q) ph_i2(q)];
    
        r11 = cos(g1(1))*cos(g1(3)) - cos(g1(2))*sin(g1(1))*sin(g1(3));
        r12 = -cos(g1(1))*sin(g1(3)) - cos(g1(2))*sin(g1(1))*cos(g1(3));
        r13 = sin(g1(2))*sin(g1(1));
        r21 = cos(g1(3))*sin(g1(1)) + cos(g1(2))*cos(g1(1))*sin(g1(3));
        r22 = cos(g1(2))*cos(g1(1))*cos(g1(3)) -  sin(g1(1))*sin(g1(3));
        r23 = -sin(g1(2))*cos(g1(1));
        r31 = sin(g1(2))*sin(g1(3));
        r32 = sin(g1(2))*cos(g1(3));
        r33 = cos(g1(2));
    
        rot1 = [r11 r12 r13;
        r21 r22 r23;
        r31 r32 r33];
        r11 = cos(g2(1))*cos(g2(3)) - cos(g2(2))*sin(g2(1))*sin(g2(3));
        r12 = -cos(g2(1))*sin(g2(3)) - cos(g2(2))*sin(g2(1))*cos(g2(3));
        r13 = sin(g2(2))*sin(g2(1));
        r21 = cos(g2(3))*sin(g2(1)) + cos(g2(2))*cos(g2(1))*sin(g2(3));
        r22 = cos(g2(2))*cos(g2(1))*cos(g2(3)) -  sin(g2(1))*sin(g2(3));
        r23 = -sin(g2(2))*cos(g2(1));
        r31 = sin(g2(2))*sin(g2(3));
        r32 = sin(g2(2))*cos(g2(3));
        r33 = cos(g2(2));
        rot2 = [r11 r12 r13;
            r21 r22 r23;
            r31 r32 r33];
        tt = rot2*rot1';
        axisangle = zeros(2,1);
        [G,H] = eig(tt);
        for j = 1:3
            if(isreal(G(:,j)))
                axisangle(1) = G(1,j);
                axisangle(2) = G(2,j);
                axisangle(3) = G(3,j);
                break
            end
        end
        a_dir = Rgb*[1;0;0];
        phia = atan2(norm(cross(a_dir,axisangle)),dot(a_dir,axisangle));
    %     th = acos((trace(tt)-1)/2); % symmetric rotation angle
        th = atan2(norm(axisangle),(trace(tt)-1)/2);
       
        if(phia < pi/18 && phia > -pi/18 )
            beta = 0;
            atgb(edge) = true;
        else
            beta = GetBeta(th);
            atgb(edge) = false;
        end
        if(beta > 2)
            beta = 0;
            atgb(edge) = true;
        end
        shearCouling(edge) = beta;
        Fpq_temp = eye(3); Fpq_temp(1,2) = beta;
        Fpq(:,:,edge) = Fpq_temp;
        Fpq(:,:,edge) = Rgb'*Fpq(:,:,edge)*Rgb;
    
        phi0(edge) = GetPhi0(phia);
    end
    fprintf('ori iter = %f\n', i)
    phi1 = phi0*0.1;
    fname = 'microstructureData_GBengSubset_'; fname = append(fname,num2str(i));
    if(bimodal)
%         fname = append(fname,append('_Bimodal_N',num2str(numSec)));
    else
%         fname = append(fname,append('_',num2str(V)));
    end
    fname = append(fname,'.mat');
    save(fname, 'Fpq', 'a', 'Grain_vol', 's', 't', 'ph_i1', 'Ph_i', 'ph_i2', 'phi1', 'phi0','Flag')  
end
%% For varying microstructure

clear all
cd 'C:\Users\icrma\OneDrive\Documents\MATLAB\NetworkPlastic\NP_stat_varMico'

% rotation of grains -- STD and Mean values from GB engineered EBSD data

% Ph_i = 0.1833*randn(15000,1) + 0.578;
% ph_i1 = 1.8178*randn(15000,1) + 3.1165;
% ph_i2 = mod(0.3181*randn(15000,1) + .8361,2*pi); %0.5681*randn(15000,1) + 0.8361;

% rotation of grains -- STD and Mean values from nonGB engineered EBSD data

    ph_i1 = 1.8017*randn(15000,1) + 3.1006;
    Ph_i = 0.1912*randn(15000,1) + 0.5735;
    ph_i2 = mod(0.3551*randn(15000,1) + 0.8592,2*pi); % 0.7151*randn(15000,1) + 0.8592;

for i = 1:100
    filename = 'C:\Users\icrma\source\repos\test\test\NP_stat\microstructure_';
    filename = append(filename,num2str(i-1));
    delimiterIn = ' ';
    headerlinesIn = 1;
    A = importdata(filename,delimiterIn,headerlinesIn);
    V = A.data(1,1); % number of grains
    bimodal = A.data(1,2);
    numSec = A.data(1,3);
    s = A.data(2:end,1) + 1;
    t = A.data(2:end,2) + 1;
    a = A.data(2:end,3);
    th_gb = A.data(2:end,4);
    Grain_vol = A.data(2:V+1,5);
    
    D = digraph(s,t,a);
    Aa = sparse(adjacency(D,a));

    edges = length(s);
    % rotation of grains -- STD and Mean values from GB engineered EBSD data
    
    % calculate Fpq and phi0
    for edge = 1:length(s)
        p = s(edge);
        q = t(edge);
        Rgb = [cos(th_gb(edge)) -sin(th_gb(edge)) 0;
        sin(th_gb(edge)) cos(th_gb(edge)) 0;
        0 0 1];
    
        g1 = [ph_i1(p) Ph_i(p) ph_i2(p)];
        g2 = [ph_i1(q) Ph_i(q) ph_i2(q)];
    
        r11 = cos(g1(1))*cos(g1(3)) - cos(g1(2))*sin(g1(1))*sin(g1(3));
        r12 = -cos(g1(1))*sin(g1(3)) - cos(g1(2))*sin(g1(1))*cos(g1(3));
        r13 = sin(g1(2))*sin(g1(1));
        r21 = cos(g1(3))*sin(g1(1)) + cos(g1(2))*cos(g1(1))*sin(g1(3));
        r22 = cos(g1(2))*cos(g1(1))*cos(g1(3)) -  sin(g1(1))*sin(g1(3));
        r23 = -sin(g1(2))*cos(g1(1));
        r31 = sin(g1(2))*sin(g1(3));
        r32 = sin(g1(2))*cos(g1(3));
        r33 = cos(g1(2));
    
        rot1 = [r11 r12 r13;
        r21 r22 r23;
        r31 r32 r33];
        r11 = cos(g2(1))*cos(g2(3)) - cos(g2(2))*sin(g2(1))*sin(g2(3));
        r12 = -cos(g2(1))*sin(g2(3)) - cos(g2(2))*sin(g2(1))*cos(g2(3));
        r13 = sin(g2(2))*sin(g2(1));
        r21 = cos(g2(3))*sin(g2(1)) + cos(g2(2))*cos(g2(1))*sin(g2(3));
        r22 = cos(g2(2))*cos(g2(1))*cos(g2(3)) -  sin(g2(1))*sin(g2(3));
        r23 = -sin(g2(2))*cos(g2(1));
        r31 = sin(g2(2))*sin(g2(3));
        r32 = sin(g2(2))*cos(g2(3));
        r33 = cos(g2(2));
        rot2 = [r11 r12 r13;
            r21 r22 r23;
            r31 r32 r33];
        tt = rot2*rot1';
        axisangle = zeros(2,1);
        [G,H] = eig(tt);
        for j = 1:3
            if(isreal(G(:,j)))
                axisangle(1) = G(1,j);
                axisangle(2) = G(2,j);
                axisangle(3) = G(3,j);
                break
            end
        end
        a_dir = Rgb*[1;0;0];
        phia = atan2(norm(cross(a_dir,axisangle)),dot(a_dir,axisangle));
    %     th = acos((trace(tt)-1)/2); % symmetric rotation angle
        th = atan2(norm(axisangle),(trace(tt)-1)/2);
       
        if(phia < pi/18 && phia > -pi/18 )
            beta = 0;
            atgb(edge) = true;
        else
            beta = GetBeta(th);
            atgb(edge) = false;
        end
        if(beta > 4)
            beta = 0;
            atgb(edge) = true;
        end
        shearCouling(edge) = beta;
        Fpq_temp = eye(3); Fpq_temp(1,2) = beta;
        Fpq(:,:,edge) = Fpq_temp;
        Fpq(:,:,edge) = Rgb'*Fpq(:,:,edge)*Rgb;
    
        phi0(edge) = GetPhi0(phia);
        if(isnan(phi0(edge)))
            error('phi0 nan')
        end
    end
    fprintf('micro iter = %f\n', i)
    phi1 = phi0*0.1;
    fname = 'microstructureData_non_'; fname = append(fname,num2str(i));
    if(bimodal)
%         fname = append(fname,append('_Bimodal_N',num2str(numSec)));
    else
%         fname = append(fname,append('_',num2str(V)));
    end
    fname = append(fname,'.mat');
    save(fname, 'Fpq', 'a', 'Grain_vol', 's', 't', 'ph_i1', 'Ph_i', 'ph_i2', 'phi1', 'phi0')  
end

%% ------------------------------------------------------------
function phi0 = GetPhi0(th)
    
    angle = mod(th,pi/2);
%     if(angle < 0.2443)
%         phi0 = 5.279396826562586*angle + 0.01;
%     elseif(angle > 1.4835)
%         phi0 = -1.031324031235482*angle + 1.63;
%     else
        tilt_angle = [0        14 19 22 25 28 31 34 38 50.5 53 57 59.5 63 64 68 71 74 77 79.5 81 83 85               90]*pi/180;
        disp_energy =[0.1      1.3 1.6 1.25 .98 .8 .7 .65 .85 .4 2.8 .3 .35 .45 .48 .55 .48 .5 .4 .25 .23 .15 .1     0.1]; % GPa; 
        phi0 = interpn(tilt_angle,disp_energy,angle,'makima'); 
    end
function beta = GetBeta(th)

 % for Al
%   if(th <= pi/6)
%       mu = 1;
%   elseif(th >= pi/6 && th < pi/3)
%       mu = 3/4;
%   else
%       mu = 1/2;
%   end
  % for Cu
  if(th <= 41*pi/200)
      mu = 1;
  elseif(th >= 41*pi/200 && th < 0.881)
      mu = 3/4;
  else
      mu = 1/2;
  end
  beta = (2*mu*tan(th/2) - 2*(1-mu)*tan(pi/2-th/2));

end