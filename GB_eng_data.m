clear all
fname = fullfile(mtexDataPath,'EBSD','Cu_GBE_EBSD_data','8july_1_1-subset1.txt'); % 8july_1_1-subset1.txt 24june_1_0.txt
CS = { loadCIF('Cu-Copper.cif')};
ebsd = loadEBSD_generic(fname,'CS',CS, 'ColumnNames', ...
{'Index' 'Phase' 'x' 'y' 'Euler1' 'Euler2' 'Euler3' 'MAD' 'OI' 'BC' 'BS'...
'Status'}, 'Bunge');
threshold = 5 *degree;
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',threshold);
ebsd(grains(grains.grainSize < 50)) = [];
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',threshold);
grains = smooth(grains,20);


% % set up the ipf coloring
% cKey = ipfColorKey(ebsd.CS.properGroup);
% cKey.inversePoleFigureDirection = xvector;
% color = cKey.orientation2color(ebsd.orientations);
% GB = grains.boundary;
% oriGB = grains.boundary.misorientation;
% 
% plot(grains,grains.meanOrientation,'faceAlpha',0.8,'figSize','large')
% plot(ebsd,ebsd.orientations)
% init override mode
% hold on
% 
% plot(GB,angle(oriGB)./degree,'linewidth',2)
% mtexColorMap blue2red
% mtexColorbar
% 
% hold off
p1 = [0,95]; p2 = [95 95]; p3 = [95,0];

% transfer this into degree and project it into the interval [0,180]
dir = 1:length(grains);

% plot the direction
% plot(grains,dir,'micronbar','off')
% mtexColorbar
plot(ebsd,ebsd.orientations)
hold on
plot(grains.boundary,'lineWidth',2)
% plot([p1(1),p2(1)],[p1(2),p2(2)],'-.b','lineWidth',2)
% plot([p2(1),p3(1)],[p2(2),p3(2)],'-.b','lineWidth',2)
% hold off
% 
% restrict to iron to iron phase transition
gB = grains.boundary;

% select CSL(3) grain boundaries
gB3 = gB(angle(gB.misorientation,CSL(3,ebsd.CS)) < 3*degree);

% figure(3)
% plot(ebsd,log(ebsd.prop.iq),'figSize','large')
% mtexColorMap black2white
% CLim(gcm,[.5,5])
% 
% the option 'FaceAlpha',0.4 makes the plot a bit transluent
% hold on
% plot(grains,grains.meanOrientation,'FaceAlpha',0.4,'linewidth',3)
% hold off
% hold on
% plot(gB3,'lineColor','gold','linewidth',3,'DisplayName','CSL 3')
% hold off

sum(gB3.segLength)/sum(gB.segLength)
%%
%Grain_id = [60 82 81 91 74 46 20 16 50 43 49 52 66 67 47 39 30 63 18 35 32 24 31 34 44 64 61 75 87 118 104 95 115 119 112 103 71 84 89 92 99 116 114 117 166 174 147 187 185 190 146 232 155 173 170 171 127 153 180 159 145 165];

% use this: ebsd(grains(grains.grainSize < 50)) = []
% Grain_id = [60 81 85 90 116 129 184 189 147 187 174 166 190 146 117 114 99 91 92 46 20 16 1 47 43 50 12 4 18 17 35 24 32 31 63 39 30 52 49 66 67 71 103 92 89 232 153 127 145 159 165 171 170 155 119 95 104 75 61 82 115 112 87 84 74        185 180 173 217 64 44 29 22 5 34];
% Grain_id = unique(Grain_id);

% 
% s = [1 1 1 2 2 2 2 2 3 4 4 5 6 6 7 7 7 8 9 10 10 11 11 11 12 12 13 13 13 14 14 15 15 15 15 15 15 15 16 16 17 18 18 18 18 19 20 20 21 21 22 23 23 24 24 25 26 27 29 29 29 30 31 32 33 33 33 33 34 34 34 35 35 36 36 36 37 37 38 38 39 39 39 39 40 41 41 41 42 42 42 43 43 43 44 45 45 45 45 46 46 47 48 48 49 49 49 50 50 51 52 52 52 52 52 52 52 52 53 53 53 54 54 55 55 56 56 57 58 58 58 59 59 60 60 60 60 61 61 62]';
% t = [2 3 7 7 10 5 13 15 2 3 2 38 9 2 8 15 15 15 10 15 12 10 9 12 13 37 9 14 11 12 37 16 17 18 18 19 20 37 17 18 18 37 25 26 27 20 21 22 23 18 23 24 18 25 18 26 29 28 31 27 33 29 30 31 28 18 34 32 32 31 52 33 36 18 33 52 10 39 13 41 14 40 36 52 36 39 4 52 4 41 43 41 44 45 52 44 42 46 52 42 49 46 47 46 48 45 50 51 45 45 35 40 33 43 51 57 56 55 51 54 50 50 55 53 56 59 58 58 52 52 60 60 58 61 62 52 52 52 58 52]';
Grain_id = 1:length(grains);
% figure(1)
% plot(grains(Grain_id))

neighbors = grains.neighbors;

[s,t] = MakeGraph(Grain_id, neighbors);

a = zeros(length(s),1); scale = 1e-6;
a_dir = zeros(length(s),2);
shift = 0; Edge_delete = zeros(length(s),1);
csl3GB = zeros(length(s),1); cslarea = []; noncslarea = [];
for edge = 1:length(s)
    
    sp = [Grain_id(t(edge)) Grain_id(s(edge)) ];
    gbs =  grains.boundary(all(grains.boundary.grainId == [sp],2));

    [empty, uuu] = size(gbs);
    if(empty == 0)
        sp = [Grain_id(s(edge)) Grain_id(t(edge)) ];
        gbs =  grains.boundary(all(grains.boundary.grainId == [sp],2));
    end
    gB3 = gbs(angle(gbs.misorientation,CSL(3,ebsd.CS)) < 3*degree);
    [empty, uuu] = size(gB3);
    if(empty ~= 0)
        csl3GB(edge) = true;
        cslarea(end+1) = edge;
    else
        noncslarea(end+1) = edge;
    end
    gbs_split = gbs;counter = 0;
%     for i = 1:length(gbs) % cut off grain boundaries
%         
%         x = gbs(i).x;
%         y = gbs(i).y;
%         if(any(abs(x) > abs(p1(2))) )
%             gbs_split(i-counter) = [];
%             counter = counter + 1;
%         elseif(any(abs(y) > abs(p1(2))) )
%             gbs_split(i-counter) = [];
%             counter = counter + 1;
%         end
%     end
   
    [empty, uuu] = size(gbs_split);
    pts = gbs.triplePoints.V;
    if(empty == 0)
%         a = min(pts)*scale;
        Edge_delete(edge) = true;
        continue
    end

    a(edge) = sum(gbs_split.segLength)*scale;
    [row, uuu] = size(pts);
    if(edge == 60)
        aaa = 0;
    end
    if(row == 1)
        pts2 = [gbs.x,gbs.y];
        a_dir(edge,:) = [(pts2(1,1)-pts2(end,1)), (pts2(1,2)-pts2(end,2))];
        a_dir(edge,:) = a_dir(edge,:)./norm(a_dir(edge,:));
    elseif(row == 0) % grain is circular 
%         s(edge) = [];
%         t(edge) = [];
%         a(edge) = [];
%         a_dir(edge,:) = [];
%         shift = shift + 1;
%         continue
        pts2 = [gbs.x,gbs.y];
        pts2 = diff(pts2);
        pts2 = mean(pts2);
        a_dir(edge,:) = null(pts2);
        a_dir(edge,:) = a_dir(edge,:)./norm(a_dir(edge,:));
    else
        pts_fin = max(pts); pts_st = min(pts);
        a_dir(edge,:) = [(pts_st(1)-pts_fin(1)), (pts_st(2)-pts_fin(2))];
        a_dir(edge,:) = a_dir(edge,:)./norm(a_dir(edge,:));
    end
    ori2 = grains(sp(2)).meanOrientation;
    ori1 = grains(sp(1)).meanOrientation;
    missorientation(edge) = angle(ori1, ori2); % rotation around axis angle

    mori(edge) = axis(ori1,ori2); % axis angle
    
end
% sp = [Grain_id(t(1)) Grain_id(s(1)) ]
% gbs =  grains.boundary(all(grains.boundary.grainId == [sp],2));

% i = 1; len = length(s); Edge_delete(110) = true; Edge_delete(111) = true;
% shift_grain = 0; shift_edge = 0;
% while i <= len 
% 
%     if(Edge_delete(i))
%         s(i-shift_edge) = [];
%         t(i-shift_edge) = [];
%         a(i-shift_edge) = [];
%         a_dir(i-shift_edge,:) = [];
%         ph1(i-shift_edge) = [];
%         Ph(i-shift_edge) = [];
%         ph2(i-shift_edge) = [];
%         
%         shift_edge = shift_edge +1;
%         for j = 2:length(s)
%             if(s(j-1) - s(j) < -1 )
%                 Grain_id(s(j) + s(j-1)) = [];
%                 shift_grain = shift_grain + 1;
%             end
%         end
%     end
%     s(i) = s(i) - shift_grain;
%     t(i) = t(i) - shift_grain;
%     i = i + 1;
%     len = length(s);
% end

% figure(2)
% plot(grains([sp]))
% hold on
% plot(gbs,'linecolor','fuchsia')
% plot(gbs.triplePoints)
% hold off

for grain = 1:length(Grain_id)

    id = Grain_id(grain);
    a_temp = 0;
%     for edge = 1:length(s)
%         
%         if(s(edge) == grain || t(edge) == grain)
%             
%             a_temp = a_temp + a(edge);
%         end
%     end
%     factor(grain) = a_temp/(sum(grains(id).boundary.segLength)*scale);
    factor(grain) = 1;
    Grain_vol(grain) = grains(id).grainSize*scale^2*factor(grain);
    ph_i1(grain) = grains(id).meanOrientation.phi1;
    Ph_i(grain) = grains(id).meanOrientation.Phi;
    ph_i2(grain) = grains(id).meanOrientation.phi2;
end

%%%%%%%%%%%%%%%%%% find the shear coupling factor

low = pi/2 - 10*pi/180; high = pi/2 + 10*pi/180;
low2 = 3*pi/2 - 10*pi/180; high2 = 3*pi/2 + 10*pi/180;
plotIPDF(ebsd.orientations,xvector)

% psi = calcKernel(grains.meanOrientation);
% odf = calcDensity(ebsd.orientations,'kernel',psi);
% ori = discreteSample(odf,15000);
for i = 1:length(s)
    
    if(isnan(a_dir(i,1)))
        a_dir(i,:) = [ceil(axisangle(1)) ceil(axisangle(2))];
    end
    axisangle = mori(i).xyz;
    perp(i,:) = null(a_dir(i,:));
%     phia = acos(dot(a_dir(i,:),axisangle));
%     th = acos((trace(tt)-1)/2); % symmetric rotation angle
    phia = atan2(norm(cross([a_dir(i,:) 0],axisangle)),dot([a_dir(i,:) 0],axisangle));
%     th = atan2(norm(axisangle),(trace(tt)-1)/2);

    th = missorientation(i);
    th_gb = atan2(norm(cross([a_dir(i,:) 0],[1,0,0])),dot([1 0 0],[a_dir(i,:) 0])); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <<<<<<<--------------------------
    ang(i) = th_gb;
    if(phia < pi/18 && phia > -pi/18 )
        beta(i) = 0;
    else
        beta(i) = GetBeta(th);
    end
    if(beta(i) > 2)
        beta(i) = 0;
    end
    
    Rgb = [cos(th_gb) -sin(th_gb) 0;
    sin(th_gb) cos(th_gb) 0;
    0 0 1];

    Fpq_temp = eye(3); Fpq_temp(1,2) = beta(i);
    Fpq(:,:,i) = Fpq_temp;
    Fpq(:,:,i) = Rgb'*Fpq(:,:,i)*Rgb;
    test(i) = [perp(i,:) 0]*Fpq(:,:,i)*[a_dir(i,:) 0]';

    phi0(i) = GetPhi0(th);
end
phi1 = phi0*0.3;
% phi0 = rand(length(s),1) + 3.5;

fpqarea = zeros(3,3); fpqareatotal = zeros(3,3);
aa = a/sum(a);
for i = 1:length(cslarea)
    fpqarea = fpqarea + aa(cslarea(i))*Fpq(:,:,cslarea(i));
end
for i = 1:length(s)
    fpqareatotal = fpqareatotal + aa(i)*Fpq(:,:,i);
end
save GrainData_GBeng.mat Fpq ph_i1 Ph_i ph_i2 s t a Grain_vol phi1 phi0

% fprintf('STGBs for GB engineered Cu: %f percent \n',sum(isEdgeSTGB)/length(s)*100)



% other method for getting parameters for graph




% a = zeros(length(s),1); scale = 1e-6;
% a_dir = zeros(length(s),2);
% counter2 = 1;
% 
% for grain = 1:length(grain_id)
% 
%     id = grain_id(grain);
%     Grain_vol(grain) = grains(id).grainSize*scale;
%     phi1(grain) = grains(id).meanOrientation.phi1;
%     Phi(grain) = grains(id).meanOrientation.Phi;
%     phi2(grain) = grains(id).meanOrientation.phi2;
% end
% 
% 
% for edge = 1:length(s)
%     
%     q = grain_id(s(edge));
%     p = grain_id(t(edge));
% 
%     x = grains(q).triplePoints.x;
%     y = grains(q).triplePoints.y;
%     tp = [x y];
%     tp1 = sortrows(tp);
%     
%     x = grains(p).triplePoints.x;
%     y = grains(p).triplePoints.y;
%     tp = [x y];
%     tp2 = sortrows(tp);
%     
%     counter = 1;
%     trip = zeros(1,2);
%     [t1, uuu] = size(tp1);
%     [t2, uuu] = size(tp2);
% 
%     for i = 1:t1
%         for j = 1:t2
%         
%         px = tp1(i,:)-tp(j,:);
%         if(px == 0)
%             trip(counter,:) = tp1(i,:);
%             counter = counter + 1;
%             break;
%         end
%         end
%     end
% %     if(edge == 22)
% %         uuu = 2;
% %     end
%     [dubs, uuu] = size(trip);
%     
%     if(dubs > 2) % for 2 triple points on a grain
% 
% %         bonds(counter2) = edge;
% %         degenerate = [s(bonds) t(bonds)];
% 
%         if(dubs == 4) % only for 4 triple points
%             a(edge) = sqrt((trip(3,1)-trip(4,1))^2 + (trip(3,2)-trip(4,2))^2)*scale;
% 
%             a_dir(edge,:) = [(trip(4,1)-trip(3,1)), (trip(4,2)-trip(3,2))];
%             a_dir(edge,:) = a_dir(edge,:)./norm(a_dir(edge,:));
%         end
%     elseif(dubs < 2) % grains at the edge of the ebsd data -- assume straight line
% %         test(counter2) = edge;
% %         counter2 = counter2 + 1; % for debug
% 
%         if(edge == 100) % errors in triple points from MTEX
%             a(edge) = sqrt((trip(1,1)-5.8)^2 + (trip(1,2)-76)^2)*scale;
%     
%             a_dir(edge,:) = [(5.8-trip(1,1)), (trip(1,2)-76)];
%             a_dir(edge,:) = a_dir(edge,:)./norm(a_dir(edge,:));
%             continue
%         elseif(edge == 102)
%             a(edge) = sqrt((trip(1,1)-5.8)^2 + (trip(1,2)-76)^2)*scale;
%     
%             a_dir(edge,:) = [(trip(1,1)-5.8), (trip(1,2)-76)];
%             a_dir(edge,:) = a_dir(edge,:)./norm(a_dir(edge,:));
%             continue
%         end
%         [B,I] = sort(trip(1,:));
%         if(I(1) == 1)
%             a_dir(edge,:) = [1, 0];
%             a(edge) = abs(min(trip(1,:)))*scale;
%         else
%             a_dir(edge,:) = [0, 1];
%             a(edge) = abs(min(trip(1,:)))*scale;
%         end
%         
%     else % 2 triple points
%         a(edge) = sqrt((trip(1,1)-trip(2,1))^2 + (trip(1,2)-trip(2,2))^2)*scale;
%         
%         a_dir(edge,:) = [(trip(2,1)-trip(1,1)), (trip(2,2)-trip(1,2))];
%         a_dir(edge,:) = a_dir(edge,:)./norm(a_dir(edge,:));
%     end
% %         counter2 = counter2 + 1; % for debug
%     if(a(edge) == 0)
%         disp('error in graph')
%     end
% end
% 
% %%%%%%%%%%%%%%%%%% find the shear coupling factor
% Fpq = zeros(3,3,length(s));
% low = pi/2 - 5*pi/180; high = pi/2 + 5*pi/180;
% low2 = 3*pi/2 - 5*pi/180; high2 = 3*pi/2 + 5*pi/180;
% for i = 1:length(s)
%     
%         p = s(i);
%     q = t(i);
%     g1 = [phi1(p) Phi(p) phi2(p)];
%     g2 = [phi1(q) Phi(q) phi2(q)];
%     
%     r11 = cos(g1(1))*cos(g1(3)) - cos(g1(2))*sin(g1(1))*sin(g1(3));
%     r12 = -cos(g1(1))*sin(g1(3)) - cos(g1(2))*sin(g1(1))*cos(g1(3));
%     r13 = sin(g1(2))*sin(g1(1));
% 
%     r21 = cos(g1(3))*sin(g1(1)) + cos(g1(2))*cos(g1(1))*sin(g1(3));
%     r22 = cos(g1(2))*cos(g1(1))*cos(g1(3)) -  sin(g1(1))*sin(g1(3));
%     r23 = -sin(g1(2))*cos(g1(1));
% 
%     r31 = sin(g1(2))*sin(g1(3));
%     r32 = sin(g1(2))*cos(g1(3));
%     r33 = cos(g1(2));
%     
%     rot1 = [r11 r12 r13;
%     r21 r22 r23;
%     r31 r32 r33];
%     
%     r11 = cos(g2(1))*cos(g2(3)) - cos(g2(2))*sin(g2(1))*sin(g2(3));
%     r12 = -cos(g2(1))*sin(g2(3)) - cos(g2(2))*sin(g2(1))*cos(g2(3));
%     r13 = sin(g2(2))*sin(g2(1));
% 
%     r21 = cos(g2(3))*sin(g2(1)) + cos(g2(2))*cos(g2(1))*sin(g2(3));
%     r22 = cos(g2(2))*cos(g2(1))*cos(g2(3)) -  sin(g2(1))*sin(g2(3));
%     r23 = -sin(g2(2))*cos(g2(1));
% 
%     r31 = sin(g2(2))*sin(g2(3));
%     r32 = sin(g2(2))*cos(g2(3));
%     r33 = cos(g2(2));
% 
%     rot2 = [r11 r12 r13;
%         r21 r22 r23;
%         r31 r32 r33];
% 
%     tt = rot1*rot2';
% 
%     axisangle = zeros(2,1);
%     [G,H] = eig(tt);
%     
%     for j = 1:3
%         
%         if(isreal(G(:,j)))
%             axisangle(1) = G(1,j);
%             axisangle(1) = G(2,j);
%             break
%         end
%     end
% 
%     phi(i) = acos(dot(a_dir(i,:),axisangle));
%     th = acos((trace(tt)-1)/2); % symmetric rotation angle
%     
% %     if(i == 86) 
% %         uuu = 0;
% %     end
%     if(th > low && th < high) % 90 degrees
%         isEdgeSTGB(i) = true;
%         beta = GetBeta(th);
%     elseif(th > low2 && th < high2) % 270 degrees
%         isEdgeSTGB(i) = true;
%         beta = GetBeta(th);
%     else
%         isEdgeSTGB(i) = false;
%         beta = 0;
%     end
%     
%     Fpq_temp = eye(3); Fpq_temp(1,2) = beta;
%     Fpq(:,:,i) = Fpq_temp;
% end
% 
% sum(isEdgeSTGB)/length(s)

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

function [s t] = MakeGraph(Grain_id, neighbors)
    [iter, uuu] = size(neighbors);
    graph_n = cell(length(Grain_id),1);
    for vertex = 1:length(Grain_id)
        temp = [];
        counter = 1;
        for i = 1:iter
            if(neighbors(i,1) == Grain_id(vertex))
                next = neighbors(i,2);
                for j = 1:length(Grain_id)
                    if(next == Grain_id(j))
                        temp(counter) = neighbors(i,2);
                        counter = counter + 1;
                    end
                end
            end
        end
    graph_n{vertex} = temp;
    end
    
    for vertex = 1:length(Grain_id)
        [row,col] = size(graph_n{vertex});
        if(row == 0)
            for i = 1:length(Grain_id)
                if(i == vertex); continue; end
                temp = graph_n{i};
                [row,col] = size(temp);
                for j = 1:length(temp)
                    if(temp(j) == Grain_id(vertex) && col ~= 1)
                        graph_n{vertex} = Grain_id(i);
                        temp(j) = [];
                        graph_n{i} = temp; 
                        break
                    end
                end
                [t1,t2] = size(graph_n{vertex});
                if(t1 > 0)
                    break
                end
            end
        end
    end
    counter = 1;
    s = []; t = [];
    for i = 1:length(Grain_id)
        connecting = graph_n{i};
        for k = 1:length(connecting)
            for j = 1:length(Grain_id)
                if(connecting(k) == Grain_id(j))
                    t(end+1) = j;
                    s(end+1) = counter;
                end
            end
        end
        counter = counter + 1;
    end
end
% function phi0 = GetPhi0(th)
%     
%     angle = mod(th,pi/2);
%     if(angle < 0.2443)
%         phi0 = 5.279396826562586*angle + 0.01;
%     elseif(angle > 1.4835)
%         phi0 = -1.031324031235482*angle + 1.63;
%     else
%         tilt_angle = [14 19 22 25 28 31 34 38 50.5 53 57 59.5 63 64 68 71 74 77 79.5 81 83 85]*pi/180;
%         disp_energy =[1.3 1.6 1.25 .98 .8 .7 .65 .85 .4 2.8 .3 .35 .45 .48 .55 .48 .5 .4 .25 .23 .15 .1]; % GPa; 
%         phi0 = interp1(tilt_angle,disp_energy,angle);
%     end
% end
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
