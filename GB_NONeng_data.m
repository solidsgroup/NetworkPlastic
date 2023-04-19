% plotting convention
clear all
plotx2east
plotzIntoPlane

% import the data
path = 'C:\Users\icrma\OneDrive\Documents\MATLAB\mtex-5.8.1\data\EBSD\Cu_ref_EBSD_data\';
ebsd = EBSD.load([path 'f_0_1-hx.ang'],'convertSpatial2EulerReferenceFrame','setting 2'); % f_5_2-hx.ang does work
CS = crystalSymmetry.load('CU-Copper');

threshold = 5 *degree;
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',threshold);
ebsd(grains(grains.grainSize < 50)) = [];
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',threshold);
grains = smooth(grains,20);
% % critical misorientation for grain reconstruction
% threshold = 10 *degree;
% 
% % first pass at reconstructing grains
% [grains, ebsd.grainId] = calcGrains(ebsd('Copper'),'angle',threshold);
% 
% % remove ebsd data that correspond to up to 4 pixel grains
% ebsd(grains(grains.grainSize < 50)) = []; % redo grain reconstruction -
% interpolate non-indexed space [grains, ebsd.grainId] =
% calcGrains(ebsd('Copper'),'angle',threshold); grains = smooth(grains,20);
% % remove all boundary grains % grains(grains.isBoundary) = [];
% 
% % remove too small irregular grains
% grains(grains.grainSize < grains.boundarySize / 2) = [];
% 
% % % set up the ipf coloring
% % cKey = ipfColorKey(ebsd.CS.properGroup);
% % cKey.inversePoleFigureDirection = xvector;
% % color = cKey.orientation2color(ebsd.orientations);
% % GB = grains.boundary;
% % oriGB = grains.boundary.misorientation;
% % 
% % plot(ebsd,color,'faceAlpha',0.5,'figSize','large')
% % 
% % % init override mode
% % hold on
% % 
% % plot(GB,angle(oriGB)./degree,'linewidth',2)
% % mtexColorMap blue2red
% % mtexColorbar
% % 
% % hold off
% p1 = [-602 -95]; p2 = [0 -95]; % horizontal 
% p3 = [-95 0]; p4 = [-95,-754]; % vertical
plot(ebsd,ebsd.orientations,'faceAlpha',0.8,'figSize','large')
hold on
plot(grains.boundary,'lineWidth',2)
% horizontal
% plot([p1(1),p2(1)],[p1(2),p2(2)],'-.b','lineWidth',3) 
% plot([p1(1),p2(1)],[p1(2)-95,p2(2)-95],'-.b','lineWidth',3) 
% plot([p1(1),p2(1)],[p1(2)-95*2,p2(2)-95*2],'-.b','lineWidth',3) 
% plot([p1(1),p2(1)],[p1(2)-95*3,p2(2)-95*3],'-.b','lineWidth',3) 
% plot([p1(1),p2(1)],[p1(2)-95*4,p2(2)-95*4],'-.b','lineWidth',3) 
% plot([p1(1),p2(1)],[p1(2)-95*5,p2(2)-95*5],'-.b','lineWidth',3) 
% plot([p1(1),p2(1)],[p1(2)-95*6,p2(2)-95*6],'-.b','lineWidth',3) 
% 
% % vertical
% plot([p3(1),p4(1)],[p3(2),p4(2)],'-.b','lineWidth',3) 
% plot([p3(1)-95,p4(1)-95],[p3(2),p4(2)],'-.b','lineWidth',3) 
% plot([p3(1)-95*2,p4(1)-95*2],[p3(2),p4(2)],'-.b','lineWidth',3) 
% plot([p3(1)-95*3,p4(1)-95*3],[p3(2),p4(2)],'-.b','lineWidth',3) 
% plot([p3(1)-95*4,p4(1)-95*4],[p3(2),p4(2)],'-.b','lineWidth',3) 
% plot([p3(1)-95*5,p4(1)-95*5],[p3(2),p4(2)],'-.b','lineWidth',3) 
% hold off

% restrict to iron to iron phase transition
gB = grains.boundary;

% select CSL(3) grain boundaries
gB3 = gB(angle(gB.misorientation,CSL(3,ebsd.CS)) < 3*degree);

figure(4)
plot(ebsd,log(ebsd.prop.iq),'figSize','large')
mtexColorMap black2white
CLim(gcm,[.5,5])

% the option 'FaceAlpha',0.4 makes the plot a bit transluent
hold on
plot(grains,grains.meanOrientation,'FaceAlpha',0.4,'linewidth',3)
hold off
hold on
plot(gB3,'lineColor','gold','linewidth',3,'DisplayName','CSL 3')
hold off
sum(gB3.segLength)/sum(gB.segLength)
%%

% Grain_id = [301 226 249 242 245 241 235 240 232 200 203 160 149 154 165 173 180 193 193 183 197 195 208 189 192 213 216 206 199 205 121 113 122 95 69 56 152 141 101 111 82 57 60 72 102 170 134 132 144 119 100 71 123 93 68 62 47 172 128 105 94 114 116 115 120 136 91 86 85 75 66];
% Grain_id = [465 348 381 370 374 367 360 366 356 311 314 322 258 236 249 266 275 284 301 290 305 303 320 296 299 325 330 317 309 316 192 184 193 160 112 88 244 228 168 182 136 92 99 115 170 271 215 208 231 190 167 114 196 155 111 102 75 274 204 174 157 185 187 186 191 220 152 144 143 119 109];

% Grain_id = [14 29 39 109 119 185 187 156 151 144 101 60 67 143 152 220 138 127 118 110 78 61 34 47 54 58 59 64 96 131 5 87 126 203 50 79 46];

% Grain_id = 1:length(grains);
% Grain_id = [47 63 69 94]; %twin file f_0_1-hx.ang
Grain_id = [553 508 586 522]; % random

% figure(1)
% plot(grains(Grain_id))
% p1 = [0,-95]; p2 = [-95 -95];
% hold on
% plot(grains.boundary,'lineWidth',2)
% plot([p1(1),p2(1)],[p1(2),p2(2)],'-.b','lineWidth',2)
% plot([p2(1),-95],[p1(2),0],'-.b','lineWidth',2)
% hold off

neighbors = grains.neighbors;

[s,t] = MakeGraph(Grain_id, neighbors);

a = zeros(length(s),1); scale = 1e-6;
a_dir = zeros(length(s),2);Edge_delete = zeros(length(s),1);
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
%     for i = 1:length(gbs)
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
    if(empty == 0)
        Edge_delete(edge) = true;
        continue
    end

    a(edge) = sum(gbs_split.segLength)*scale;
    pts = gbs.triplePoints.V;
    [row, uuu] = size(pts);
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

%     mill(edge) = round(Miller(mori(end),CS));
end
% sp = [Grain_id(t(1)) Grain_id(s(1)) ]
% gbs =  grains.boundary(all(grains.boundary.grainId == [sp],2));

% i = 1; len = length(s);
% shift = 0;
% while i <= len 
% 
%     if(Edge_delete(i))
%         s(i) = [];
%         t(i) = [];
%         a(i) = [];
%         a_dir(i,:) = [];
%         ph1(i) = [];
%         Ph(i) = [];
%         ph2(i) = [];
%         
%         shift = shift + 1;
%         for j = 2:length(s)
%             if(s(j-1) - s(j) < -1 )
%                 Grain_id(s(j) + s(j-1)) = [];
%             end
%         end
%     end
%     s(i) = s(i) - shift;
%     t(i) = t(i) - shift;
%     i = i + 1;
%     len = length(s);
% end
% Grain_id(8) = [];
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
for i = 1:length(s)

    axisangle = mori(i).xyz;

%     phia = acos(dot(a_dir(i,:),axisangle));
%     th = acos((trace(tt)-1)/2); % symmetric rotation angle
%     th_gb = acos(dot(a_dir(i,:),[0 1])); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <<<<<<<--------------------------
    perp(i,:) = null(a_dir(i,:));
    phia = atan2(norm(cross([a_dir(i,:) 0],axisangle)),dot([a_dir(i,:) 0],axisangle));
    th = missorientation(i);
    th_gb = atan2(norm(cross([a_dir(i,:) 0],[1,0,0])),dot([1 0 0],[a_dir(i,:) 0]));
    thh(i) = th_gb;
    if(phia < pi/18 && phia > -pi/18 )
        beta = 0;
    else
        beta = GetBeta(th);
    end
    if(beta > 2)
        beta = 0;
    end
    bbb(i) = beta;
    Rgb = [cos(th_gb) -sin(th_gb) 0;
    sin(th_gb) cos(th_gb) 0;
    0 0 1];

    Fpq_temp = eye(3); Fpq_temp(1,2) = beta;
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
save GrainData_GB_NONeng_random.mat Fpq ph_i1 Ph_i ph_i2 s t a Grain_vol phi1 phi0

% fprintf('STGBs for GB non-engineered Cu: %f percent \n',sum(isEdgeSTGB)/length(s)*100)

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
                    if(temp(j) == Grain_id(vertex) && col > 1)
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
