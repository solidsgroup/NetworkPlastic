% plotting convention
clear all
plotx2east
plotzIntoPlane

% import the data
path = 'C:\Users\icrma\OneDrive\Documents\MATLAB\mtex-5.8.1\data\EBSD\Cu_ref_EBSD_data\';
ebsd = EBSD.load([path 'f_0_1-hx.ang'],'convertSpatial2EulerReferenceFrame','setting 2');

grains = ebsd.calcGrains;
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

%%

% Grain_id = [301 226 249 242 245 241 235 240 232 200 203 160 149 154 165 173 180 193 193 183 197 195 208 189 192 213 216 206 199 205 121 113 122 95 69 56 152 141 101 111 82 57 60 72 102 170 134 132 144 119 100 71 123 93 68 62 47 172 128 105 94 114 116 115 120 136 91 86 85 75 66];
Grain_id = [465 348 381 370 374 367 360 366 356 311 314 322 258 236 249 266 275 284 301 290 305 303 320 296 299 325 330 317 309 316 192 184 193 160 112 88 244 228 168 182 136 92 99 115 170 271 215 208 231 190 167 114 196 155 111 102 75 274 204 174 157 185 187 186 191 220 152 144 143 119 109];
figure(1)
plot(grains(Grain_id))

neighbors = grains.neighbors;

[s,t] = MakeGraph(Grain_id, neighbors);

a = zeros(length(s),1); scale = 1e-6;
a_dir = zeros(length(s),2);
for edge = 1:length(s)
    
    sp = [Grain_id(t(edge)) Grain_id(s(edge)) ];
    gbs =  grains.boundary(all(grains.boundary.grainId == [sp],2));
    
    [empty, uuu] = size(gbs);
    if(empty == 0)
        sp = [Grain_id(s(edge)) Grain_id(t(edge)) ];
        gbs =  grains.boundary(all(grains.boundary.grainId == [sp],2));
    end

    a(edge) = sum(gbs.segLength)*scale;
    pts = gbs.triplePoints.V;
    a_dir(edge,:) = [(pts(1,1)-pts(2,1)), (pts(1,2)-pts(2,2))];
    a_dir(edge,:) = a_dir(edge,:)./norm(a_dir(edge,:));
%     angle = mean(gbs.misorientation);
%     phi1(edge) = angle.phi1;
%     Phi(edge) = angle.Phi;
%     phi2(edge) = angle.phi2;
end
% sp = [Grain_id(t(1)) Grain_id(s(1)) ]
% gbs =  grains.boundary(all(grains.boundary.grainId == [sp],2));

% figure(2)
% plot(grains([sp]))
% hold on
% plot(gbs,'linecolor','fuchsia')
% plot(gbs.triplePoints)
% hold off
for grain = 1:length(Grain_id)

    id = Grain_id(grain);
    Grain_vol(grain) = grains(id).grainSize*scale;
    phi1(grain) = grains(id).meanOrientation.phi1;
    Phi(grain) = grains(id).meanOrientation.Phi;
    phi2(grain) = grains(id).meanOrientation.phi2;
end

%%%%%%%%%%%%%%%%%% find the shear coupling factor

low = pi/2 - 10*pi/180; high = pi/2 + 10*pi/180;
low2 = 3*pi/2 - 10*pi/180; high2 = 3*pi/2 + 10*pi/180;
for i = 1:length(s)
    
    p = s(i);
    q = t(i);
    g1 = [phi1(p) Phi(p) phi2(p)];
    g2 = [phi1(q) Phi(q) phi2(q)];
    
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

    tt = rot1*rot2';

    axisangle = zeros(2,1);
    [G,H] = eig(tt);
    
    for j = 1:3
        
        if(isreal(G(:,j)))
            axisangle(1) = G(1,j);
            axisangle(1) = G(2,j);
            break
        end
    end
    phi(i) = acos(dot(a_dir(i,:),axisangle));
    th = acos((trace(tt)-1)/2); % symmetric rotation angle
    

    if(phi(i) > low && phi(i) < high) % 90 degrees
        isEdgeSTGB(i) = true;
        beta = GetBeta(th);
    elseif(phi(i) > low2 && phi(i) < high2) % 270 degrees
        isEdgeSTGB(i) = true;
        beta = GetBeta(phi(i));
    else
        isEdgeSTGB(i) = false;
        beta = 0;
    end
    
    Fpq_temp = eye(3); Fpq_temp(1,2) = beta;
    Fpq(:,:,i) = Fpq_temp;

end

fprintf('STGBs for GB non-engineered Cu: %f percent \n',sum(isEdgeSTGB)/length(s)*100)

function beta = GetBeta(th)


  if(th <= pi/6)
      mu = 1;
  elseif(th >= pi/6 && th < pi/3)
      mu = 3/4;
  else
      mu = 1/2;
  end
  beta = (2*mu*tan(th/2) - 2*(1-mu)*tan(pi/4-th/2));

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