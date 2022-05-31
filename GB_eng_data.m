clear all
fname = fullfile(mtexDataPath,'EBSD','Cu_GBE_EBSD_data','8july_2_0-Subset.txt');
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
% % plot(grains,grains.meanOrientation,'faceAlpha',0.8,'figSize','large')
% plot(ebsd,ebsd.orientations)
% % init override mode
% hold on
% 
% plot(GB,angle(oriGB)./degree,'linewidth',2)
% mtexColorMap blue2red
% mtexColorbar
% 
% hold off

%%
Grain_id = [60 82 81 91 74 46 20 16 50 43 49 52 66 67 47 39 30 63 18 35 32 24 31 34 44 64 61 75 87 118 104 95 115 119 112 103 71 84 89 92 99 116 114 117 166 174 147 187 185 190 146 232 155 173 170 171 127 153 180 159 145 165];
% 
% s = [1 1 1 2 2 2 2 2 3 4 4 5 6 6 7 7 7 8 9 10 10 11 11 11 12 12 13 13 13 14 14 15 15 15 15 15 15 15 16 16 17 18 18 18 18 19 20 20 21 21 22 23 23 24 24 25 26 27 29 29 29 30 31 32 33 33 33 33 34 34 34 35 35 36 36 36 37 37 38 38 39 39 39 39 40 41 41 41 42 42 42 43 43 43 44 45 45 45 45 46 46 47 48 48 49 49 49 50 50 51 52 52 52 52 52 52 52 52 53 53 53 54 54 55 55 56 56 57 58 58 58 59 59 60 60 60 60 61 61 62]';
% t = [2 3 7 7 10 5 13 15 2 3 2 38 9 2 8 15 15 15 10 15 12 10 9 12 13 37 9 14 11 12 37 16 17 18 18 19 20 37 17 18 18 37 25 26 27 20 21 22 23 18 23 24 18 25 18 26 29 28 31 27 33 29 30 31 28 18 34 32 32 31 52 33 36 18 33 52 10 39 13 41 14 40 36 52 36 39 4 52 4 41 43 41 44 45 52 44 42 46 52 42 49 46 47 46 48 45 50 51 45 45 35 40 33 43 51 57 56 55 51 54 50 50 55 53 56 59 58 58 52 52 60 60 58 61 62 52 52 52 58 52]';

figure(1)
plot(grains(Grain_id))

neighbors = grains.neighbors;

[s,t] = MakeGraph(Grain_id, neighbors);

a = zeros(length(s),1); scale = 1e-6;
a_dir = zeros(length(s),2);
shift = 0;
for edge = 1:length(s)
    edge = edge - shift;
    sp = [Grain_id(t(edge)) Grain_id(s(edge)) ];
    gbs =  grains.boundary(all(grains.boundary.grainId == [sp],2));
    
    [row, uuu] = size(gbs);
    if(row == 0)
        sp = [Grain_id(s(edge)) Grain_id(t(edge)) ];
        gbs =  grains.boundary(all(grains.boundary.grainId == [sp],2));
    end
%     if(edge == 100)
%         uuu = 1;
%     end
    a(edge) = sum(gbs.segLength)*scale;
    pts = gbs.triplePoints.V;
    [row, uuu] = size(pts);
    if(row == 1)
        pts2 = [gbs.x,gbs.y];
        a_dir(edge,:) = [(pts2(1,1)-pts2(end,1)), (pts(1,2)-pts(end,2))];
        a_dir(edge,:) = a_dir(edge,:)./norm(a_dir(edge,:));
    elseif(row == 0) % grain is circular 
%         s(edge) = [];
%         t(edge) = [];
%         a(edge) = [];
%         a_dir(edge,:) = [];
%         shift = shift + 1;
%         continue
    else
        a_dir(edge,:) = [(pts(1,1)-pts(2,1)), (pts(1,2)-pts(2,2))];
        a_dir(edge,:) = a_dir(edge,:)./norm(a_dir(edge,:));
    end
%     angle = mean(gbs.misorientation);
%     phi1(edge) = angle.phi1;
%     Phi(edge) = angle.Phi;
%     phi2(edge) = angle.phi2;
end
% sp = [Grain_id(t(1)) Grain_id(s(1)) ]
% gbs =  grains.boundary(all(grains.boundary.grainId == [sp],2));
% 
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
    
    if(i == 53)
        uuu = 2;
    end
    if(phi(i) > low && phi(i) < high) % 90 degrees
        isEdgeSTGB(i) = true;
        beta = GetBeta(th);
    elseif(phi(i) > low2 && phi(i) < high2) % 270 degrees
        isEdgeSTGB(i) = true;
        beta = GetBeta(th);
    else
        isEdgeSTGB(i) = false;
        beta = 0;
    end
    
    Fpq_temp = eye(3); Fpq_temp(1,2) = beta;
    Fpq(:,:,i) = Fpq_temp;

end

fprintf('STGBs for GB engineered Cu: %f percent \n',sum(isEdgeSTGB)/length(s)*100)



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