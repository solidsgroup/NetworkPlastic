clear all
% x = [0.0734;
%     0.8714;
%     0.6947;
%     0.8414;
%     0.6762;
%     0.7014;
%     0.8954;
%     0.5874;
%     0.3252;
%     0.9949];
% y = [0.2075;
%     0.2863;
%     0.8331;
%     0.7702;
%     0.2006;
%     0.7648;
%     0.6159;
%     0.9643;
%     0.5037;
%     0.4787];
pnts = [];
numSec = 5; numpnts_fine = 100; numpnts_corse = 10; iter = 0;
for sec = 1:numSec
    offset = 1/numSec*iter;
    if(mod(sec,2) == 1)
        section = [rand(numpnts_fine,1),rand(numpnts_fine,1)*1/numSec + offset];
    else
        section = [rand(numpnts_corse,1),rand(numpnts_corse,1)*1/numSec + offset];
    end
    pnts = [pnts;section];
    iter = iter + 1;
end
[VX,VY] = voronoi(pnts(:,1),pnts(:,2));
% h = plot(VX,VY,'-b',x,y,'.r');
plot(VX,VY,'-b'); xlim([0 1]); ylim([0 1])
dt = delaunayTriangulation([pnts(:,1),pnts(:,2)]);
[V,R] = voronoiDiagram(dt);
neigh = cell(length(pnts(:,1)),1);
a = []; th_gb = [];
neighbors = [];

for i = 1:length(pnts(:,1)) % vertex values
    counter = 0;
    for j = 1:length(R{i}) % delete any vertex at infinity
        j = j - counter;
        if(R{i}(j) == 1)
            R{i}(j) = [];
            counter = counter + 1;
        end
    end
    % get volume
    xy = V(R{i},:);
    counter2 = 0;
    for t = 1:length(xy(:,1))
        if(xy(t-counter2,1) < -0.2 || xy(t-counter2,1) > 1.2)
            xy(t-counter2,:) = [];
            counter2 = counter2 + 1;
        elseif(xy(t-counter2,2) < -0.2 || xy(t-counter2,2) > 1.2)
            xy(t-counter2,:) = [];
            counter2 = counter2 + 1;
        end
    end
    xy = sortrows(xy); sum1 = 0; sum2 = 0;
    [row,col] = size(xy);
    if(row > 2)
        for k = 1:length(xy(:,1))
            ll = mod(k+1,length(xy(:,1))+1);
            if(ll == 0); ll = 1; end
            sum1 = sum1 + xy(k,1)*xy(ll,2);
            sum2 = sum2 + xy(k,2)*xy(ll,1);
        end
    %     vol2(i) = polyarea(xy(:,1),xy(:,2));
        Grain_vol(i) = abs(sum1 - sum(sum2))/2;
    else
        Grain_vol(i) = 1/(numpnts_fine+numpnts_corse);
    end
    ph_i1(i) = randn(1);
    Ph_i(i) = randn(1);
    ph_i2(i) = randn(1);
end
counter = 1;
for grain_id = 1:length(pnts(:,1)) % edge values
    tripnts = R{grain_id};
    for neighbor = 1:length(pnts(:,1))
        if(grain_id == neighbor)
            continue
        end
        l1 = [];
        nei_tripnts = R{neighbor};
        for i = 1:length(tripnts)
            for j = 1:length(nei_tripnts)
                if(tripnts(i) == nei_tripnts(j))
                    l1(end+1) = nei_tripnts(j);
                end
            end
        end
        size = length(l1);
        if(size == 2)
            neigh{grain_id}(end+1) = neighbor;
            v1 = V(l1(1),:); v2 = V(l1(2),:);
            a(end+1) = sqrt((v1(1)-v2(1))^2 + (v1(2)-v2(2))^2);
            if(a(end) > 0.2) a(end) = a(end)/100; end
            vv = sortrows([v1; v2]);
            vv = [vv(1,1)-vv(2,1), vv(1,2)-vv(2,2)];
            a_dir(counter,:) = vv./norm(vv);
            th_gb(end+1) = acosd(dot(a_dir(counter,:),[0,1]));
            counter = counter + 1;
        end
    end
end
grain_id = 1:length(pnts(:,1));
[s,t,graph_n, del] = MakeGraph(grain_id,neigh);
counter = 0;
low = pi/2 - 10*pi/180; high = pi/2 + 10*pi/180;
low2 = 3*pi/2 - 10*pi/180; high2 = 3*pi/2 + 10*pi/180;
for e = 1:length(del)
    del_edge = del(e);
    a(del_edge-counter) = [];
    th_gb(del_edge-counter) = [];
    a_dir(del_edge-counter,:) = [];
    counter = counter + 1;
end
for edge = 1:length(s)
    p = s(edge);
    q = t(edge);
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
    phia = acos(dot(a_dir(edge,:),axisangle));
    th = acos((trace(tt)-1)/2); % symmetric rotation angle
    
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
    vv = max(Grain_vol(p),Grain_vol(q));
    if(vv > 1/(numSec*numpnts_corse*10))
        beta = 0;
        atgb(edge) = true;
    end
    
    Rgb = [cos(th_gb(edge)) -sin(th_gb(edge)) 0;
    sin(th_gb(edge)) cos(th_gb(edge)) 0;
    0 0 1];

    Fpq_temp = eye(3); Fpq_temp(1,2) = beta;
    Fpq(:,:,edge) = Fpq_temp;
    Fpq(:,:,edge) = Rgb'*Fpq(:,:,edge)*Rgb;
end

phi1 = rand(length(s),1)*0.08 + .001;
phi0 = rand(length(s),1) + 3.5;

save BimodalData_N2.mat Fpq a Grain_vol s t ph_i1 Ph_i ph_i2 phi1 phi0 pnts
% hold on
% triplot(dt,'-r');
% hold off
%%
% pnts = [];
% numSec = 3; numpnts_fine = 5000; numpnts_corse = 10; iter = 0;
% for sec = 1:numSec
%     offset = 1/numSec*iter;
%     if(mod(sec,2) == 1)
%         section = [rand(numpnts_fine,1),rand(numpnts_fine,1)*1/numSec + offset];
%     else
%         section = [rand(numpnts_corse,1),rand(numpnts_corse,1)*1/numSec + offset];
%     end
%     pnts = [pnts;section];
%     iter = iter + 1;
% end
% voronoi(pnts(:,1),pnts(:,2))
function [s t graph_n, del] = MakeGraph(Grain_id, neigh)
    graph_n = neigh; neigh_edge = neigh;
    counter = 1;

    for vertex = 1:length(neigh)
        temp = neigh{vertex};
        for i = 1:length(temp)
            neigh_edge{vertex}(i) = counter;
            counter = counter + 1;
        end
    end
    del = zeros(neigh_edge{end}(end),1); test = [];
    for vertex = 1:length(Grain_id)
        for i = 1:length(graph_n{vertex})
            connect = graph_n{vertex}(i);
            for next_v = 1:length(Grain_id)
                if(vertex == next_v)
                    continue
                end
                counter = 0;
                for j = 1:length(graph_n{next_v})
                    if(vertex == graph_n{next_v}(j-counter) && connect == next_v)
                        graph_n{next_v}(j-counter) = [];
                        test(end+1) = neigh_edge{next_v}(j-counter);
                        del(neigh_edge{next_v}(j-counter)) = true;
                        counter = counter + 1;
                    end

                end
            end
        end
    end

    for vertex = 1:length(Grain_id) %%%%%%%%%%%%%%%%%%%%%%%%%%%% <-------
        [row,col] = size(graph_n{vertex});
        if(col == 0)
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
                if(t2 > 0)
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
    del = test;
end
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