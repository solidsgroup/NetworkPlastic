% plotting convention
clear all
plotx2east
plotzIntoPlane

% import the data
path = 'C:\Users\icrman\Documents\MATLAB\mtex-5.8.0\data\EBSD\Cu_ref_EBSD_data\';
ebsd = EBSD.load([path 'f_0_1-hx.ang'],'convertSpatial2EulerReferenceFrame','setting 2');

% critical misorientation for grain reconstruction
threshold = 10 *degree;

% first pass at reconstructing grains
[grains, ebsd.grainId] = calcGrains(ebsd('Copper'),'angle',threshold);

% remove ebsd data that correspond to up to 4 pixel grains
ebsd(grains(grains.grainSize < 500)) = [];

% redo grain reconstruction - interpolate non-indexed space
[grains, ebsd.grainId] = calcGrains(ebsd('Copper'),'angle',threshold);
grains = smooth(grains,500);
% remove all boundary grains
% grains(grains.isBoundary) = [];

% remove too small irregular grains
grains(grains.grainSize < grains.boundarySize / 2) = [];
pairs = grains.neighbors;
% plot the result
plot(ebsd,ebsd.orientations)
%%
scale = 50*1e-6; % 50 micro meters
% incidence_matrix = zeros(length(grains),length(pairs));
% adjacency_matrix = zeros(length(grains),length(grains));
% for i = 1:length(pairs)
%     counter = 1;
%     for j = 1:length(pairs)
%         if(pairs(i,1) == pairs(j,1))
%             incidence_matrix(pairs(i,1),counter) = 1;
%             incidence_matrix(pairs(j,2),counter) = 1;
%             adjacency_matrix(pairs(i,1),pairs(j,2)) = 1;
%         end
%         
%         counter = counter + 1;
%     end
% end
x = grains(1).triplePoints.x;
y = grains(1).triplePoints.y;
tp = [x y];
tp = sortrows(tp);

a12 = sqrt((tp(1,1)--131)^2+(tp(1,2)-0)^2)*scale; a21 = a12;
th_gb12 =  atan(abs(tp(1,2)-0)/abs(tp(1,1)--131)); th_gb21 = th_gb12;
a13 = sqrt((tp(1,1)-tp(2,1))^2+(tp(1,2)-tp(2,2))^2)*scale; a31 = a13;
th_gb13 = atan(-abs(tp(1,2)-tp(2,2))/abs(tp(1,1)-tp(2,1))); th_gb31 = th_gb13;
a14 = sqrt((tp(2,1)-tp(3,1))^2+(tp(2,2)-tp(3,2))^2)*scale; a41 = a14;
th_gb14 = atan(abs(tp(2,2)-tp(3,2))/abs(tp(2,1)-tp(3,1))); th_gb41 = th_gb14;
a15 = sqrt((tp(3,1)--71)^2 + (tp(3,2)-0)^2)*scale; a51 = a15;
th_gb15 = atan(abs(tp(3,2)-0)/abs(tp(3,1)--71)); th_gb51 = th_gb15;

x = (grains(8).triplePoints.x);
y = (grains(8).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a23 = sqrt((tp(3,1)-tp(4,1))^2+(tp(3,2)-tp(4,2))^2)*scale; a32 = a23;
th_gb23 = atan(abs(tp(3,2)-tp(4,2))/abs(tp(3,1)-tp(4,1))); th_gb32 = th_gb23;

x = (grains(9).triplePoints.x);
y = (grains(9).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a34 = sqrt((tp(1,1)-tp(5,1))^2+(tp(1,2)-tp(5,2))^2)*scale; a43 = a34;
th_gb34 = atan(abs(tp(1,2)-tp(5,2))/abs(tp(1,1)-tp(5,1))); th_gb43 = th_gb34;

x = (grains(17).triplePoints.x);
y = (grains(17).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a45 = sqrt((tp(4,1)-tp(6,1))^2+(tp(4,2)-tp(6,2))^2)*scale; a54 = a45;
th_gb45 = atan(abs(tp(4,2)-tp(6,2))/abs(tp(4,1)-tp(6,1))); th_gb54 = th_gb45;
a411 = sqrt((tp(2,1)-tp(5,1))^2+(tp(2,2)-tp(5,2))^2)*scale; a114 = a411;
th_gb411 = atan(-abs(tp(2,2)-tp(5,2))/abs(tp(2,1)-tp(5,1))); th_gb114 = th_gb411;

x = sort(grains(3).triplePoints.x);
y = sort(grains(3).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a56 = sqrt((tp(1,1)--30)^2 + (tp(1,2) -0)^2)*scale; a65 = a56;
th_gb56 = atan(-abs(tp(1,2) -0)/abs(tp(1,1)--30)); th_gb65 = th_gb56;
a67 = sqrt((tp(1,1)-0)^2 + (tp(1,2)-6.2)^2)*scale; a76 = a67;
th_gb67 = atan(-abs(tp(1,2)-6.2)/abs(tp(1,1)-0)); th_gb76 = th_gb67;

x = sort(grains(13).triplePoints.x);
y = sort(grains(13).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a57 = sqrt((tp(1,1)-tp(2,1))^2+(tp(1,2) - tp(2,2))^2)*scale;
a78 = sqrt((tp(1,1)-tp(3,1))^2+(tp(1,2) - tp(3,2))^2)*scale;

x = sort(grains(18).triplePoints.x);
y = sort(grains(18).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a58 = sqrt((tp(2,1)-tp(4,1))^2+(tp(2,2) - tp(4,2))^2)*scale; a85 = a58;
a48 = sqrt((tp(1,1)-tp(2,1))^2+(tp(1,2) - tp(2,2))^2)*scale; a84 = a48;
a118 = sqrt((tp(1,1)-tp(3,1))^2+(tp(1,2) - tp(3,2))^2)*scale; a811 = a118;

x = sort(grains(35).triplePoints.x);
y = sort(grains(35).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a1110 = sqrt((tp(1,1)-tp(2,1))^2+(tp(1,2) - tp(2,2))^2)*scale; a1011 = a1110;
a1210 = sqrt((tp(1,1)-tp(3,1))^2+(tp(1,2) - tp(3,2))^2)*scale; a1012 = a1210;
a1310 = sqrt((tp(3,1)-tp(4,1))^2+(tp(3,2) - tp(4,2))^2)*scale; a1013 = a1310;
a1410 = sqrt((tp(4,1)-tp(5,1))^2+(tp(4,2) - tp(5,2))^2)*scale; a1014 = a1410;
a1510 = sqrt((tp(5,1)-tp(6,1))^2+(tp(5,2) - tp(6,2))^2)*scale; a1015 = a1510;
a1610 = sqrt((tp(6,1)-tp(7,1))^2+(tp(6,2) - tp(7,2))^2)*scale; a1016 = a1610;
a910 = sqrt((tp(7,1)-tp(9,1))^2+(tp(7,2) - tp(9,2))^2)*scale; a109 = a910;
a710 = sqrt((tp(9,1)-tp(8,1))^2+(tp(9,2) - tp(8,2))^2)*scale; a107 = a710;
a810 = sqrt((tp(8,1)-tp(2,1))^2+(tp(8,2) - tp(2,2))^2)*scale; a108 = a810;

x = (grains(31).triplePoints.x);
y = (grains(31).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a1711 = sqrt((tp(2,1)-tp(3,1))^2+(tp(2,2) - tp(3,2))^2)*scale; a1117 = a1711;
a2011 = sqrt((tp(2,1)-tp(1,1))^2+(tp(2,2) - tp(1,2))^2)*scale; a1120 = a2011;
a1211 = sqrt((tp(4,1)-tp(5,1))^2+(tp(4,2) - tp(5,2))^2)*scale; a1112 = a1211;

x = (grains(52).triplePoints.x);
y = (grains(52).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a1213 = sqrt((tp(3,1)-tp(5,1))^2+(tp(3,2) - tp(5,2))^2)*scale; a1312 = a1213;

x = (grains(57).triplePoints.x);
y = (grains(57).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a1314 = sqrt((tp(3,1)-tp(5,1))^2+(tp(3,2) - tp(5,2))^2)*scale; a1413 = a1314;

x = (grains(63).triplePoints.x);
y = (grains(63).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a1415 = sqrt((tp(3,1)-tp(5,1))^2+(tp(3,2) - tp(5,2))^2)*scale; a1514 = a1415;

x = (grains(56).triplePoints.x);
y = (grains(56).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a1516 = sqrt((tp(4,1)-tp(5,1))^2+(tp(4,2) - tp(5,2))^2)*scale; a1615 = a1516;

x = (grains(34).triplePoints.x);
y = (grains(34).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a916 = sqrt((tp(1,1) - 0)^2+(tp(1,2) - -128)^2)*scale; a169 = a916;
a716 = sqrt((tp(2,1) - 0)^2+(tp(2,2) - -71)^2)*scale; a167 = a716;

x = (grains(23).triplePoints.x);
y = (grains(23).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a1719 = sqrt((tp(1,1)-tp(2,1))^2+(tp(1,2) - tp(2,2))^2)*scale; a1917 = a1719;
a173 = sqrt((tp(4,1)-tp(3,1))^2+(tp(4,2) - tp(3,2))^2)*scale; a317 = a173;
a1718 = sqrt((tp(2,1)-tp(3,1))^2+(tp(2,2) - tp(3,2))^2)*scale; a1817 = a1718;
a174 = sqrt((tp(4,1)-tp(6,1))^2+(tp(4,2) - tp(6,2))^2)*scale; a417 = a174;
a1720 = sqrt((tp(1,1)-tp(5,1))^2+(tp(1,2) - tp(5,2))^2)*scale; a2017 = a1720;

x = (grains(10).triplePoints.x);
y = (grains(10).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a1918 = sqrt((tp(1,1)-tp(2,1))^2+(tp(1,2) - tp(2,2))^2)*scale; a1819 = a1918;
a318 = sqrt((tp(3,1)-tp(4,1))^2+(tp(3,2) - tp(4,2))^2)*scale; a183 = a318;
a218 = sqrt((tp(1,1)-tp(4,1))^2+(tp(1,2) - tp(4,2))^2)*scale; a182 = a218;

x = (grains(26).triplePoints.x);
y = (grains(26).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a192 = sqrt((tp(4,1)-tp(5,1))^2+(tp(4,2) - tp(5,2))^2)*scale; a219 = a192;
a1920 = sqrt((tp(6,1)-tp(7,1))^2+(tp(6,2) - tp(7,2))^2)*scale; a2019 = a1920;

x = (grains(13).triplePoints.x);
y = (grains(13).triplePoints.y);
tp = [x y];
tp = sortrows(tp);

a79 = sqrt((tp(4,1)-0)^2+(tp(4,2) - -70)^2)*scale; a97 = a79;
hold on
plot(grains.boundary)
plot(grains(1).boundary,'LineWidth',4,'linecolor','b') % 1
plot(grains(8).boundary,'LineWidth',4,'linecolor','b') % 2
plot(grains(9).boundary,'LineWidth',4,'linecolor','b') % 3
plot(grains(17).boundary,'LineWidth',4,'linecolor','b') % 4
plot(grains(6).boundary,'LineWidth',4,'linecolor','b') % 5
plot(grains(3).boundary,'LineWidth',4,'linecolor','b')  % 6
plot(grains(13).boundary,'LineWidth',4,'linecolor','b') % 7
plot(grains(18).boundary,'LineWidth',4,'linecolor','b') % 8
plot(grains(34).boundary,'LineWidth',4,'linecolor','b')  % 9
plot(grains(35).boundary,'LineWidth',4,'linecolor','b')  % 10
plot(grains(31).boundary,'LineWidth',4,'linecolor','b')  % 11
plot(grains(52).boundary,'LineWidth',4,'linecolor','b')  % 12
plot(grains(57).boundary,'LineWidth',4,'linecolor','b')  % 13
plot(grains(63).boundary,'LineWidth',4,'linecolor','b')  % 14
plot(grains(56).boundary,'LineWidth',4,'linecolor','b')  % 15
plot(grains(55).boundary,'LineWidth',4,'linecolor','b')  % 16
plot(grains(23).boundary,'LineWidth',4,'linecolor','b')  % 17
plot(grains(10).boundary,'LineWidth',4,'linecolor','b')  % 18
plot(grains(26).boundary,'LineWidth',4,'linecolor','b')  % 19
plot(grains(33).boundary,'LineWidth',4,'linecolor','b')  % 20
plot(grains(13).triplePoints)
hold off

%%
s = [1 1 1 1 2 2 3 3 3 4 4 4 5 5 5 6 7 7 8 8 9 10 10 10 10 11 11 12 13 13 14 14 15 16 17 17 17 17 18 19 19 20]';
t = [2 3 4 5 3 18 4 17 18 5 8 11 6 8 7 7 8 10 10 11 7 12 15 16 9 10 12 13 10 14 10 15 16 9 4 11 19 20 17 18 20 11]';
a = [a12 a13 a14 a15 a23 a218 a34 a317 a318 a45 a48 a411 a56 a58 a57 a67 a78 a710 a810 a811 a97 a1012 a1015 a1016 a109 a1110 a1112 a1213 a1310 a1314 a1410 a1415 a1516 a169 a174 a1711 a1719 a1720 a1817 a1918 a1920 a2011]';
D = digraph(s,t,a);
Aa = full(adjacency(D,a)); 
Grain_id = [1;8;9;17;6;3;13;18;34;35;31;52;57;63;56;55;23;10;26;33];
Grain_vol = zeros(length(Grain_id),1);
for i = 1:length(Grain_id)
    id = Grain_id(i);
    Grain_vol(i) = grains(id).grainSize*scale;
    phi1(i) = grains(id).meanOrientation.phi1;
    Phi(i) = grains(id).meanOrientation.Phi;
    phi2(i) = grains(id).meanOrientation.phi2;
end
Fpq = zeros(3,3,length(s));
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
    
    th = acos((trace(tt)-1)/2); % symmetric rotation angle
%     if(abs(abs(th)-pi) < 0.15)
%         beta = 2*tan(pi/4-th/2);
%     else
%         beta = 2*tan(th/2);
%     end
%     beta = 2*0.01*tan(th/2) - 2*(1-0.01)*tan((2*pi/5-th)/2);
    beta = 2*0.1*tan(th/2) - 2*(1-0.1)*tan((2*pi/5-th)/2);
    Fpq_temp = eye(3); Fpq_temp(1,2) = beta;
    Fpq(:,:,i) = tt*Fpq_temp*tt'
end
for i = 1:length(a)
   
    phi0(i) = rand(1)/2+.1;
    phi11(i) = (rand(1)+.5)/100;
end 
save GrainData2.mat Aa phi0 phi11 phi1 phi2 Phi Fpq s t Grain_vol a