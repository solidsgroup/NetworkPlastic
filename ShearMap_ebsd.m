clear all
fname = fullfile(mtexDataPath,'EBSD','Cu_GBE_EBSD_data','8july_2_0-Subset.txt');
CS = { loadCIF('Cu-Copper.cif')};
ebsd = loadEBSD_generic(fname,'CS',CS, 'ColumnNames', ...
{'Index' 'Phase' 'x' 'y' 'Euler1' 'Euler2' 'Euler3' 'MAD' 'OI' 'BC' 'BS'...
'Status'}, 'Bunge');

threshold = 5 *degree;
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',threshold);
ebsd(grains(grains.grainSize < 100)) = [];
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',threshold);
grains = smooth(grains,20);


% get colorbar scale
shr = 1:.002:2;
BlueRed = zeros(length(shr),3);
for i = 1:length(shr)
    BlueRed(i,:) = GetColor(shr(i),1 ,2);
end

% p = plot(nan,nan);
% x0=10;
% y0=10;
% width=550*2;
% height=400*2;
% set(gcf,'position',[x0,y0,width,height])
% axis normal

newmap = BlueRed;                %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(1/2 * ncol);    %1/2 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
colorbar('Ticks',[]);
endt = 1000;
% V0 = V(1,:)';
% V01 = repmat(V0,1,5001)';
% Dv = V-V01;

for t = 1:100:endt
    
    newMtexFigure('figSize','large')
    val = rand(1,length(grains))*2 - 1;
    mi = min(val);
    ma = max(val);
    for i = 1:length(grains)
        c = GetColor(val(i),mi ,ma);

        hold on
        plot(grains(i),'FaceColor',c,'micronbar','off')
       
    end
        hold on
        plot(grains.boundary)
        hold off
        mtexColorbar
        mtexColorMap blue2red
        setColorRange([mi ma])
        exportgraphics(gca,"Shear_ebsd.gif","Append",true)
        fprintf('t = %f\n',t)
end

function [c] = GetColor(shear,mi ,ma)
    if(shear == 0)
        shear_norm = 0;
        c = [.82 .82 .82];
        return
    else
        shear_norm = (shear-mi)/(ma-mi);
    end
    if(shear_norm < 0.5)
        percent = shear_norm*2;
        c = [percent percent 1];
    else
        percent = -2*shear_norm + 2;
        c = [1 percent percent];
    end
end
