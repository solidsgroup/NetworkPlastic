clear all
load shear_map_test.mat

% get colorbar scale
shr = 1:.002:2;
BlueRed = zeros(length(shr),3);
for i = 1:length(shr)
    if(i == 250)
        aaa = 5;
    end
    BlueRed(i,:) = GetColor(shr(i),1 ,2);
end

figure(1)
p = plot(nan,nan);

newmap = BlueRed;                %starting map
ncol = size(newmap,1);           %how big is it?
zpos = 1 + floor(1/2 * ncol);    %1/2 of way through
newmap(zpos,:) = [1 1 1];        %set that position to white
colormap(newmap);                %activate it
colorbar('Ticks',[]);
endt = 5000;
V0 = V(1,:)';
V01 = repmat(V0,1,5001)';
Dv = V-V01;
for t = 1:50:endt
    shear = log(1+Dv(t,:));
    mi = min(shear);
    ma = max(shear);
    annotation('textbox',[.87 .1 .5 .05],'String',num2str(mi),'EdgeColor','none','BackgroundColor','w')
    annotation('textbox',[.87 .92 .5 .05],'String',num2str(ma),'EdgeColor','none','BackgroundColor','w')
    annotation('textbox',[.87 .5 .5 .05],'String',num2str((ma+mi)/2),'EdgeColor','none','BackgroundColor','w')

    for i = 1:length(Connection_Vertices)
    
        vert = vertices(Connection_Vertices{i},:);
        if(length(vert(:,1)) <= 2) continue; end
        
        c = GetColor(shear(i),mi ,ma);
        p1 = polyshape(vert(:,1),vert(:,2));
        hold on
        p = plot(p1,'FaceAlpha',1,'FaceColor', c);
        xlim([0 1])
        ylim([0 1])
    end
    exportgraphics(gca,"Shear.gif","Append",true)
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