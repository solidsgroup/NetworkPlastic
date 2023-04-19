clear all
filename = 'C:\Users\icrma\source\repos\test\test\NP_stat\poly_microstructure_0';
opts = detectImportOptions(filename);
opts.Delimiter = ' ';
opts.DataLines = [1 inf];
C = readcell(filename,opts);
C = cellfun(@rmmissing, C, 'UniformOutput', false);

filename = 'C:\Users\icrma\source\repos\test\test\NP_stat\microstructure_0';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);
V = A.data(1,1); % number of grains
bimodal = A.data(1,2);
numSec = A.data(1,3);

[row,col] = size(C);
vert = cell(length(C),2);

for i = 1:row
    if(i == 100)
        aaa = 1;
    end
    tempx = [];
    tempy = [];
    for j = 1:col
        [r,c] = size(C{i,j});
        if(r == 0); break; end

        if(mod(j,2) == 1)
            tempx(end+1) = C{i,j};
        else
            tempy(end+1) = C{i,j};
        end
    end
    vert{i,1} = tempx;
    vert{i,2} = tempy;
end

for i = 1:row

    xlen = length(vert{i,1});
    ylen = length(vert{i,2});
    if(ylen ~= xlen)
        error('x and y have different number of vertices')
    end
end

fname = 'poly_';
if(bimodal)
    fname = append(fname,append('N',num2str(numSec)));
    fname = append(fname,append('_g',num2str(V)));
else
    fname = append(fname,append('_',num2str(V)));
end
fname = append(fname,'.mat');
save(fname, 'vert')  