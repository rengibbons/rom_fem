%%
function write_mesh()

format long
xDim   = 0.1;
yDim   = 0.1;
nXele  = 20;
nYele  = 20;
folder = 'input';
file   = 'square20x20';

write_coordinates(xDim,yDim,nXele,nYele,folder,file);
write_elements(nXele,nYele,folder,file);
fprintf(strcat('Written mesh files for:\t',folder,'/',file,'\n'))
end

%%
function write_coordinates(xDim,yDim,nXele,nYele,folder,file)
nCoord = (nXele+1)*(nYele+1);
coordinates = zeros(nCoord,4);

nXnodes = nXele + 1;
for ii = 1 : nCoord
   coordinates(ii,1) = ii;
   coordinates(ii,3) = (mod(ii-1,nXnodes)) / nXele * xDim;
   coordinates(ii,4) = floor((ii-1)/nXnodes) / nYele * yDim;
end
dlmwrite(strcat(folder,'/',file,'.nodes'),coordinates,'delimiter','\t','precision',10)
end

%%
function write_elements(nXele,nYele,folder,file)
nEle = nXele * nYele;
elements = zeros(nEle,7);
elements(:,3) = 1;

count = 1;
for y = 1 : nYele
    for x = 1 : nXele
        elements(count,1) = count;
        elements(count,4) = x   + (nXele+1)*(y-1);
        elements(count,5) = x+1 + (nXele+1)*(y-1);
        elements(count,6) = x+1 + (nXele+1)*y;
        elements(count,7) = x   + (nXele+1)*y;
        count = count + 1;
    end
end
dlmwrite(strcat(folder,'/',file,'.elems'),elements,'delimiter','\t','precision',10)
end
