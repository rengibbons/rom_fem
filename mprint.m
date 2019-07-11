function []=mprint(A,name,file)

fprintf(file,'\n--%s--  ', name);
n = size(A);
fprintf(file,'(%d x %d)\n', n(1), n(2));
fprintf(file,['\t' repmat('%5.4e\t',1,n(2)-1) '%5.4e\n'],A.');

end