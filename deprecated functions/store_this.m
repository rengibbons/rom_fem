function [cell_stored] =store_this(a,b,location,cell_stored)

cell_stored{location,1}=a;
cell_stored{location,2}=b;

end