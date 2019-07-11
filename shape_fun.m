function [ N, dNdxi ] = shape_fun(xi,sys_vars)
%GSHAPEF Summary of this function goes here

element = sys_vars.element_type;

% xi(3,1) = [xi eta rho]

if (strcmp(element,'BRICK8'))

    N(1) = (1-xi(1))*(1-xi(2))*(1-xi(3))/8;
    N(2) = (1+xi(1))*(1-xi(2))*(1-xi(3))/8;
    N(3) = (1+xi(1))*(1+xi(2))*(1-xi(3))/8;
    N(4) = (1-xi(1))*(1+xi(2))*(1-xi(3))/8;
    N(5) = (1-xi(1))*(1-xi(2))*(1+xi(3))/8;
    N(6) = (1+xi(1))*(1-xi(2))*(1+xi(3))/8;
    N(7) = (1+xi(1))*(1+xi(2))*(1+xi(3))/8;
    N(8) = (1-xi(1))*(1+xi(2))*(1+xi(3))/8;


    dNdxi(1,1) = -(1-xi(2))*(1-xi(3))/8;
    dNdxi(1,2) = -(1-xi(1))*(1-xi(3))/8;
    dNdxi(1,3) = -(1-xi(1))*(1-xi(2))/8;


    dNdxi(2,1) =  (1-xi(2))*(1-xi(3))/8;
    dNdxi(2,2) = -(1+xi(1))*(1-xi(3))/8;
    dNdxi(2,3) = -(1+xi(1))*(1-xi(2))/8;


    dNdxi(3,1) =  (1+xi(2))*(1-xi(3))/8;
    dNdxi(3,2) =  (1+xi(1))*(1-xi(3))/8;
    dNdxi(3,3) = -(1+xi(1))*(1+xi(2))/8;


    dNdxi(4,1) = -(1+xi(2))*(1-xi(3))/8;
    dNdxi(4,2) =  (1-xi(1))*(1-xi(3))/8;
    dNdxi(4,3) = -(1-xi(1))*(1+xi(2))/8;


    dNdxi(5,1) = -(1-xi(2))*(1+xi(3))/8;
    dNdxi(5,2) = -(1-xi(1))*(1+xi(3))/8;
    dNdxi(5,3) =  (1-xi(1))*(1-xi(2))/8;

    dNdxi(6,1) =  (1-xi(2))*(1+xi(3))/8;
    dNdxi(6,2) = -(1+xi(1))*(1+xi(3))/8;
    dNdxi(6,3) =  (1+xi(1))*(1-xi(2))/8;

    dNdxi(7,1) =  (1+xi(2))*(1+xi(3))/8;
    dNdxi(7,2) =  (1+xi(1))*(1+xi(3))/8;
    dNdxi(7,3) =  (1+xi(1))*(1+xi(2))/8;

    dNdxi(8,1) = -(1+xi(2))*(1+xi(3))/8;
    dNdxi(8,2) =  (1-xi(1))*(1+xi(3))/8;
    dNdxi(8,3) =  (1-xi(1))*(1+xi(2))/8;

elseif (strcmp(element,'QUAD4'))
    N(1) = (1-xi(1))*(1-xi(2))/4;

    dNdxi(1,1) = -(1-xi(2))/4;
    dNdxi(1,2) = -(1-xi(1))/4;

    N(2) = (1+xi(1))*(1-xi(2))/4;

    dNdxi(2,1) =  (1-xi(2))/4;
    dNdxi(2,2) = -(1+xi(1))/4;

    N(3) = (1+xi(1))*(1+xi(2))/4;

    dNdxi(3,1) = (1+xi(2))/4;
    dNdxi(3,2) = (1+xi(1))/4;

    N(4) = (1-xi(1))*(1+xi(2))/4;

    dNdxi(4,1) = -(1+xi(2))/4;
    dNdxi(4,2) =  (1-xi(1))/4;
else
    
    disp('EROOR REOROREOR EOROE ROOEROEROEOR OERR OERROR ')

end

