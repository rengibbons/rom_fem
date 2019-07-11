function     [d,sigma,sig]=mate01(F,J,mat_para,debug,file01)
%MATE01 Summary of this function goes here
%   Detailed explanation goes here

g12 = mat_para(1);
lambda = mat_para(2);

b = F*F';

kirch = g12.*b+(-g12+lambda/2*J*J-lambda/2).*eye(3);

sigma(1) = kirch(1,1)/J;
sigma(2) = kirch(2,2)/J;
sigma(3) = kirch(3,3)/J;
sigma(4) = kirch(1,2)/J;
sigma(5) = kirch(2,3)/J;
sigma(6) = kirch(1,3)/J;

sig=kirch./J;

Id2 = eye(3);

Nota= [ 1 4 6
        4 2 5
        6 5 3 ];

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                II4(i,j,k,l)=(0.5)*(Id2(i,k)*Id2(j,l) + Id2(i,l)*Id2(j,k));
                CmatE(i,j,k,l)= lambda.*J*J*Id2(i,j)*Id2(k,l)/J + (2*g12 + lambda - lambda*J*J)*II4(i,j,k,l)/J;
                d(Nota(i,j),Nota(k,l))=CmatE(i,j,k,l);
            end
        end
    end
end

if(debug)
    mprint(sigma,'CAUCHY STRESS',file01);
    mprint(d,'CONSISTANT TANGENT',file01);
end


end

