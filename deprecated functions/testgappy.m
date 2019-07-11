J = [3 11 14 17];
N = 200;
A = rand(N);

P = zeros(N,size(J,2));

for ii = 1 : size(J,2)
    P(J(ii),ii) = 1;
end

A1 = P'*A;
A2 = A(J,:);
A1
norm_diff = norm(A1-A2)




