function [evec,eval] = eigs3_fast(A)

if ~issymmetric(A)
    error('MATRIX FOR WHICH EVALS ARE QUERIED IS NOT SYMMETRIC')
elseif size(A,1)~=3 && size(A,1)~=3
    error('SIZE OF MATRIX FOR WHICH EVALS ARE QUERIED IS NOT 3X3')
end

tol = 1e-12;

% p1 = A(1,2)^2 + A(1,3)^2 + A(2,3)^2;
diag_test = abs(A(1,2)) + abs(A(1,3)) + abs(A(2,3));
eval = zeros(3);
evec = zeros(3);

if diag_test<tol && A(1,1)==A(2,2)==A(3,3)
    % A is diagonal isotropic
%     fprintf('In eig3_fast: diagona iso\n')
    eval(1,1) = A(1,1);
    eval(2,2) = A(2,2);
    eval(3,3) = A(3,3);
    
    evec = eye(3);
elseif diag_test<tol
    % A is diagonal and anisotropic
%     fprintf('In eig3_fast: diagonal aniso\n')
    min_A = abs(A(1,1));
    max_A = abs(A(1,1));
    min_ind = 1;
    max_ind = 1;
    for ii = 2 : 3
        if abs(A(ii,ii))>max_A
            max_A = abs(A(ii,ii));
            max_ind = ii;
        elseif abs(A(ii,ii))<min_A
            min_A = abs(A(ii,ii));
            min_ind = ii;
        end
    end
    mid_ind = setdiff(1:3,[min_ind max_ind]);
    
    eval(1,1) = A(min_ind,min_ind);
    eval(2,2) = A(mid_ind,mid_ind);
    eval(3,3) = A(max_ind,max_ind);
    
    evec(min_ind,1) = 1;
    evec(mid_ind,2) = 1;
    evec(max_ind,3) = 1;

else
    % A is symmetric but not diagonal
%     fprintf('In eig3_fast: not diagonal\n')
    q  = trace(A)/3;
    p1 = A(1,2)^2 + A(1,3)^2 + A(2,3)^2;
    p2 = (A(1,1) - q)^2 + (A(2,2) - q)^2 + (A(3,3) - q)^2 + 2 * p1;
    p  = sqrt(p2 / 6);
    B  = (1 / p) * (A - q * eye(3));
    r  = det(B) / 2;

    % In exact arithmetic for a symmetric matrix  -1 <= r <= 1
    % but computation error can leave it slightly outside this range.
    if r<=-1
        phi = pi / 3;
    elseif r>=1
        phi = 0;
    else
        phi = acos(r) / 3;
    end

    % the eigenvalues satisfy eval3 <= eval2 <= eval1
    eval(1,1) = q + 2 * p * cos(phi);
    eval(3,3) = q + 2 * p * cos(phi + (2*pi/3));
    eval(2,2) = 3 * q - eval(1,1) - eval(3,3);
   
%     % Reorder to match order of eigs for 2D plane strain
%     evaltemp(1,1) = eval(3,3);
%     evaltemp(2,2) = eval(2,2);
%     evaltemp(3,3) = eval(1,1);
%    
%     eval = evaltemp;
   
%     %% Basis full
%     ebas = zeros(3,3,3);
%     I = eye(3);
%     for kk = 1 : 3
%         B = rem(kk,3) + 1;
%         C = rem(B, 3) + 1;
%         for ii = 1 : 3
%             for jj = 1 : 3
%                 for ll = 1 : 3
%                     ebas(ii,jj,kk) = ebas(ii,jj,kk) ...
%                                    + (A(ii,ll) - eval(B,B)*I(ii,ll)) / (eval(B,B) - eval(kk,kk)) ...
%                                    * (A(ll,jj) - eval(C,C)*I(ll,jj)) / (eval(C,C) - eval(kk,kk));
%                 end
%             end
%         end
%     end

    %% Basis top row only
    ebascol = zeros(3,3);
    I = eye(3);
    for kk = 1 : 3
        B = rem(kk,3) + 1;
        C = rem(B, 3) + 1;
        for ii = 1 : 3
            for ll = 1 : 3
                ebascol(ii,kk) = ebascol(ii,kk) ...
                               + (A(ii,ll) - eval(B,B)*I(ii,ll)) / (eval(B,B) - eval(kk,kk)) ...
                               * (A(ll,1)  - eval(C,C)*I(ll,1))  / (eval(C,C) - eval(kk,kk));
            end
        end
    end
    
    [eval_one,~] = find(abs(eval-1)<tol);
    eval_others = setdiff(1:3,eval_one);
    
    if ~isempty(eval_one)
        evec(3,eval_one) = 1;
    end
    
    for ii = 1 : size(eval_others,2)
        ind = eval_others(ii);
        evec(:,ind) = ebascol(:,ind) / ebascol(1,ind);
        evec(:,ind) = evec(:,ind) / norm(evec(:,ind));
    end

end

end