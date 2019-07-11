clc, clear all, close all

% causes error
% A = [-1.1 0 0;
%      0 1.1 0;
%      0 0 1];

A = [-1.1 0 0;
     0 1 0;
     0 0 -1.1];

N = 1;
format long 

for ii = 1 : N
    A = zeros(3);
    A(1,1) = (-1)^round(rand) * rand;
    A(2,2) = (-1)^round(rand) * rand;
    A(3,3) = (-1)^round(rand) * rand;
    
    A(1,2) = (-1)^round(rand) * rand;
    A(1,3) = (-1)^round(rand) * rand;
    A(2,3) = (-1)^round(rand) * rand;
    
    A(2,1) = A(1,2);
    A(3,1) = A(1,3);
    A(3,2) = A(2,3);

    [vec1,val1] = eigs(A);
    [vec2,val2] = eigs3_fast(A);
    

    vecdiff = norm(vec1-vec2);
    valdiff = norm(val1-val2);
    
    fprintf('vecdiff = %.8f    valdiff = %.8f\n',[vecdiff valdiff])
    
%     for jj = 1 : 3
%         e_test1 = (A-val1(jj,jj)*eye(3)) * vec1(:,jj);
%         e_test2 = (A-val2(jj,jj)*eye(3)) * vec2(:,jj);
%         fprintf('test1 = %.12f    test2 = %.12f\n',[e_test1 sum(e_test2])
%     end
    
    if vecdiff~=0
        fprintf('BBBBBBBLLLLLLLLLLLLLLLLAAAAAAAAAAAAAAAAAHHHHHHHHHHHHH\n')
        val1
        val2
%         A
%         val1
%         val2
%         vec2
%         vec1
%         vec2
    end
end


 
% N = 10000;
% format long 
% 
% total_eigs_time      = 0;
% total_eig3_fast_time = 0;
% 
% for ii = 1 : N
%     A = zeros(3);
%     A(1,1) = rand;
%     A(2,2) = rand;
%     A(3,3) = 1;
%     A(1,2) = rand;
%     A(2,1) = A(1,2);
% %     A(1,3) = rand;
% %     A(2,3) = rand;
% %     A(3,1) = A(1,3);
% %     A(3,2) = A(2,3);
% %     A = [-5 0 0; 
% %          0 1 0;
% %          0 0 2];
%     tic
%     [evec1,eval1] = eigs(A);
%     t1 = toc;
%     
%     tic
%     [evec2,eval2] = eig3_fast(A);
%     t2 = toc;
% 
%     total_eigs_time      = total_eigs_time + t1;
%     total_eig3_fast_time = total_eig3_fast_time + t2;
%     
% %     tol = 1e-8;
% %     if abs(eval1(1,1)-eval2(1,1))>tol || abs(eval1(2,2)-eval2(2,2))>tol || abs(eval1(3,3)-eval2(3,3))>tol
% %         fprintf('evals not ordered same. skip to next test.\n')
% %         continue
% %     end
% %     
% %     for jj = 1 : 3
% %         if evec1(1,jj)+evec2(1,jj)<tol && evec1(1,jj)~=0
% %             evec2(:,jj) = -evec2(:,jj);
% %         end
% %     end
% %     
% %     normdiff = norm(evec1-evec2);
% %     fprintf('diff = %.10f\n',normdiff)
% 
% end
% 
% total_eigs_time
% total_eig3_fast_time
% times_fast = total_eigs_time / total_eig3_fast_time
