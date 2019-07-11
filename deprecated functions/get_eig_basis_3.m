function [A_bas,A_val] = get_eig_basis_3(A)

[A_vec,A_val] = eig(A);

A_bas = zeros(3,3,3);
for ii = 1 : 3
    A_bas(:,:,ii) = A_vec(:,ii) * A_vec(:,ii)';
end

end