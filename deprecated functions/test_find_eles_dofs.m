function [eles,dofs] = test_find_eles_dofs(J,connect_list,element_dof)

for ii = 1 : size(J,2)
    
    
    [eles,~] = find(element_dof==J(ii));
    eles     = sort(eles);

    dofs = element_dof(eles,:);
    dofs = sort(unique(dofs(:)));
    
    fprintf('j = %.0f\n',J(ii))
    fprintf('\t eles = {')
    for jj = 1 : size(eles,1)
        fprintf(' %.0f',eles(jj))
    end
    fprintf(' }\n\t dofs = {')
    for jj = 1 : size(dofs,1)
        fprintf(' %.0f',dofs(jj))
    end
    fprintf(' }\n')
    
end

eles=0;
dofs=0;
end