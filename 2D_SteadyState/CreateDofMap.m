function CreateDofMap
% create an array that defines the map from local to
% global degrees of freedom

Globals2D;
n_vars = 1;
% define DofMap to be a cell array, which will store a matrix for
% each variable
K = size(EToN,1);
DofMap = cell(n_vars);
next_dof = 1;


if p==0
    n_local_dofs = 1;
else
    n_local_dofs = 3*p;
end

for var=1:n_vars
    DofMap{var} = zeros(K,n_local_dofs);
    for i=1:K
        DofMap{var}(i,:) = next_dof:(next_dof+n_local_dofs-1);%DofMap{variable}(element,local dofs)
        next_dof = next_dof + n_local_dofs;
    end
end