function CreateDofMap
% create an array that defines the map from local to
% global degrees of freedom

Globals1D;

DofMap = zeros(K,p+1);
next_dof = 1;

for i=1:K
    DofMap(i,:) = next_dof:(next_dof+p);
    next_dof = next_dof + p + 1;
end

