function AssembleAdvectionMatrix

%Should be "diagonal" with no coupling between elements and "stamps of
%smaller n_local_dofsxn_local_dofs matrices

global EToN
global p
global n_vars
global DofMap
global S

if p==0
    n_local_dofs = 1;
    N_quad_pts = 1;
else
    n_local_dofs = 3*p;
    N_quad_pts = 3;
end
n_total_dofs = n_local_dofs*size(EToN,1);
[x_quad, w] = GaussQuad2D(N_quad_pts);
phi = GetPhi2D(x_quad);
dphi = GetDPhi2D(x_quad);

% allocate space for the sparse mass matrix
% the matrix will have a (n_loc_dofs)x(n_loc_dofs) block
% for each element


AdvectionMatrix = cell(2,1);

AdvectionMatrix{1} = spalloc(n_total_dofs, n_total_dofs, ...
    n_vars* size(EToN,1)*n_local_dofs*n_local_dofs);
AdvectionMatrix{2} = spalloc(n_total_dofs, n_total_dofs, ...
    n_vars* size(EToN,1)*n_local_dofs*n_local_dofs);

for elem_id=1:size(EToN,1)
    
    ElemJacobian = GetElemJacobian2D(elem_id);              %OUTPUT: the element jacobian, a 2x2 matrix
    invJac_T = inv(ElemJacobian');                          %OUTPUT: the inverse of the transpose of the element jacobian, a 2x2 matrix
    
    LocalAdvectionMatrix = zeros(n_local_dofs,n_local_dofs,2);
    
    for dim = 1:2
        
        for i=1:n_local_dofs
            for j=1:n_local_dofs
%                 LocalAdvectionMatrix(i,j,dim) =  (invJac_T(dim,:)*dphi(:,:,i)')*(w.*phi(:,j))*det(ElemJacobian);
                LocalAdvectionMatrix(i,j,dim) =  (invJac_T(dim,:)*dphi(:,:,i)')*(w.*phi(:,j))*det(ElemJacobian);

            end
        end

        for i = 1:n_vars
            dofs_qi = DofMap{i}(elem_id,:);
            AdvectionMatrix{dim}(dofs_qi,dofs_qi) = LocalAdvectionMatrix(:,:,dim);
        end
        
    end
end

S = AdvectionMatrix;