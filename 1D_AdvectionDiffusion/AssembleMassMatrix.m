function M = AssembleMassMatrix
% Assemble the mass matrix

Globals1D;

% NumQuadPts is the number of points in the Gauss quadrature
% rule that we use. Need to be able to integrate a 2*p degree
% polynomial exactly
NumQuadPts = p+1;
[x_quad, w_quad] = GaussQuad(p+1);

% retrieve the basis functions at the quadrature points
% on the reference element
phi = GetPhi(p, x_quad);

% allocate space for the sparse mass matrix
% the matrix will have a (p+1)x(p+1) block
% for each element
n_local_dofs = p+1;
n_total_dofs = K*n_local_dofs;
M = spalloc(n_total_dofs, n_total_dofs, ...
                     K*n_local_dofs*n_local_dofs);
for elem_id=1:K

    % get the Jacobian of the mapping from reference
    % to physical element
    ElemJacobian = GetElemJacobian(elem_id);
    
    % assemble the (dense) local mass matrix
    LocalMassMatrix = zeros(n_local_dofs,n_local_dofs);
    for i=1:n_local_dofs
        for j=1:n_local_dofs
            % perform quadrature
            LocalMassMatrix(i,j) = w_quad * (phi(:,i) .* phi(:,j)) * ...
                                    ElemJacobian;
        end
    end
    
    % add the local mass matrix to the global mass matrix
    local_dofs = DofMap(elem_id,:);
    M(local_dofs,local_dofs) = LocalMassMatrix;
end

