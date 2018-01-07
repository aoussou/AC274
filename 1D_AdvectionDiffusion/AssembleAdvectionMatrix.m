function S = AssembleAdvectionMatrix
% Assemble the right-hand side vector

Globals1D;

% NumQuadPts is the number of points in the Gauss quadrature
% rule that we use. Need to be able to integrate a 2*p degree
% polynomial exactly
NumQuadPts = p+1;
[x_quad, w_quad] = GaussQuad(p+1);

% retrieve the basis functions at the quadrature points
% on the reference element
phi = GetPhi(p, x_quad);

% get the derivative of the shape functions at the quadrature
% points
dphi = GetDPhi(p, x_quad);

% allocate space for the sparse mass matrix
% the matrix will have a (p+1)x(p+1) block
% for each element
n_local_dofs = p+1;
n_total_dofs = K*n_local_dofs;
S = spalloc(n_total_dofs, n_total_dofs, ...
                          K*n_local_dofs*n_local_dofs);
for elem_id=1:K
    local_dofs = DofMap(elem_id,:);

    % assemble the local RHS vector
    LocalAdvectionMatrix = zeros(n_local_dofs,n_local_dofs);
    for i=1:n_local_dofs
        for j=1:n_local_dofs
            % perform quadrature
            LocalAdvectionMatrix(i,j) = w_quad * ...
                         (phi(:,j) .* dphi(:,i));
        end
    end
    
    % add the local mass matrix to the global mass matrix
    S(local_dofs,local_dofs) = LocalAdvectionMatrix;
end

