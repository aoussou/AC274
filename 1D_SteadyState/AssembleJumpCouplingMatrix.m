function JumpCouplingMatrix = AssembleJumpCouplingMatrix
% Assemble JumpCouplingMatrix, which provides the jump part of the
% numerical fluxes

Globals1D;

% retrieve the basis functions at the edges of the reference element
boundary_phi = GetPhi(p, [-1;1]);

% Assemble JumpCouplingMatrix

% allocate space for the sparse coupling matrix
% the matrix will have a three (p+1)x(p+1) blocks
% for each element
n_local_dofs = p+1;
n_total_dofs = K*n_local_dofs;
JumpCouplingMatrix = spalloc(n_total_dofs, n_total_dofs, ...
    K*3*n_local_dofs*n_local_dofs);
for elem_id=1:K
    
    % get the local degree of freedom indices
    dofs_e = DofMap(elem_id,:);
    
    % loop over the faces of this element
    for face_id=1:2
        % get the element ID of the neighbor on the current side
        neighbor_id = ElemNeighbors(elem_id,face_id);
        
        % which face of the neighbor element are we looking at?
        % also, get the normal vectors (either +1 or -1 in 1D)
        % on the local and neighbor elements
        if face_id == 1
            neighbor_face_id = 2;
            normal_e         = -1;
            normal_n         = 1;
        else
            neighbor_face_id = 1;
            normal_e         = 1;
            normal_n         = -1;
        end
        
        % get the shape function values on this element
        boundary_phi_e = boundary_phi(face_id,:);
        
        
        % build the dense matrix for the element-element coupling
        JumpCouplingMatrix_ee = zeros(n_local_dofs,n_local_dofs);
        for i=1:n_local_dofs
            for j=1:n_local_dofs
                JumpCouplingMatrix_ee(i,j) = ...
                    normal_e * normal_e * ...
                    boundary_phi_e(i) * boundary_phi_e(j);
            end
        end
        
        % make sure we're not on a boundary, this matrix has no role
        % in assembling boundary fluxes
        if neighbor_id ~= -1
            % also we only do assembly if elem_id > neighbor_id
            % to avoid assembling terms twice
            if elem_id < neighbor_id
                continue;
            end
            

            
            % which face of the neighbor element are we looking at?
            % also, get the normal vectors on the local and neighbor elements
            if face_id == 1
                neighbor_face_id = 2;
                normal_e         = -1;
                normal_n         = 1;
            else
                neighbor_face_id = 1;
                normal_e         = 1;
                normal_n         = -1;
            end
            
            % get the shape function values on the neighbor
            boundary_phi_n = boundary_phi(neighbor_face_id,:);
            
            % build the dense matrix for the element-neighbor coupling
            JumpCouplingMatrix_en = zeros(n_local_dofs,n_local_dofs);
            for i=1:n_local_dofs
                for j=1:n_local_dofs
                    JumpCouplingMatrix_en(i,j) = ...
                        normal_e * normal_n * ...
                        boundary_phi_e(i) * boundary_phi_n(j);
                end
            end
            
            % build the dense matrix for the neighbor-neighbor coupling
            JumpCouplingMatrix_nn = zeros(n_local_dofs,n_local_dofs);
            for i=1:n_local_dofs
                for j=1:n_local_dofs
                    JumpCouplingMatrix_nn(i,j) = ...
                        normal_n * normal_n * ...
                        boundary_phi_n(i) * boundary_phi_n(j);
                end
            end
            
            % build the dense matrix for the neighbor-element coupling
            JumpCouplingMatrix_ne = zeros(n_local_dofs,n_local_dofs);
            for i=1:n_local_dofs
                for j=1:n_local_dofs
                    JumpCouplingMatrix_ne(i,j) = ...
                        normal_n * normal_e * ...
                        boundary_phi_n(i) * boundary_phi_e(j);
                end
            end
            
            % add the dense matrices to the sparse matrix!
            dofs_n = DofMap(neighbor_id,:);
            JumpCouplingMatrix(dofs_e,dofs_e) = JumpCouplingMatrix(dofs_e,dofs_e) + ...
                JumpCouplingMatrix_ee;
            JumpCouplingMatrix(dofs_e,dofs_n) = JumpCouplingMatrix(dofs_e,dofs_n) + ...
                JumpCouplingMatrix_en;
            JumpCouplingMatrix(dofs_n,dofs_n) = JumpCouplingMatrix(dofs_n,dofs_n) + ...
                JumpCouplingMatrix_nn;
            JumpCouplingMatrix(dofs_n,dofs_e) = JumpCouplingMatrix(dofs_n,dofs_e) + ...
                JumpCouplingMatrix_ne;
        end
        
    end
end

