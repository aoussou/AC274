function AssembleBoundaryConditionMatrix
% Assemble BoundaryConditionMatrix, which allows us to impose the BCs

Globals1D;
n_vars=1;
% retrieve the basis functions at the edges of the reference element
boundary_phi = GetPhi(p, [-1;1]);

% first assemble AvgCouplingMatrix

% allocate space for the sparse boundary condition matrix
n_local_dofs = p+1;
n_total_dofs = n_vars*K*n_local_dofs;
BoundaryConditionMatrix = spalloc(n_total_dofs, n_total_dofs, ...
                            n_vars*3*2*n_local_dofs*n_local_dofs);
for elem_id=[1,K]
    
    % loop over the faces of this element
    for face_id=1:2
        % get the element ID of the neighbor on the current side
        neighbor_id = ElemNeighbors(elem_id,face_id);
        
        % only treat boundary terms here
        if neighbor_id == -1

            if face_id == 1
                normal_e = -1;
            else
                normal_e = 1;
            end

            % get the shape function values on this element
            boundary_phi_e = boundary_phi(face_id,:);

            % build the dense matrix for the element-element coupling
            BoundaryConditionMatrix_ee = zeros(n_local_dofs,n_local_dofs);
            for i=1:n_local_dofs
                for j=1:n_local_dofs
                    BoundaryConditionMatrix_ee(i,j) = ...
                        normal_e * 0.5 * boundary_phi_e(i) * boundary_phi_e(j);
                end
            end

            for i = 1:n_vars
%                 dofs_qi_e = DofMap{i}(elem_id,:);
                dofs_qi_e = DofMap(elem_id,:);

                % outflow BC on both boundaries
                BoundaryConditionMatrix(dofs_qi_e,dofs_qi_e) = ...
                    BoundaryConditionMatrix(dofs_qi_e,dofs_qi_e) + 2*BoundaryConditionMatrix_ee;
            end
        end
    end
end
BC = BoundaryConditionMatrix;
