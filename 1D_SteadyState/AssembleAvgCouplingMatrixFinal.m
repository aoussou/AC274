function Cavg = AssembleAvgCouplingMatrixFinal
% Assemble AvgCouplingMatrix, which provides the average part of the
% numerical fluxes

Globals1D;

% retrieve the basis functions at the edges of the reference element
boundary_phi = GetPhi(p, [-1;1]);

% first assemble AvgCouplingMatrix

% allocate space for the sparse coupling matrix
% the matrix will have a three (p+1)x(p+1) blocks
% for each element
n_local_dofs = p+1;
n_total_dofs = K*n_local_dofs;
AvgCouplingMatrix = spalloc(n_total_dofs, n_total_dofs, ...
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
        
%         if problem_2_flag
%             if elem_id==1
%                 if face_id ==1
%                     continue;
%                 end
%             end
%         end
        
        % build the dense matrix for the element-element coupling
        AvgCouplingMatrix_ee = zeros(n_local_dofs,n_local_dofs);
        for i=1:n_local_dofs
            for j=1:n_local_dofs
                AvgCouplingMatrix_ee(i,j) = ...
                    normal_e * 0.5 * boundary_phi_e(i) * boundary_phi_e(j);
            end
        end
        
        % assemble boundary condition terms
        if neighbor_id == -1
            
            if normal_e * advection_velocity > 0
                % outflow BC
                % impose a one-sided flux, f(q^-), from inside the domain
                % using AvgCouplingMatrix_ee
                % We need to double AvgCouplingMatrix_ee here since we're
                % not "averaging" here
                AvgCouplingMatrix(dofs_e,dofs_e) = ...
                    AvgCouplingMatrix(dofs_e,dofs_e) + 2*AvgCouplingMatrix_ee;
                %inflow BC
                
                
                
            else
                % inflow BC
                % impose one-sided flux from outside the domain, hence
                % do nothing here: outflow BCs need to be imposed on the
                % RHS vector once it is assembled
            end
            
        else
            
            % we only do assembly if elem_id > neighbor_id
            % to avoid assembling terms twice
            if elem_id < neighbor_id
                continue;
            end
            
            
            
            % get the shape function values on the neighbor
            boundary_phi_n = boundary_phi(neighbor_face_id,:);
            
            % build the dense matrix for the element-neighbor coupling
            AvgCouplingMatrix_en = zeros(n_local_dofs,n_local_dofs);
            for i=1:n_local_dofs
                for j=1:n_local_dofs
                    AvgCouplingMatrix_en(i,j) = ...
                        normal_e * 0.5 * boundary_phi_e(i) * boundary_phi_n(j);
                end
            end
            
            % build the dense matrix for the neighbor-neighbor coupling
            AvgCouplingMatrix_nn = zeros(n_local_dofs,n_local_dofs);
            for i=1:n_local_dofs
                for j=1:n_local_dofs
                    AvgCouplingMatrix_nn(i,j) = ...
                        normal_n * 0.5 * boundary_phi_n(i) * boundary_phi_n(j);
                end
            end
            
            % build the dense matrix for the neighbor-element coupling
            AvgCouplingMatrix_ne = zeros(n_local_dofs,n_local_dofs);
            for i=1:n_local_dofs
                for j=1:n_local_dofs
                    AvgCouplingMatrix_ne(i,j) = ...
                        normal_n * 0.5 * boundary_phi_n(i) * boundary_phi_e(j);
                end
            end
            
            % add the dense matrices to the sparse matrix!
            dofs_n = DofMap(neighbor_id,:);
            AvgCouplingMatrix(dofs_e,dofs_e) = AvgCouplingMatrix(dofs_e,dofs_e) + ...
                AvgCouplingMatrix_ee;
            AvgCouplingMatrix(dofs_e,dofs_n) = AvgCouplingMatrix(dofs_e,dofs_n) + ...
                AvgCouplingMatrix_en;
            AvgCouplingMatrix(dofs_n,dofs_n) = AvgCouplingMatrix(dofs_n,dofs_n) + ...
                AvgCouplingMatrix_nn;
            AvgCouplingMatrix(dofs_n,dofs_e) = AvgCouplingMatrix(dofs_n,dofs_e) + ...
                AvgCouplingMatrix_ne;
            
        end
    end
end

% if problem_2_flag
%     [inflow,outflow] = problem2_BC;
%     AvgCouplingMatrix(end-p-1:end,end-p-1:end) = outflow;
%     AvgCouplingMatrix(1:p+1,1:p+1) = zeros(p+1);
% end
    
    Cavg = AvgCouplingMatrix;
    
