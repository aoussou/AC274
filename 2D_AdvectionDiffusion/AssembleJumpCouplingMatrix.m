function AssembleJumpCouplingMatrix
% Assemble JumpCouplingMatrix, which provides the jump part of the
% numerical fluxes

global ElemNeighbors
global p
global EToN
global VX
global DofMap
global n_vars
global JumpCouplingMatrix


% Assemble JumpCouplingMatrix

% allocate space for the sparse coupling matrix
% the matrix will have a three (p+1)x(p+1) blocks
% for each element
if p==0
    n_local_dofs = 1;
    nodes_per_face = 1;
else
    n_local_dofs = 3*p;
    nodes_per_face = p+1;
end

n_faces = 3;
calN = n_local_dofs*size(EToN,1);
n_total_dofs = n_vars*calN;

[x_quad, w] = GaussQuad(nodes_per_face);
 w=w';
ref_nodes = RefNodeLocations;

JumpCouplingMatrix = spalloc(n_total_dofs, n_total_dofs, ...
    n_vars*3*size(EToN,1)*n_local_dofs*n_local_dofs);
for elem_id=1:size(EToN,1)
    
    % loop over the faces of this element
    for face_id=1:n_faces
        % get the element ID of the neighbor on the current side
        neighbor_id = ElemNeighbors(elem_id,face_id);
        
        if neighbor_id ~= -1
            if elem_id < neighbor_id
                continue;
            end
            [v, ~] = GetFaceLocalDofs(face_id);

            neighbor_face_id = find(ElemNeighbors(neighbor_id,:)==elem_id);
            [vn, ~] = GetFaceLocalDofs(neighbor_face_id);
            
            if p == 0 
            phi_e = GetPhi2D([1/3, 1/3]);%phi in columns, points in rows
            phi_n = GetPhi2D([1/3, 1/3]);%phi in columns, points in rows
            else
            G = segmap(x_quad,ref_nodes(v,:));
            phi_e = GetPhi2D(G);%phi in columns, points in rows
            Gn = segmap(x_quad,ref_nodes(vn,:));
            phi_n = GetPhi2D(Gn);%phi in columns, points in rows
            end
            
            w1 = VX(EToN(elem_id,v(1)),:);
            w2 = VX(EToN(elem_id,v(2)),:);
            LineElemJacobian = 0.5*sqrt((w2-w1)*(w2-w1)');
            % calculate the unit normal
            T = (w2-w1)/norm(w2-w1,2);
            normal_e = [T(2), -T(1)];
            normal_n = -[T(2), -T(1)];
            
            % build the dense matrix for the element-element coupling
            JumpCouplingMatrix_ee = zeros(n_local_dofs,n_local_dofs);
            JumpCouplingMatrix_ne = zeros(n_local_dofs,n_local_dofs);
            JumpCouplingMatrix_en = zeros(n_local_dofs,n_local_dofs);
            JumpCouplingMatrix_nn = zeros(n_local_dofs,n_local_dofs);


            for i=1:n_local_dofs
                for j=1:n_local_dofs
                    
                    JumpCouplingMatrix_ee(i,j) = JumpCouplingMatrix_ee(i,j)+...
                        normal_e * normal_e' * (w'*(phi_e(:,i) .* phi_e(:,j)))*LineElemJacobian;

                    JumpCouplingMatrix_en(i,j) = JumpCouplingMatrix_en(i,j)+...
                        normal_e * normal_n' * ...
                        (w'*(phi_e(:,i) .* flipud(phi_n(:,j))))*LineElemJacobian;

                    JumpCouplingMatrix_nn(i,j) = JumpCouplingMatrix_nn(i,j) +...
                        normal_n * normal_n' * ...
                        (w'*(flipud(phi_n(:,i)) .* flipud(phi_n(:,j))))*LineElemJacobian;

                    JumpCouplingMatrix_ne(i,j) = JumpCouplingMatrix_ne(i,j)+...
                        normal_n * normal_e' * ...
                        (w'*(flipud(phi_n(:,i)) .* phi_e(:,j)))*LineElemJacobian;
                end
            end
            
            % add the dense matrices to the sparse matrix!
            for i = 1:n_vars
                dofs_qi_e = DofMap{i}(elem_id,:);
                dofs_qi_n = DofMap{i}(neighbor_id,:);
                
                JumpCouplingMatrix(dofs_qi_e,dofs_qi_e) = JumpCouplingMatrix(dofs_qi_e,dofs_qi_e) + ...
                    JumpCouplingMatrix_ee;
                JumpCouplingMatrix(dofs_qi_e,dofs_qi_n) = JumpCouplingMatrix(dofs_qi_e,dofs_qi_n) + ...
                    JumpCouplingMatrix_en;
                JumpCouplingMatrix(dofs_qi_n,dofs_qi_n) = JumpCouplingMatrix(dofs_qi_n,dofs_qi_n) + ...
                    JumpCouplingMatrix_nn;
                JumpCouplingMatrix(dofs_qi_n,dofs_qi_e) = JumpCouplingMatrix(dofs_qi_n,dofs_qi_e) + ...
                    JumpCouplingMatrix_ne;
            end
            
        end

        
    end
end
