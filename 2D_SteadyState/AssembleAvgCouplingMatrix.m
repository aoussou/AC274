function AssembleAvgCouplingMatrix
% Assemble AvgCouplingMatrix, which provides the average part of the
% numerical fluxes 1/2 (q- + q+)

global ElemNeighbors
global p
global EToN
global VX
global n_vars
global DofMap
global Cavg
% allocate space for the sparse coupling matrix
% the matrix will have a three (n-local_dofs)x(n-local_dofs) blocks
% for each element

n_faces = 3;
if p==0
    n_local_dofs = 1;
    nodes_per_face = 1;
else
    n_local_dofs = 3*p;
    nodes_per_face = p+1;
end

calN = n_local_dofs*size(EToN,1);
n_total_dofs = n_vars*calN;

ref_nodes = RefNodeLocations;
[x_quad, w] = GaussQuad(nodes_per_face);
w=w';

AvgCouplingMatrix = cell(2,1);
for dim = 1:2
    AvgCouplingMatrix{dim} = spalloc(n_total_dofs, n_total_dofs, ...
        n_vars*3*size(EToN,1)*n_local_dofs*n_local_dofs);
end

for elem_id=1:size(EToN,1)
    for face_id=1:n_faces
        
        % get the element ID of the neighbor on the current side
        neighbor_id = ElemNeighbors(elem_id,face_id);
        
        if neighbor_id ~= -1 %Skip if boundary face
            
            if elem_id < neighbor_id %Avoid repeats by always looking from larger elem id to smaller over an element face
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
            
            % calculate the unit normal
            T = (w2-w1)/norm(w2-w1,2);
            normal_e = [T(2), -T(1)];
            normal_n = -[T(2), -T(1)];
            
            %calculate the line element jac
            LineElemJacobian = 0.5*sqrt((w2-w1)*(w2-w1)');
            
            
            % build the dense matrix for the element-element coupling
            AvgCouplingMatrix_ee = zeros(n_local_dofs,n_local_dofs,2);
            AvgCouplingMatrix_en = zeros(n_local_dofs,n_local_dofs,2);
            AvgCouplingMatrix_nn = zeros(n_local_dofs,n_local_dofs,2);
            AvgCouplingMatrix_ne = zeros(n_local_dofs,n_local_dofs,2);
            
            
            for dim = 1:2
                for i=1:n_local_dofs
                    for j=1:n_local_dofs
                        
                        AvgCouplingMatrix_ee(i,j,dim) = AvgCouplingMatrix_ee(i,j,dim)+...
                            normal_e(dim) * 0.5 * (w'*(phi_e(:,i).* (phi_e(:,j))))*LineElemJacobian;
                        
                        AvgCouplingMatrix_en(i,j,dim) = AvgCouplingMatrix_en(i,j,dim)+...
                            normal_e(dim) * 0.5 * (w'*(phi_e(:,i).* flipud(phi_n(:,j))))*LineElemJacobian;
                        
                        AvgCouplingMatrix_nn(i,j,dim) = AvgCouplingMatrix_nn(i,j,dim)+...
                            normal_n(dim) * 0.5 *  (w'*(flipud(phi_n(:,i)).* flipud(phi_n(:,j))))*LineElemJacobian;
                        
                        AvgCouplingMatrix_ne(i,j,dim) = AvgCouplingMatrix_ne(i,j,dim)+...
                            normal_n(dim) * 0.5 *  (w'*(flipud(phi_n(:,i)) .* (phi_e(:,j))))*LineElemJacobian;
                        
                    end
                end
                
                % add the dense matrices to the sparse matrix!
                for i = 1:n_vars
                    
                    dofs_qi_e = DofMap{i}(elem_id,:);
                    dofs_qi_n = DofMap{i}(neighbor_id,:);
                    
                    
                    AvgCouplingMatrix{dim}(dofs_qi_e,dofs_qi_e) = AvgCouplingMatrix{dim}(dofs_qi_e,dofs_qi_e) + ...
                        AvgCouplingMatrix_ee(:,:,dim);
                    AvgCouplingMatrix{dim}(dofs_qi_e,dofs_qi_n) = AvgCouplingMatrix{dim}(dofs_qi_e,dofs_qi_n) + ...
                        AvgCouplingMatrix_en(:,:,dim);
                    AvgCouplingMatrix{dim}(dofs_qi_n,dofs_qi_n) = AvgCouplingMatrix{dim}(dofs_qi_n,dofs_qi_n) + ...
                        AvgCouplingMatrix_nn(:,:,dim);
                    AvgCouplingMatrix{dim}(dofs_qi_n,dofs_qi_e) = AvgCouplingMatrix{dim}(dofs_qi_n,dofs_qi_e) + ...
                        AvgCouplingMatrix_ne(:,:,dim);
                    
                end
            end
        end
    end
end

Cavg = AvgCouplingMatrix;