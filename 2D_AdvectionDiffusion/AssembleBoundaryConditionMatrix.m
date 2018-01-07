function AssembleBoundaryConditionMatrix
% Assemble BoundaryConditionMatrix, which allows us to impose the BCs

global p
global EToN
global VX
global ElemNeighbors
global n_vars
global DofMap
global BC

BoundaryConditionMatrix = cell(2,1);


n_faces = 3;
if p==0
    n_local_dofs = 1;
    nodes_per_face = 1;
else
    n_local_dofs = 3*p;
    nodes_per_face = p+1;
end
calN = n_local_dofs*size(EToN,1);
n_total_dofs=calN*n_vars;
[x_quad ,w] = GaussQuad(nodes_per_face);
w=w';
ref_nodes = RefNodeLocations;


for dim = 1:2
    BoundaryConditionMatrix{dim} = spalloc(n_total_dofs, n_total_dofs, ...
        n_vars*3*2*n_local_dofs*n_local_dofs);
end

for elem_id=1:size(EToN,1)
    for face_id=1:n_faces
        neighbor_id = ElemNeighbors(elem_id,face_id);
        
        % only treat boundary terms here
        if neighbor_id == -1
            %             vec_of_ele = [vec_of_ele;elem_id];
            [v, ~] = GetFaceLocalDofs(face_id);
            if p ==0
                boundary_phi_e = GetPhi2D([1/3 1/3]);%phi in columns, points in rows
                
            else
                G = segmap(x_quad,ref_nodes(v,:));
                boundary_phi_e = GetPhi2D(G);%phi in columns, points in rows
            end
            w1 = VX(EToN(elem_id,v(1)),:);
            w2 = VX(EToN(elem_id,v(2)),:);
            
            % calculate the unit normal
            T = (w2-w1)/norm(w2-w1,2);
            normal_e = [T(2), -T(1)];
            
            LineElemJacobian = 0.5*sqrt((w2-w1)*(w2-w1)');

            % build the dense matrix for the element-element coupling
            BoundaryConditionMatrix_ee = zeros(n_local_dofs,n_local_dofs,2);
            for dim = 1:2
                for i=1:n_local_dofs
                    for j=1:n_local_dofs
                        BoundaryConditionMatrix_ee(i,j,dim) = ...
                            normal_e(dim) * 0.5 * (w'*(boundary_phi_e(:,i).* boundary_phi_e(:,j)))*LineElemJacobian;
                    end
                end
                
                for i = 1:n_vars
                    dofs_qi_e = DofMap{i}(elem_id,:);
                    BoundaryConditionMatrix{dim}(dofs_qi_e,dofs_qi_e) = ...
                        BoundaryConditionMatrix{dim}(dofs_qi_e,dofs_qi_e) + 2*BoundaryConditionMatrix_ee(:,:,dim);
                end
            end
        end
    end
end
BC = BoundaryConditionMatrix;