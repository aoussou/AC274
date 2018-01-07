function jumpMatrix2
% Assemble AvgCouplingMatrix, which provides the average part of the
% numerical fluxes

Globals2D;

% retrieve the basis functions at the edges of the reference element


% first assemble AvgCouplingMatrix

% allocate space for the sparse coupling matrix
% the matrix will have a three (p+1)x(p+1) blocks
% for each element
switch n_dim
    case 1
        n_faces  =2;
        n_local_dofs = p+1;
        n_total_dofs = n_vars*K*n_local_dofs;
        boundary_phi = GetPhi(p, [-1;1]);
    case 2
        n_faces = 3;
        if p==0
            n_local_dofs = 1;
            calN = n_local_dofs*size(EToN,1);
            nodes_per_face = 1;
        else
            n_local_dofs = 3*p;
            nodes_per_face = p+1;
            calN = n_local_dofs*size(EToN,1);
        end
        n_total_dofs = n_vars*calN;
        ref_nodes = RefNodeLocations;
%         phi = GetPhi2D(ref_nodes);%phi in columns, points in rows
        
end

    JumpCoupMat= spalloc(n_total_dofs, n_total_dofs, ...
        n_vars*3*K*n_local_dofs*n_local_dofs);

for elem_id=1:K
    
    % loop over the faces of this element
    for face_id=1:n_faces

        % get the element ID of the neighbor on the current side
        neighbor_id = ElemNeighbors(elem_id,face_id);
        
        % which face of the neighbor element are we looking at?
        % also, get the normal vectors (either +1 or -1 in 1D)
        % on the local and neighbor elements
        
        if n_dim == 2
            [v, V_e] = GetFaceLocalDofs(face_id);
            [x_quad, w] = GaussQuad(3);
            w=w';
            ref_nodes = RefNodeLocations;
            
%             phys_pts = GetPhysicalPoints2D(ref_nodes(v));
            G = segmap(x_quad,ref_nodes(v,:));
            phi = GetPhi2D(G);%phi in columns, points in rows
            
            
        elseif n_dim ==1
            if face_id == 1
                neighbor_face_id = 2;
                normal_e         = -1;
                normal_n         = 1;
            else
                neighbor_face_id = 1;
                normal_e         = 1;
                normal_n         = -1;
            end
        end
        if neighbor_id ~= -1
            
            neighbor_id = ElemNeighbors(elem_id,face_id);
            neighbor_face_id = find(ElemNeighbors(neighbor_id,:)==elem_id);
            [vn, V_n] = GetFaceLocalDofs(neighbor_face_id);
        

        end
        
        
        w1 = VX(EToN(elem_id,v(1)),:);
        w2 = VX(EToN(elem_id,v(2)),:);
        
        % calculate the unit normal
        T = (w2-w1)/norm(w2-w1,2);
        normal_e = [T(2), -T(1)];
        normal_n = -[T(2), -T(1)];
        
        %calculate the line element jac
        LineElemJacobian = 0.5*sqrt((w2-w1)*(w2-w1)');
        
        % get the shape function values on this element
        if n_dim ==1
%             boundary_phi_n = boundary_phi(neighbor_face_id,:);
        else
            
            boundary_phi_e = phi;
        end
        %         boundary_phi_e = boundary_phi(face_id,:);
        
        % build the dense matrix for the element-element coupling
        AvgCouplingMatrix_ee = zeros(n_local_dofs,n_local_dofs);
        
        
%             for node = 1:nodes_per_face
                for i=1:n_local_dofs
                    for j=1:n_local_dofs
                        
                        AvgCouplingMatrix_ee(i,j) = AvgCouplingMatrix_ee(i,j)+...
                            normal_e*normal_e' * (w'*(boundary_phi_e(:,i).* (boundary_phi_e(:,j))))*LineElemJacobian;
                    end
                end
%             end
        
        % only treat interior terms here, we do the BCs elsewhere
        if neighbor_id ~= -1
            [v, V_n] = GetFaceLocalDofs(neighbor_face_id);
            
            if n_dim ==2
                boundary_phi_n = phi(:,:);
            else
                boundary_phi_n = boundary_phi(neighbor_face_id,:);
            end
            %         % assemble boundary condition terms
            %
            %             % outflow BC on q2 on both boundaries
            %             % impose a one-sided flux, f(q^-), from inside the domain
            %             % using AvgCouplingMatrix_ee
            %             % We need to double AvgCouplingMatrix_ee here since we're
            %             % not "averaging" here
            %             AvgCouplingMatrix(dofs_q1_e,dofs_q1_e) = ...
            %                 AvgCouplingMatrix(dofs_q1_e,dofs_q1_e) + 2*AvgCouplingMatrix_ee;
            %
            %         else
            %         % assemble element interior terms
            
            % we only do assembly if elem_id > neighbor_id
            % to avoid assembling terms twice
            if elem_id < neighbor_id
                continue;
            end
            
            % get the shape function values on the neighbor
 
            
            % build the dense matrix for the element-neighbor coupling
            AvgCouplingMatrix_en = zeros(n_local_dofs,n_local_dofs);
            
%                 for node = 1:nodes_per_face
                    for i=1:n_local_dofs
                        for j=1:n_local_dofs
                            AvgCouplingMatrix_en(i,j) = AvgCouplingMatrix_en(i,j)+...
                                normal_n*normal_e' * 0.5 * (w'*(boundary_phi_e(:,i).* flipud(boundary_phi_n(:,j))))*LineElemJacobian;
                        end
                    end
%                 end
            
            
            % build the dense matrix for the neighbor-neighbor coupling
            AvgCouplingMatrix_nn = zeros(n_local_dofs,n_local_dofs);

                    for i=1:n_local_dofs
                        for j=1:n_local_dofs
                            AvgCouplingMatrix_nn(i,j) = AvgCouplingMatrix_nn(i,j)+...
                                normal_n*normal_n' * 0.5 *  (w'*(flipud(boundary_phi_n(:,i)).* flipud(boundary_phi_n(:,j))))*LineElemJacobian;
                        end
                    end
            % build the dense matrix for the neighbor-element coupling
            AvgCouplingMatrix_ne = zeros(n_local_dofs,n_local_dofs);
                    for i=1:n_local_dofs
                        for j=1:n_local_dofs % previously all i's and j's were 1:n_local_dofs
                            AvgCouplingMatrix_ne(i,j) = AvgCouplingMatrix_ne(i,j)+...
                                normal_n*normal_e' *  (w'*(flipud(boundary_phi_n(:,i)) .* flipud(boundary_phi_e(:,j))))*LineElemJacobian;
                        end
                    end
            
            % add the dense matrices to the sparse matrix!
                for i = 1:n_vars
                    
                    dofs_qi_e = DofMap{i}(elem_id,:);
                    dofs_qi_n = DofMap{i}(neighbor_id,:);
                    
                        
                        
                    JumpCoupMat(dofs_qi_e,dofs_qi_e) = JumpCoupMat(dofs_qi_e,dofs_qi_e) + ...
                        AvgCouplingMatrix_ee(:,:);
                    JumpCoupMat(dofs_qi_e,dofs_qi_n) = JumpCoupMat(dofs_qi_e,dofs_qi_n) + ...
                        AvgCouplingMatrix_en(:,:);
                    JumpCoupMat(dofs_qi_n,dofs_qi_n) = JumpCoupMat(dofs_qi_n,dofs_qi_n) + ...
                        AvgCouplingMatrix_nn(:,:);
                    JumpCoupMat(dofs_qi_n,dofs_qi_e) = JumpCoupMat(dofs_qi_n,dofs_qi_e) + ...
                        AvgCouplingMatrix_ne(:,:);

                end
        end
    end
end
end