function M = AssembleMassMatrix
% Assemble the mass matrix

%Should be a "diagonal" matrix s.t. no coupling bewteen elements for DG
Globals2D;
global p
global EToN
global n_vars
global DofMap
n_vars = 1;
if p==0
    n_local_dofs = 1;
    N_quad_pts = 1;%number of GQ intgration points
else
    n_local_dofs = 3*p;
    N_quad_pts = p+2;
end
n_total_dofs = n_local_dofs*size(EToN,1);
[x_quad,w] = GaussQuad2D(N_quad_pts);       %OUTPUT: w = N_quad_ptsx1, 
                                            % x_quad = N_quad_ptsx2 for [xi, eta]
phi = GetPhi2D(x_quad);                     %OUTPUT: N_quad_pts x n_local_dofs 
                                            % row i <=> point i, col j <=> phi j

MassMatrix = spalloc(n_total_dofs, n_total_dofs, ...
    n_vars*size(EToN,1)*n_local_dofs*n_local_dofs);

for elem_id = 1:size(EToN,1)
    ElemJacobian = det(GetElemJacobian2D(elem_id));
    phys_points = GetPhysicalPoints2D(elem_id, x_quad);
    
    
    for i = 1:n_local_dofs
        for j = 1:n_local_dofs
            LocalMassMatrix(i,j) = w'*(phi(:,i).*phi(:,j))*ElemJacobian; %Gaussian integration
        end
    end
    
    for i = 1:n_vars
        dofs_qi = DofMap{i}(elem_id,:);
        MassMatrix(dofs_qi,dofs_qi) = LocalMassMatrix;
    end

end
M = MassMatrix;