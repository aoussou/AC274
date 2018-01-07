function AssembleLoadVector

global DOF
global VX
global p
global LoadVector
global ElemNeighbors
global bndry_condition_string
global ndofs_p


global dofs_t bele_t
global dofs_r bele_r
global dofs_b bele_b
global dofs_l bele_l
global dofs_c1 bele_c1
global dofs_c2 bele_c2
global dofs_c3 bele_c3
global dofs_c4 bele_c4


BExt = findDomainBoundaries;
n_total_dofs = size(VX,1);

n_local_dofs = 3*p;
n_total_dofs = size(VX,1);

[x_quad2D, w2D] = GaussQuad2D(4);
phi = GetPhi2D(x_quad2D);%d/d_xi in 1st col, x/d_eta in 2nd col, diff phi per 3D card

ElemLoadVector = zeros(n_local_dofs,1);
LoadVector = zeros(n_total_dofs,1);

gN = zeros(n_total_dofs,1);

[kb sb] = find(ElemNeighbors == -1);%the element on the boundary kb and the side of that element sb
BExt=[kb sb];% 2D matrix with [boundary element, face id] in each row

switch bndry_condition_string % create the f vector and the gN vector
    case '1b'
        f = ones(size(VX,1),1);
        
    case '1c'
        f = -10*ones(size(VX,1),1);
        
    case '1d'
        f = ones(size(VX,1),1);
        
        dofs_Nb = [dofs_b];
        gN(dofs_Nb) = 1;
        dofs_Nc = [ dofs_c1;dofs_c2;dofs_c3;dofs_c4];
        gN(dofs_Nc) = 0;
        dofs_N = [dofs_Nb ; dofs_Nc];
        BExt = [bele_b; bele_c1; bele_c2; bele_c3; bele_c4];
        
    case '1e'
        gR = zeros(n_total_dofs,1);
        
        dofs_Nc = [ dofs_c1;dofs_c2;dofs_c3;dofs_c4];
        gN(dofs_Nc) = 0;
        dofs_Nr = [ dofs_r];
        gR(dofs_Nr) = 1;
        f = zeros(size(VX,1),1);
        BExt = [bele_r; bele_c1; bele_c2; bele_c3; bele_c4];
        
        [v(1,:), V_e(1,:)] = GetFaceLocalDofs(1);
        [v(2,:), V_e(2,:)] = GetFaceLocalDofs(2);
        [v(3,:), V_e(3,:)] = GetFaceLocalDofs(3);
        
        [x_quad, w] = GaussQuad(p+1);
        w=w';
        ref_nodes = RefNodeLocations;
        
        for belem_id = 1:size(bele_r,1)
            ElemLoadVector = zeros(3,1);
            face = BExt(belem_id,2);
            
            dofs_e = DOF(bele_r(belem_id,1),:);
            G = segmap(x_quad,ref_nodes(v(face,:),:));
            phi_N = GetPhi2D(G);%phi in columns, points in rows
            w1 = VX(DOF(bele_r(belem_id),v(bele_r(belem_id,2),1)),:);% get the local nodes on each face to get the coordinates
            w2 = VX(DOF(bele_r(belem_id),v(bele_r(belem_id,2),2)),:);
            LineElemJacobian = 0.5*sqrt((w2-w1)*(w2-w1)');
            for i = 1:3
                %             for i = 1:n_local_dofs
                ElemLoadVector(i) = (gR(dofs_e(V_e(bele_r(belem_id,2),i))))*phi_N(:,(V_e(bele_r(belem_id,2),i)))'*w*LineElemJacobian;
            end
            LoadVector(dofs_e(V_e(bele_r(belem_id,2),:))) = ElemLoadVector + LoadVector(dofs_e(V_e(bele_r(belem_id,2),:)));
        end
        
        
    case '2a'
        %         f = exp(-100*((VX(:,1)-0.25).^2+(VX(:,2)+0.25).^2));
        
    case '2b'
        f = exp(-100*((VX(:,1)-0.25).^2+(VX(:,2)+0.25).^2));
        dofs_N = [dofs_l;dofs_t];
        gN(dofs_N) = 0;
        BExt = [bele_l; bele_t];
end



%% integral of f over Tk =================================================

for elem_id = 1:size(DOF,1)
    ElemLoadVector = zeros(6,1);
    
    elemJac = GetElemJacobian2D(elem_id);
    dofs_e = DOF(elem_id,:);
    
    phys_pts = GetPhysicalPoints2D(elem_id,x_quad2D);
    for ii = 1:n_local_dofs
        switch bndry_condition_string
 
            case {'2a','2b'}  
                
                ElemLoadVector(ii) = (exp(-100*((phys_pts(1,:)-0.25).^2+(phys_pts(2,:)+0.25).^2))).*phi(:,ii)'*w2D*det(elemJac);
                
            otherwise
                ElemLoadVector(ii) = f(dofs_e(ii))*phi(:,ii)'*w2D*det(elemJac);
        end
    end
    LoadVector(dofs_e) = ElemLoadVector + LoadVector(dofs_e);
end

%% Neumann BC load vector =================================================

[v(1,:), V_e(1,:)] = GetFaceLocalDofs(1);
[v(2,:), V_e(2,:)] = GetFaceLocalDofs(2);
[v(3,:), V_e(3,:)] = GetFaceLocalDofs(3);

[x_quad, w] = GaussQuad(p+1);
w=w';
ref_nodes = RefNodeLocations;

for belem_id = 1:size(BExt,1)
    ElemLoadVector = zeros(3,1);
    face = BExt(belem_id,2);
    %     if BExt(belem_id,1) == 40
    %         keyboard
    %     end
    
    dofs_e = DOF(BExt(belem_id,1),:);
    G = segmap(x_quad,ref_nodes(v(face,:),:));
    phi_N = GetPhi2D(G);%phi in columns, points in rows
    w1 = VX(DOF(BExt(belem_id),v(BExt(belem_id,2),1)),:);% get the local nodes on each face to get the coordinates
    w2 = VX(DOF(BExt(belem_id),v(BExt(belem_id,2),2)),:);
    LineElemJacobian = 0.5*sqrt((w2-w1)*(w2-w1)');
    for i = 1:3
        %             for i = 1:n_local_dofs
        ElemLoadVector(i) = (gN(dofs_e(V_e(BExt(belem_id,2),i))))*phi_N(:,(V_e(BExt(belem_id,2),i)))'*w*LineElemJacobian;
    end
    LoadVector(dofs_e(V_e(BExt(belem_id,2),:))) = ElemLoadVector + LoadVector(dofs_e(V_e(BExt(belem_id,2),:)));
end



