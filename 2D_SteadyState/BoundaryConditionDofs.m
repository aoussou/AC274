function BoundaryConditionDofs


global DofMap
global EToN
global VX
global bndry_condition_string
global dofs_N
global dofs_D
global p 
global normal_vector
global lineElemJac

% dofs_D => dofs with Dirichlet BCs
% dofs_N => dofs with Dirichlet BCs


        dofs_N = [];
                global dofs_l dofs_b dofs_t dofs_r

        [BExt] = findDomainBoundaries;% first col element, second col side
        tol = 1e-6;
        for b_elem = 1:size(BExt,1)
            [~, face_local_dofs]=GetFaceLocalDofs(BExt(b_elem,2));
            w1 = VX(EToN(BExt(b_elem,1),face_local_dofs(1)),:);
            w2 = VX(EToN(BExt(b_elem,1),face_local_dofs(end)),:);
            midPoint = 0.5*[w1(1)+w2(1);w1(2)+w2(2)];
            X = [w1(1) midPoint(1) w2(1)];
            Y = [w1(2) midPoint(2) w2(2)];
            tol = 1e-6;
            

            if abs(w1(1)+0.5)<tol && abs(w2(1)+0.5)<tol
                dofs_l = [dofs_l;(DofMap{1}(BExt(b_elem,1),face_local_dofs))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs))'];

            elseif abs(w1(2)+0.5)<tol && abs(w2(2)+0.5)<tol
                dofs_b = [dofs_b;(DofMap{1}(BExt(b_elem,1),face_local_dofs))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs))'];

            elseif abs(w1(1)-0.5)<tol && abs(w2(1)-0.5)<tol
                dofs_r = [dofs_r;(DofMap{1}(BExt(b_elem,1),face_local_dofs))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs))'];

            elseif abs(w1(2)-0.5)<tol && abs(w2(2)-0.5)<tol
                dofs_t = [dofs_t;(DofMap{1}(BExt(b_elem,1),face_local_dofs))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs))'];

            elseif abs(w1(1)+0.5)<tol
                dofs_l = [dofs_l;(DofMap{1}(BExt(b_elem,1),face_local_dofs(1)))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs(1)))'];

            elseif abs(w1(2)+0.5)<tol 
                dofs_b = [dofs_b;(DofMap{1}(BExt(b_elem,1),face_local_dofs(1)))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs(1)))'];

            elseif abs(w2(1)+0.5)<tol
                dofs_l = [dofs_l;(DofMap{1}(BExt(b_elem,1),face_local_dofs(end)))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs(end)))'];

            elseif abs(w2(2)+0.5)<tol 
                dofs_b = [dofs_b;(DofMap{1}(BExt(b_elem,1),face_local_dofs(end)))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs(end)))'];

            elseif abs(w1(1)-0.5)<tol
                dofs_r = [dofs_r;(DofMap{1}(BExt(b_elem,1),face_local_dofs(1)))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs(1)))'];

            elseif abs(w1(2)-0.5)<tol 
                dofs_t = [dofs_t;(DofMap{1}(BExt(b_elem,1),face_local_dofs(1)))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs(1)))'];

            elseif abs(w2(1)-0.5)<tol
                dofs_r = [dofs_r;(DofMap{1}(BExt(b_elem,1),face_local_dofs(end)))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs(end)))'];

            elseif abs(w2(2)-0.5)<tol 
                dofs_t = [dofs_t;(DofMap{1}(BExt(b_elem,1),face_local_dofs(end)))'];
                dofs_N = [dofs_N;(DofMap{1}(BExt(b_elem,1),face_local_dofs(end)))'];


            end
        end
        
            
end




