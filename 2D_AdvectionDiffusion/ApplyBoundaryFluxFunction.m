function f_boundary = ApplyBoundaryFluxFunction(q,time)
% compute the flux for various problems

Globals2D;


if p==0
    n_local_dofs = 1;
    calN = n_local_dofs*size(EToN,1);
    n_total_dofs = calN*n_vars;
else
    n_local_dofs = 3*p;
    calN = n_local_dofs*size(EToN,1);
    
    n_total_dofs = calN*n_vars;
end


% reshape q so we can apply the flux function column by column
q_reordered = zeros(n_vars,calN);

for var = 1:n_vars
    first_index = calN*(var-1) + 1;
    last_index = var*calN;
    q_reordered(var,:) = q(first_index:last_index)';
end

% initialize flux to 0
f_reordered = 0*q_reordered;


% calculate boudary flux for various problems
switch bndry_condition_string
    
    case 'None'
        % do nothing
        
    case 'Problem1a'
        %zero inflow bcs on the inflow, outflow, do nothing.
        
        
%         
%         q_reordered(1,dof_to_0)=0;
%         f_reordered = q_reordered;
%         % finally put f into the standard ordering
        f_boundary = zeros(n_total_dofs,2);
%         for var = 1:n_vars
%             first_index = calN*(var-1) + 1;
%             last_index = var*calN;
%             f_boundary(first_index:last_index) = f_reordered(var,:);
%         end
        
    case 'Problem1d'
        
        [BExt] = findDomainBoundaries;% first col element, second col side
        %         u * n < 0 => inflow
        for b_elem = 1:size(BExt,1)
            [~, face_local_dofs]=GetFaceLocalDofs(BExt(b_elem,2));
            w1 = VX(EToN(BExt(b_elem,1),face_local_dofs(1)),:);
            w2 = VX(EToN(BExt(b_elem,1),face_local_dofs(end)),:);
            T = (w2-w1)/norm(w2-w1,2);
            normal_e = [T(2), -T(1)];
            midPoint = 0.5*[w1(1)+w2(1);w1(2)+w2(2)];
            
            %                 if (normal_e*[wave_speed{1}(midPoint(1),midPoint(2));wave_speed{2}(midPoint(1),midPoint(2))])<0 %inflow
            if (normal_e*[midPoint(2);-midPoint(1)])<0 %inflow
                dof_to_0 = DofMap{i}(BExt(b_elem,1),face_local_dofs);
                q_reordered(i,dof_to_0)=0;
                
            end
        end
        f_reordered = q_reordered;
        % finally put f into the standard ordering
        f_boundary = zeros(n_total_dofs,1);
        for var = 1:n_vars
            first_index = calN*(var-1) + 1;
            last_index = var*calN;
            f_boundary(first_index:last_index) = f_reordered(var,:);
        end
    case 'Problem2'
        %Vector comes out in REORDERED format

% apply fluxfunction to the q's first THEN apply the flux function to that to get the right vector

        rho0 = 1;K0 = 1;
%         I have some dof and I want the element
% do x and y seperately because I know the direction of the normal vectors
% in this case
% 
%         [b_ele, b_ele_dof] = find(DofMap{1}==dof_NoFlux_x);
% 
%         [~, face_local_dofs]=GetFaceLocalDofs(BExt(b_elem,2));
%         w1 = VX(EToN(BExt(b_elem,1),face_local_dofs(1)),:);
%         w2 = VX(EToN(BExt(b_elem,1),face_local_dofs(end)),:);
%         T = (w2-w1)/norm(w2-w1,2);
%         normal = [T(2), -T(1)];
% 
%         q_b(1,dof_NoFlux) = 0;
%         q_b(2,dof_NoFlux) = (1/rho0)*q(1,dof_NoFlux)*normal(1);
%         q_b(3,dof_NoFlux) = (1/rho0)*q(1,dof_NoFlux)*normal(2);
% 
% 
% 
% 
%         [BExt] = findDomainBoundaries;% first col element, second col side
%         f_noFlux = spalloc(3,calN,n_vars*length(dof_NoFlux));
%         K0 = 1; rho0 = 1;
% fb_x1 = zeros(3,calN);        
% fb_x2 = zeros(3,calN);        

[BExt] = findDomainBoundaries;% first col element, second col side

fb_x1 = zeros(n_vars,calN);
fb_x2 = zeros(n_vars,calN);


% for b_elem = 1:size(BExt,1)
% %compute the normal vector
%         [~, face_local_dofs]=GetFaceLocalDofs(BExt(b_elem,2));
%         w1 = VX(EToN(BExt(b_elem,1),face_local_dofs(1)),:);
%         w2 = VX(EToN(BExt(b_elem,1),face_local_dofs(end)),:);
%         T = (w2-w1)/norm(w2-w1,2);
%         normal = [T(2), -T(1)];
% %compute the boundary condition vector
%         belem_g_dofs = DofMap{1}(b_elem,face_local_dofs);
% 
%         fb_x1(1,belem_g_dofs) = 0;
%         fb_x2(1,belem_g_dofs) = 0;
% 
%         fb_x1(2,belem_g_dofs) = (1/rho0)*q_reordered(1,belem_g_dofs)*(normal(1));
%         fb_x2(2,belem_g_dofs) = 0;
% 
%         fb_x1(3,belem_g_dofs) = 0;
%         fb_x2(3,belem_g_dofs) = (1/rho0)*q_reordered(1,belem_g_dofs)*(normal(2));
% 
%         fb_x1(2,dof_NoFlux) = (1/rho0)*q_reordered(1,dof_NoFlux)*(normal(i));
%         fb_x2(3,dof_NoFlux) = (1/rho0)*q_reordered(1,dof_NoFlux)*(normal(i));
% 
% 
% end

% for i = 1:length(dof_NoFlux)
% 
%         fb_x1(2,dof_NoFlux) = (1/rho0)*q_reordered(1,dof_NoFlux)*(normal_noFlux(i,1));
%         fb_x2(3,dof_NoFlux) = (1/rho0)*q_reordered(1,dof_NoFlux)*(normal_noFlux(i,2));
% end




f_boundary = [fb_x1(1,:)' fb_x2(1,:)' ; fb_x1(2,:)' fb_x2(2,:)' ; fb_x1(3,:)' fb_x2(3,:)'];

%{

        dof_to_0   = [];
        dof_NoFlux_x = [];
        dof_NoFlux_y = [];
        [BExt] = findDomainBoundaries;% first col element, second col side
        %         u * n < 0 => inflow
        tol = 1e-6;
        %         inc = 1;
        for i = 1%I'm pretty sure that this works to repeat for the reordered variables
            for b_elem = 1:size(BExt,1)
                [~, face_local_dofs]=GetFaceLocalDofs(BExt(b_elem,2));
                w1 = VX(EToN(BExt(b_elem,1),face_local_dofs(1)),:);
                w2 = VX(EToN(BExt(b_elem,1),face_local_dofs(end)),:);
                T = (w2-w1)/norm(w2-w1,2);
                normal = [T(2), -T(1)];
                if abs(w1(1)+1)<tol && abs(w2(1)+1)<tol
                    dof_to_0 = [dof_to_0;(DofMap{1}(BExt(b_elem,1),face_local_dofs))'];
                    %I know all my boundaries have normal vectors [+/- 1 (0), 0 (+/-1)]
                else
                    dof_NoFlux = [dof_NoFlux_y;(DofMap{1}(BExt(b_elem,1),face_local_dofs))'];
                    q_noFlux(1,dof_NoFlux(dof)) = 0;
                    q_noFlux(2,dof_NoFlux(dof)) = (1/rho0)*q_reordered(1,dof_NoFlux(dof))*normal(1);%row vector
                    q_noFlux(3,dof_NoFlux(dof)) = (1/rho0)*q_reordered(1,dof_NoFlux(dof))*normal(2); 
                end
            end
        end
%}
        
        %         for dof = 1:length(dof_NoFlux)
        %                 %find the element that this dof corresponds to
        %
        %                 %I need to do an inverse DofMap
        %                 [b_ele, b_ele_loc_dof] = find(DofMap{1}==dof_NoFlux(dof));
        %
        %                 w1 = VX(EToN(b_ele,b_ele_loc_dof),:);
        %                 w2 = VX(EToN(BExt(b_elem,1),face_local_dofs(end)),:);
        %                 T = (w2-w1)/norm(w2-w1,2);
        %                 normal = [T(2), -T(1)];
        %                 %store the variables in columns so that your col-wise indexing
        %                 %works
        %                 %             q_reordered = q_reordered';
        %
        %                 f_noFlux(1,dof_NoFlux(dof)) = 0;
        % % K0*(q_reordered(2,dof_NoFlux(b_elem))*normal(1)+q_reordered(3,dof_NoFlux(b_elem))*normal(2));
        %                 f_noFlux(2,dof_NoFlux(dof)) = (1/rho0)*q_reordered(1,dof_NoFlux(dof))*normal(1);%row vector
        %                 f_noFlux(3,dof_NoFlux(dof)) = (1/rho0)*q_reordered(1,dof_NoFlux(dof))*normal(2);
        % %             end
        %         end
        
        
        
        
        
        
        
        %         for b_elem = 1:size(BExt,1)
        % %             if ~ismember(dof_NoFlux(b_elem),dof_to_0)
        %
        %                 [~, face_local_dofs]=GetFaceLocalDofs(BExt(b_elem,2));
        %                 w1 = VX(EToN(BExt(b_elem,1),face_local_dofs(1)),:);
        %                 w2 = VX(EToN(BExt(b_elem,1),face_local_dofs(end)),:);
        %                 T = (w2-w1)/norm(w2-w1,2);
        %                 normal = [T(2), -T(1)];
        %                 %store the variables in columns so that your col-wise indexing
        %                 %works
        %                 %             q_reordered = q_reordered';
        %
        %                 f_noFlux(1,dof_NoFlux(b_elem)) = 0;
        % % K0*(q_reordered(2,dof_NoFlux(b_elem))*normal(1)+q_reordered(3,dof_NoFlux(b_elem))*normal(2));
        %                 f_noFlux(2,dof_NoFlux(b_elem)) = (1/rho0)*q_reordered(1,dof_NoFlux(b_elem))*normal(1);%row vector
        %                 f_noFlux(3,dof_NoFlux(b_elem)) = (1/rho0)*q_reordered(1,dof_NoFlux(b_elem))*normal(2);
        % %             end
        %         end
        %build up the no flux vector then set the dirichlet to zero
        %with dof_to_zero
        %         f_boundary = zeros(n_total_dofs,1);
        %         for var = 1:n_vars
        %             first_index = calN*(var-1) + 1;
        %             last_index = var*calN;
        %             f_boundary(first_index:last_index) = f_reordered(var,:);
        %         end
        %         f_boundary = [f_noFlux(1,:)' f_noFlux(2,:)' f_noFlux(3,:)'];
        
        %{
%         [~, ~, ~, ~,BExt] = findDomainBoundaries;% first col element, second col side
%         %first make all of the boundary zero, then add the pressure forcing
%         inc = 1;%counter for dofs to set with dirichlet-ish conditions
%         f = cell(n_vars,1);
%         for i = 1:n_vars
%             f{i} = zeros(1,calN);
%         end
%
%         for b_elem = 1:size(BExt,1)
%             [~, face_local_dofs]=GetFaceLocalDofs(BExt(b_elem,2));
%             w1 = VX(EToN(BExt(b_elem,1),face_local_dofs(1)),:);
%             w2 = VX(EToN(BExt(b_elem,1),face_local_dofs(end)),:);
%             T = (w2-w1)/norm(w2-w1,2);
%             normal_e = [T(2), -T(1)];
%
%             global_dofs(1,:) = DofMap{1}(BExt(b_elem,1),face_local_dofs);
%             global_dofs(2,:)  = DofMap{2}(BExt(b_elem,1),face_local_dofs);
%             global_dofs(3,:)  = DofMap{3}(BExt(b_elem,1),face_local_dofs);
%
%             f{1} = f{1}+K0*(q_reordered(2,global_dofs(1,:))*normal_e(1)+q_reordered(3,global_dofs(1,:))*normal_e(2));
%             f{2} = f{2}+1/rho0*q_reordered(1,global_dofs(2,:))*normal_e(1);
%             f{3} = f{3}+1/rho0*q_reordered(1,global_dofs(3,:))*normal_e(2);
%             %             for i = 2:3%the u and v compoonants
%             %                 q_reordered(i,dof_to_0)=0;
%             %             end
%
%             if abs(w1(1)+1)<1e-6  && abs(w2(1)+1)<1e-6
%                 dof_to_set(length(face_local_dofs)*(inc-1)+1:length(face_local_dofs)*(inc),1) = DofMap{1}(BExt(b_elem,1),face_local_dofs);
%                 inc = inc+1;
%             end
%             keyboard
%         end
%         f{1}(dof_to_set)=0.01*sin(20*time);%we are using the dof_to_set dofs from the first loop
%         for var = 1:n_vars
%             first_index = calN*(var-1) + 1;
%             last_index = var*calN;
%
%             f_boundary(first_index:last_index) = f{var};
%         end
        %}
    case 'Problem2b'
        
        [BExt] = findDomainBoundaries;% first col element, second col side
        f_noFlux = zeros(3,calN);
        K0 = 1; rho0 = 1;
        %
        %         for b_elem = 1:size(BExt,1)
        %             if ~ismember(dof_NoFlux_x(b_elem),dof_to_0)
        %                 if ~ismember(dof_NoFlux_y(b_elem),dof_to_0)
        %
        %                     [~, face_local_dofs]=GetFaceLocalDofs(BExt(b_elem,2));
        %                     w1 = VX(EToN(BExt(b_elem,1),face_local_dofs(1)),:);
        %                     w2 = VX(EToN(BExt(b_elem,1),face_local_dofs(end)),:);
        %                     T = (w2-w1)/norm(w2-w1,2);
        %                     normal = [T(2), -T(1)];
        %                     %store the variables in columns so that your col-wise indexing
        %                     %works
        %                     %             q_reordered = q_reordered';
        f_noFlux(2,dof_NoFlux_x) = 0;
        f_noFlux(3,dof_NoFlux_y) = 0;
        %                 end
        %             end
        %         end
        %build up the no flux vector then set the dirichlet to zero
        %with dof_to_zero
        %         f_boundary = zeros(n_total_dofs,1);
        %         for var = 1:n_vars
        %             first_index = calN*(var-1) + 1;
        %             last_index = var*calN;
        %             f_boundary(first_index:last_index) = f_reordered(var,:);
        %         end
        %         f_boundary = [f_noFlux(1,:)' f_noFlux(2,:)' f_noFlux(3,:)'];
        f_boundary = f_noFlux;
    case 'Problem3'
        [BExt] = findDomainBoundaries;% first col element, second col side
        f_noFlux = spalloc(3,calN,n_vars*length(dof_NoFlux));
        K0 = 1; rho0 = 1;
        Pressure = (gamma - 1)*(q_reordered(4,:)-0.5*(q_reordered(2,:).^2+q_reordered(3,:).^2)./q_reordered(1,:));
        
        for b_elem = 1:size(BExt,1)
            if ~ismember(dof_NoFlux(b_elem),dof_to_0)
                
                [~, face_local_dofs]=GetFaceLocalDofs(BExt(b_elem,2));
                w1 = VX(EToN(BExt(b_elem,1),face_local_dofs(1)),:);
                w2 = VX(EToN(BExt(b_elem,1),face_local_dofs(end)),:);
                T = (w2-w1)/norm(w2-w1,2);
                normal = [T(2), -T(1)];
                %store the variables in columns so that your col-wise indexing
                %works
                %             q_reordered = q_reordered';
                
                f_noFlux(1,dof_NoFlux(b_elem)) = zeros(calN,1);
                f_noFlux(2,dof_NoFlux(b_elem)) = p(b_elem)*normal(1);
                f_noFlux(3,dof_NoFlux(b_elem)) = p(b_elem)*normal(2);
                f_noFlux(4,dof_NoFlux(b_elem)) = zeros(calN,1);
                
            end
        end
        
        
end
end






