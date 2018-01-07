function plot_DG_solution(solution_coeffs,time,c)
% plot the DG solution specified by solution_coeffs
if nargin<3
    c='b';
end

Globals1D;

reference_pts = linspace(-1,1,10)';
n_reference_pts = length(reference_pts);
phi = GetPhi(p, reference_pts);

n_local_dofs = p+1;

%account for ghost elements
if ghost_cell_flag
    elem1 = 2;     elemK = K-1; 
else
    elem1 = 1;     elemK = K; 
end

for elem_id=elem1:elemK 
    local_dofs = DofMap(elem_id,:);
    
    local_qh = zeros(n_reference_pts,1);
    for j=1:n_local_dofs
        local_dof_j = local_dofs(j);
        local_qh = local_qh + ...
            phi(:,j)*solution_coeffs(local_dof_j);
    end
    
    local_x = GetPhysicalPoints(elem_id, reference_pts);

    
    plot(local_x, local_qh,c,'linewidth',2);

    
    if(elem_id==elem1)
%         ylim([-0.5,1.5]);
        hold on
    end
end


    
% hold off
% 
% title(['$t=',num2str(time),'$'],'interpreter','latex','fontsize',18)
% xlabel('$x$','interpreter','latex','fontsize',18)
