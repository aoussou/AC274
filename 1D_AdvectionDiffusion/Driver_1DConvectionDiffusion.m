%% Final Project AC 274
% John Aoussou and Sydney Sroka
% May 15, 2014

% This driver solves a diffusion problem for advection_velocity = 0 and an
% advection problem for epsilon (diffusion coefficient) =0. The
% advection-diffusion problem is solved for advection_velocity and
% epsilon ~= 0


clear all; clc;close all
Globals1D;

%% set problem parameters

a                        = -1;        % mesh min x
b                        = 1;         % mesh max x
K                        = 75;        % mesh number of elements
advection_velocity       = 0;
delta_t                  = 0.0001;    % time step size
n_tsteps                 = 10000;    % number of time steps
p                        = 2;         % basis function polynomial order 0,1,2
epsilon                  = .1;       % diffusion coefficient
periodic_bcs             = 1;
ghost_cell_flag          = 1;         % Use ghost cells at either edge of the boundary

initial_condition_string = 'Gaussian';  % initial condition (e.g. Gaussian, Square)
bndry_condition_string   = 'NoPenetration';  % boundary condition for diffusion equation, NoPenetration or Outflow
flux_function_string     = 'advection'; %diffusion or advection
%% Check CFL violation based on input parameters
dx = (b-a)/K;
CFL = advection_velocity*delta_t/dx;
if CFL > .9999
    disp('Warning: CFL VIOLATED')
    return
end
CFL_eps = epsilon*delta_t/(dx^2);
if CFL_eps > .49999
    disp('Warning: CFL EPS VIOLATED')
    return
end

%% initialization ---------------------------------------------------------
if ghost_cell_flag
    a = a-dx; b=b+dx;K = K+2;
end

% generate a mesh with K equally sized elements
[VX, EToN] = MeshGen1D(a, b, K);

% generate the Neighbor array that defines element connectivity
ElemNeighbors = GetConnectivity(EToN);

if periodic_bcs
    ElemNeighbors(1) = K;
    ElemNeighbors(end) = 1;
end

% generate the degree-of-freedom map, DofMap
CreateDofMap;

% assemble our discretization matrices
AssembleMassMatrix;
AssembleAdvectionMatrix;
AssembleAvgCouplingMatrix;
AssembleJumpCouplingMatrix;
AssembleBoundaryConditionMatrix;

% set the initial conditions
initial_solution_coeffs = InitializeSolutionCoeffs(initial_condition_string);

% initialize solution_coeffs vector
solution_coeffs = initial_solution_coeffs;
% 
time = 0;
figure(1)
plot_DG_solution(solution_coeffs,time);
hold off

%% integrate in time ------------------------------------------------------
tic
% for epsilon = [0.01 0.1 1]
    
    for tstep = 1:n_tsteps
        
        time = tstep*delta_t;
        
        
        %RK4
        % for an RHS that is time independant
        k1 = M\AssembleRHS(solution_coeffs,time);
        k2 = M\AssembleRHS(solution_coeffs + 0.5*delta_t*k1,time);
        k3 = M\AssembleRHS(solution_coeffs + 0.5*delta_t*k2,time);
        k4 = M\AssembleRHS(solution_coeffs + delta_t*k3,time);
        
        solution_coeffs = solution_coeffs + ...
            delta_t/6 * (k1+2*k2+2*k3+k4);
        
        % Impose Dirichlet Conditions
        solution_coeffs= ImposeDirichlet(solution_coeffs);
        %
                T = n_tsteps*delta_t;
%                 if time == T
                if mod(tstep,500) == 0
%                     if p ==0 ; c = 'b' ; elseif p ==1 ; c = 'r' ; elseif p ==2 ; c = 'g' ; end
                    figure(1)
                    plot_DG_solution(solution_coeffs,time);
                    title(sprintf('Advection of a Gaussian t=%g   (K=%g,dt = 0.0001,epsilon = %g)',time,K-2,epsilon))
                    ylim([-.1 1])
                    drawnow
                    hold off
                end
%         
        %{
        if epsilon == 1
            c= 'b';
            if norm(double(time==[T]))
                h = figure(2);
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
                    
                    
                    h3 = plot(local_x, local_qh,c,'linewidth',2);
                    
                    
                    if(elem_id==elem1)
                        ylim([-0.5,1.5]);
                        hold on
                    end
                end
                
                h_legend =legend([h1 h2 h3],'epsilon = 0.01','epsilon = 0.1','epsilon = 1');
                set(h_legend,'FontSize',14)
                title(sprintf('Diffusion of a Gaussian t=%g (p=%g,K=%g,dt =%g)',time,p,K-2,delta_t))
                drawnow
                hf = gcf;
                saveas(hf,sprintf('gaussian_t15_p0_3eps_a1_bneg1.png'))
            end
        elseif epsilon ==0.1
            c = 'r';
            if norm(double(time==[T]))
                h = figure(2);
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
                    
                    
                    h2 = plot(local_x, local_qh,c,'linewidth',2);
                    
                    
                    if(elem_id==elem1)
                        ylim([-0.5,1.5]);
                        hold on
                    end
                end
                
                drawnow
            end
        elseif epsilon == 0.01
            c= 'm';
            if norm(double(time==[T]))
                h = figure(2);
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
                    
                    
                    h1 = plot(local_x, local_qh,c,'linewidth',2);
                    
                    
                    if(elem_id==elem1)
                        ylim([-0.5,1.5]);
                        hold on
                    end
                end
                
                drawnow
                %         hold all
            end
        end
        
        %}
        
    end

% end
toc





