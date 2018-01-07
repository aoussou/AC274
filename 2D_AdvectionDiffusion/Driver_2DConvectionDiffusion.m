%% Final Project AC 274
% John Aoussou and Sydney Sroka
% May 16, 2014

% This driver solves a diffusion problem in 2D for advection_velocity = 0 and an 
% advection problem for epsilon (diffusion coefficient) =0. The
% advection-diffusion problem is solved for advection_velocity and 
% epsilon ~= 0


clear all; clc;close all
addpath('./Plotting/');
addpath('./distmesh/');

Globals2D;

%% set problem parameters
n_vars = 1;
advection_velocity       = [0 0];
delta_t                  = 0.00001;      % time step size
n_tsteps                 = 100000;       % number of time steps
p                        = 2;           % basis function polynomial order 0,1,2
epsilon                  = 1;       % diffusion coefficient

initial_condition_string = 'Gaussian';  % initial condition (e.g. Gaussian, Square)
bndry_condition_string   = 'NoPenetration';      % boundary condition, NoPenetration or Outflow
flux_function_string     = 'advection'; %diffusion or advection



%% initialization =========================================================

%-----------initialize domain----------------------------------------------
% fd=@(p)(drectangle(p,0,1,0,1));
% [VX,EToN]=distmesh2d(fd,@huniform,.05,[0 0;1,1],[0,0; 0 1; 1,1; 1 0]);
h = 0.1;
     fd=@(p) drectangle(p,-.5,0.5,-0.5,0.5);
     [VX,EToN]=distmesh2d(fd,@huniform,h,[-0.5,-0.5;0.5,0.5],[-0.5,-0.5;-0.5,0.5;0.5,-0.5;0.5,0.5]);
% CFL Condition

if epsilon*delta_t/(h^2)>.249999
    disp('Warning: CFL Violated')
    return
end
      
%-----------initialize Matrix Assembly-------------------------------------

ElemNeighbors = GetConnectivity_2D;         %Generates a Kx3 matrix; row i <=> element, 
                                            %col j <=> neighbor element on the jth face

CreateDofMap;                               % generate the degree-of-freedom map, DofMap
                                            % global dofs = DofMap{var}(elem_id,elementlocal dofs)
                                    
AssembleMassMatrix;                         %OUTPUT: sparse n_tot_dofsxn_total_dofs

AssembleAdvectionMatrix;                    %OUTPUT: sparse n_tot_dofsxn_total_dofs, {Advec_x_1 ; Advec_x_1 }

AssembleAvgCouplingMatrix;                  %OUTPUT: sparse n_tot_dofsxn_total_dofs, {AvgCoup_x_1 ; AvgCoup_x_1 }

AssembleJumpCouplingMatrix;                 %OUTPUT: sparse n_tot_dofsxn_total_dofs

AssembleBoundaryConditionMatrix;            %OUTPUT: sparse n_tot_dofsxn_total_dofs, {AvgCoup_x_1 ; AvgCoup_x_1 }

BoundaryConditionDofs;                      %OUTPUT: dofs_D [num Dirichlet dofs,1] and dofs_N [num Neumann dofs,1]

% set the initial conditions 
initial_solution_coeffs = InitializeSolutionCoeffs;

% initialize solution_coeffs vector
solution_coeffs = initial_solution_coeffs;

time = 0;
plot_q;
title('100*exp(-100*((x1).^2+(x2).^2))')
drawnow
colorbar

%% integrate in time ------------------------------------------------------
tic
for tstep = 1:n_tsteps
    clf
    time = tstep*delta_t;
    
    %RK4
    k1 = M\AssembleRHS(solution_coeffs,time);
    k2 = M\AssembleRHS(solution_coeffs + 0.5*delta_t*k1,time);
    k3 = M\AssembleRHS(solution_coeffs + 0.5*delta_t*k2,time);
    k4 = M\AssembleRHS(solution_coeffs + delta_t*k3,time);
    
    solution_coeffs = solution_coeffs + ...
        delta_t/6 * (k1+2*k2+2*k3+k4);
    
    % Impose Dirichlet Conditions
    solution_coeffs= ImposeDirichlet(solution_coeffs);
    
    if mod(tstep,50) == 0

        plot_q;
        title(sprintf('DG Implementation of Diffusion Equation at time = %g \n(p=2, dt=0.001,eps=0.01,h=0.05)',time))
        colorbar
        drawnow
        
    end
end
toc





