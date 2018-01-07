


addpath('./distmesh/')

%% Driver  Repeat steady state advection-diffusion

clear all; clc;
close all

Globals2D;
global dofs_l dofs_b dofs_t dofs_r

%% set problem parameters

advection_velocity       = [0,0];
p                        = 2;        % basis function polynomial order

tau                      = 1;      % Penalty term coefficient ( >0)
epsilon                  = 1;        % Diffusion coefficient

% f                        = @(x,y)-sin(y); % Forcing Function / Source Term (fxn of x)
f                      = @(x,y)-sin(x); % Forcing Function / Source Term (fxn of x)
Fdirac                   = 0;
% f                        = @(x,y) -exp(-100*(x.^2+y.^2));

ex = 1; %Diffusion in the x direction
ey = 1; %Diffusion in the x direction
 

%% initialization =========================================================

%-----------initialize domain----------------------------------------------
% fd=@(p)(drectangle(p,0,1,0,1));
% [VX,EToN]=distmesh2d(fd,@huniform,.05,[0 0;1,1],[0,0; 0 1; 1,1; 1 0]);
h = 0.06;
     fd=@(p) drectangle(p,-.5,0.5,-0.5,0.5);
%           [VX,EToN]=distmesh2d(fd,@huniform,h,[-0.5,-0.5;0.5,0.5],[-0.5,-0.5;-0.5,0.5;0.5,-0.5;0.5,0.5]);
figure(3)
% adding a point at zero to apply the dirac
     [VX,EToN]=distmesh2d(fd,@huniform,h,[-0.5,-0.5;0.5,0.5],[-0.5,-0.5;-0.5,0.5;0.5,-0.5;0.5,0.5;0 0]);

% CFL Condition
if epsilon*delta_t/(h^2)>.249999
    disp('Warning: CFL Violated')
    return
end

% VX = [VX(:,2),VX(:,1)];
      
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


%%

n_local_dofs = 3*p;
n_total_dofs = size(EToN,1)*n_local_dofs;

[x_quad2D, w2D] = GaussQuad2D(4);
phi = GetPhi2D(x_quad2D);%d/d_xi in 1st col, x/d_eta in 2nd col, diff phi per 3D card

F = zeros(n_total_dofs,1);

[kb sb] = find(ElemNeighbors == -1);%the element on the boundary kb and the side of that element sb
BExt=[kb sb];% 2D matrix with [boundary element, face id] in each row


x_hat = RefNodeLocations;
for elem_id = 1:size(EToN,1)
phys_pts = GetPhysicalPoints2D(elem_id, x_hat);
dofs_e = DofMap{1}(elem_id,:);
F(dofs_e) = f(phys_pts(1,:)',phys_pts(2,:)');
end


C = norm(advection_velocity);

RHS = ex*(-S{1} + Cavg{1}+BC{1})*inv(M)*(Cavg{1}-S{1}) + ey*(-S{2} + Cavg{2}+BC{2})*inv(M)*(Cavg{2}-S{2}) - epsilon*tau*JumpCouplingMatrix... %DIFFUSION
      + (S{1} - Cavg{1})*advection_velocity(1) + (S{2} - Cavg{2})*advection_velocity(2) + C/2*JumpCouplingMatrix; %ADVECTION


if Fdirac
    F(:) = 0;
    N = sqrt(sum(abs(VX).^2,2));
    dirac_DOF = find(N<1e-2);
    F(dirac_DOF) = 1;
end


% Solve for U
U = RHS\(M*F);
solution_coeffs = U;
figure(2)
plot_q
title(sprintf('Steady State Advection Diffusion for f = %s\n u = [%g %g] p=%g h = %g epsilon = %g tau = %g',char(f),advection_velocity(1),advection_velocity(2),p,h,epsilon,tau))
colorbar




%{
Outline small issue with standard mesh from distmesh

the mseh aligns itself with vertical lines predominantly running
horizontally

for purely horizontal diffusion, no tau is req for stability with sin(x)
for purely vertical diffusion with sin(y) we need a tau for stability

similarly
after flipping hte mesh about the x=y line with 
VX = [VX(:,2) VX(:,1)]
for purely horizontal diffusion, a tau is req for stability with sin(x)
ex = 1, ey = 0

for purely vertical diffusion with sin(y) no tau is required for stability
ex = 0 ey = 1

%}











