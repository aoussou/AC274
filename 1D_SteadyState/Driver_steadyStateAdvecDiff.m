

%% Driver  for steady state advection-diffusion 

clear all; clc;close all
Globals1D;

%% set problem parameters

a                        = 0;        % mesh min x
b                        = pi*2;     % mesh max x
K                        = 100;       % number of elements
advection_velocity       = 0;
p                        = 4;        % basis function polynomial order

tau                      = 1;      % Penalty term coefficient ( >0)
epsilon                  = 1;        % Diffusion coefficient
BL                       = 0;   %Dirichlet Boundary Conditions left node
BR                       = 0;   %Dirichlet Boundary Conditions right node
f                        = @(x)-sin(x); % Forcing Function / Source Term (fxn of x)


dx = (b-a)/K;
CFL = 2.5*delta_t/dx;
if CFL > .9999
    disp('Warning: CFL VIOLATED')
    return
end

%% initialization

% generate a mesh with K equally sized elements
[VX, EToN] = MeshGen1D(a, b, K);

% generate the Neighbor array that defines element connectivity
ElemNeighbors = GetConnectivity(EToN);

% generate the degree-of-freedom map, DofMap
CreateDofMap;

% assemble our discretization matrices
AssembleMassMatrix;
AssembleAdvectionMatrix;
AssembleAvgCouplingMatrixFinal;
AssembleJumpCouplingMatrix;
AssembleBoundaryConditionMatrix;
Cj =JumpCouplingMatrix;

n_local_dofs = p+1;
calN = n_local_dofs*K;

nGQ_pts = 20;

[x_hat,w] = GaussQuad(nGQ_pts);
w = w';
phi = GetPhi(p, x_hat);
F = zeros(calN,1);
LineElemJacobian = GetElemJacobian(1);% all elements are the same length

Fdirac = F;
Fdirac(floor(calN/2))=1;



x_hat = RefNodeLocations(p);
for elem_id = 1:K
    phys_pts = GetPhysicalPoints(elem_id, x_hat);
    dofs_e = DofMap(elem_id,:);
    
    F(dofs_e) = f(phys_pts);
end


BCin = zeros(size(BC,1),1);
BCin(1) = BL;
BCin(end) = -BR;

RHS = ((epsilon*(-S + Cavg+BC)*inv(M)*(Cavg-S) + (S - Cavg)*advection_velocity - epsilon*tau*JumpCouplingMatrix));


U = RHS\(M*F+epsilon*(-S + Cavg + BC)*inv(M)*BCin);

U_reshape = reshape(U,p+1,K);

plot(linspace(a,b,(p+1)*K),U,'b','LineWidth',1)
title(sprintf('f(x) = -sin(x), tau = %g  u(%g) = %g u(2*pi) = %g',tau,a,BL,BR))
hold on








