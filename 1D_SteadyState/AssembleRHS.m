function RHS = AssembleRHS(solution_coeffs,time)
% Assemble the right-hand side vector

% Globals1D;
global M
global S
global Cavg
global BC
global JumpCouplingMatrix
global advection_velocity
global epsilon
global bndry_condition_string


C = abs(advection_velocity);

f = FluxFunction(solution_coeffs);

% boundary flux function
fb = BoundaryFluxFunction(solution_coeffs);

RHS1 = (Cavg-S)*solution_coeffs;

Q = M\RHS1;

%Impose Neumann
Q = ImposeNeumann(Q);


RHS = epsilon*(-S + Cavg)*Q + (S  - Cavg)*f - ...
                C/2*JumpCouplingMatrix*solution_coeffs+BC*fb;
  
