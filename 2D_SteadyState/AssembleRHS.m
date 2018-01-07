function RHS = AssembleRHS(solution_coeffs,time)
% Assemble the right-hand side vector

Globals2D;

global M
global S
global Cavg
global BC
global JumpCouplingMatrix
global advection_velocity
global epsilon
global bndry_condition_string
global dofs_l dofs_b dofs_t dofs_r

%------------ RHS Assembly Procedure --------------------------------------

%{
1. Set the solution coefficient vector [n_tot_dofs,1] to the q's in cells
so that q{mew} = solution_coeffs( dofs corresponding to the mew-th variable

2. Apply the flux function to the q's to create the F vector

     ***The FluxFunction accepts a {n_vars,1} cell array q, where each element of the
        cell is a [calN,1] vector
        The OUTPUT of the Flux Function is a {n_vars,1} cell array f where each element
        of the f cell is a [calN,2] vector for [x1,x2]

3. Enforce the boundary conditions weakly through the flux by enforcing the
relevant values of q{i} to create a boundary cell array qb{i} for the ith
variable

4. Apply the flux function to the qb's to create the Fb vector

5. Assemble the RHS

%}

if p==0
    n_local_dofs = 1;
else
    n_local_dofs = 3*p;
end
calN = n_local_dofs*size(EToN,1);



C = norm(advection_velocity);

f = FluxFunction(solution_coeffs);
fb = f;
% fb(dofs_N,:) = [solution_coeffs(dofs_N,:),solution_coeffs(dofs_N,:)];
fb(dofs_N,:) = 0;

RHS1_1 = (Cavg{1}-S{1})*solution_coeffs;

RHS1_2 = (Cavg{2}-S{2})*solution_coeffs;

Q(:,1) = M\RHS1_1;
Q(:,2) = M\RHS1_2;

% Impose Neumann
Q = ImposeNeumann(Q);

RHS = epsilon*(-S{1} + Cavg{1})*Q(:,1) + epsilon*(-S{2} + Cavg{2})*Q(:,2)+...
       (S{1}  - Cavg{1})*f(:,1) + (S{2}  - Cavg{2})*f(:,2) - ...
         C/2*JumpCouplingMatrix*solution_coeffs-BC{1}*fb(:,1)-BC{2}*fb(:,2);






