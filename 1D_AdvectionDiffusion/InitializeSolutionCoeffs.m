function solution_coeffs = InitializeSolutionCoeffs(initial_condition_string)
% initialize the solution coefficients by projecting
% the specified function into the finite element space

Globals1D;

switch initial_condition_string
    
    case 'Gaussian'
        q0 = @(x)(exp(-100*(x).^2));
%         q0 = @(x)exp(-pi*(x-1).^2);
%         q0 = @(x)(exp(-.1*(x).^2));%area under here from -1 to 1 should be 2
%           q0 = @(x)exp(-200*(1-x).^2);

        
   case 'SquarePulse'
                q0 = @(x)heaviside(x+.5).*heaviside(0.5-x);

    otherwise
       error('Invalid initial_condition_string')   
        
end

% Get the node locations on the reference element
ref_nodes = RefNodeLocations(p);

% retrieve the basis functions at the quadrature points
% on the reference element
phi = GetPhi(p, ref_nodes);

% allocate space for the sparse mass matrix
% the matrix will have a (p+1)x(p+1) block
% for each element
n_local_dofs = p+1;
n_total_dofs = K*n_local_dofs;
solution_coeffs = zeros(n_total_dofs, 1);
for elem_id=1:K
    local_dofs = DofMap(elem_id,:);
    
    % get the image of the node locations points in physical space
    phys_pts = GetPhysicalPoints(elem_id, ref_nodes);

    % evaluate q0 at phys_points in order to interpolate
    solution_coeffs(local_dofs) = q0(phys_pts);
end

return


% helper functions for defining initial conditions

function q0 = square_wave(x)

% set to 1 on [0.2+eps,0.4-eps] to make sure we get the right looking
% initial condition
epsilon = 0;

q0 = 0*x;
indices = (0.2 <= x) & (x <= 0.4);
q0(indices) = 1;
return
