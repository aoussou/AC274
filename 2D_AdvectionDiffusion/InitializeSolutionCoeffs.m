function solution_coeffs = InitializeSolutionCoeffs
% initialize the solution coefficients by projecting
% the specified function into the finite element space

Globals2D;

switch initial_condition_string
    
%     case 'Sinusoid2'
%         L = b-a;
%         k = 9;
%         init_qi = {@(x) 0*x , ...
%             @(x) -pi*k/L*cos(pi*k*x/L)};
%         
%     case 'Gaussian2'
%         % zero velocity
%         % u0(x) = exp(-100*(x-1).^2)
%         init_qi = {@(x) 0*x, ...
%             @(x) 200*exp(-100*(x-1).^2).*(x-1) };
    case 'Gaussian'
        % zero velocity
        % u0(x) = exp(-100*(x-1).^2)
        init_qi = {@(x1,x2) 100*exp(-100*((x1).^2+(x2).^2))};
%     case 'Problem2'
%         % zero velocity
%         % u0(x) = exp(-100*(x-1).^2)
%         init_qi = {@(x1,x2) 0*x1;
%                    @(x1,x2) 0*x1;
%                    @(x1,x2) 0*x1};
%     case 'Problem2b'
%         
%         % q = [ h ; hu ; hv ] 
%         
%         init_qi = {@(x1,x2) ones(length(x1),1);
%                    @(x1,x2) 0*x1;
%                    @(x1,x2) 0*x1};
%     case 'Problem3'
%         
%         %q = [rho rho*u rho*v E+Pressure];
%         
%         init_qi = {@(x1,x2) ones(length(x1),1);
%                    @(x1,x2) ones(length(x1),1);
%                    @(x1,x2) zeros(length(x1),1);
%                    @(x1,x2) ones(length(x1),1)};
%                
%     case 'testProblem2a'
% 
%         init_qi = {@(x1,x2) exp(-100*(x1.^2+x2.^2));
%                    @(x1,x2) zeros(length(x1),1);
%                    @(x1,x2) zeros(length(x1),1)};

    otherwise
        error('Invalid initial_condition_string')
        
end

% Get the node locations on the reference element
ref_nodes = RefNodeLocations;

% retrieve the basis functions at the quadrature points
% on the reference element

        phi = GetPhi2D(ref_nodes);

% calculate initial solution coefficients

if p==0
    n_local_dofs = 1;
else
    n_local_dofs = 3*p;
end
calN = n_local_dofs*size(EToN,1);
n_total_dofs = n_vars*calN;
ref_nodes = RefNodeLocations;

solution_coeffs = zeros(n_total_dofs, 1);

for elem_id=1:size(EToN,1)
     phys_pts = GetPhysicalPoints2D(elem_id, ref_nodes);  % get the image of the node locations points in physical space
                                                          %  ref_nodes = 2xnum_points 
                                                          % OUTPUT:physical points 2xN_points s.t.[x;y]
    % evaluate qi at phys_points in order to interpolate
    for i = 1:n_vars
        local_dofs_qi = DofMap{i}(elem_id,:);
        solution_coeffs(local_dofs_qi) = init_qi{i}(phys_pts(1,:),phys_pts(2,:));
    end
end





% helper functions for defining initial conditions
