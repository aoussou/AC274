function f = FluxFunction(u)
% Define the flux function for this problem

% Globals1D;

global advection_velocity
global flux_function_string
switch flux_function_string
    
    case 'advection'
        f = advection_velocity*u;
    case 'diffusion'
        f = u;
    otherwise
        error('global variable flux_function_string is not specified correctly')
        
end