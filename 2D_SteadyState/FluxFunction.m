function f = FluxFunction(q)
% Define the flux function for this problem

% Globals2D;
global flux_function_string
global advection_velocity

switch flux_function_string
    
    case 'advection'
        f = [advection_velocity(1)*q,advection_velocity(2)*q];

    otherwise
        error('global variable flux_function_string is not specified correctly')
        
end