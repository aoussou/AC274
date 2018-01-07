function [solution_coeffs] = ImposeDirichlet(solution_coeffs)

global bndry_condition_string
switch bndry_condition_string
    case 'Outflow'
        solution_coeffs(1) = solution_coeffs(2);
        solution_coeffs(end) = solution_coeffs(end-1);

    otherwise
        %do nothing
        
        
end