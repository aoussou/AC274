function f_boundary = BoundaryFluxFunction(solution_coeffs)

global bndry_condition_string

switch bndry_condition_string
    case 'Outflow'
        f_boundary = 0*solution_coeffs;
    otherwise
        %do nothing
        f_boundary = 0*solution_coeffs;
end

        