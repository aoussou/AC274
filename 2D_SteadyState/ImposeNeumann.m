
function [Q] = ImposeNeumann(Q)
% Globals2D;
global bndry_condition_string

switch bndry_condition_string
    
    case 'NoPenetration'

        global dofs_l dofs_b dofs_t dofs_r


        Q(dofs_l,1) =0;
        Q(dofs_r,1) =0;
        Q(dofs_b,2) =0;
        Q(dofs_t,2) =0;


    otherwise
        %do nothing

        
end

