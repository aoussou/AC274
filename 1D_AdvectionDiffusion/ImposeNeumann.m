
function [Q] = ImposeNeumann(Q)

global bndry_condition_string

switch bndry_condition_string
    
    case 'NoPenetration'
%         Q(1) = Q(2);
%         Q(end) = Q(end-1);
        Q(1) = 0;
        Q(end) = 0;

    otherwise
        %do nothing
        
        
end

