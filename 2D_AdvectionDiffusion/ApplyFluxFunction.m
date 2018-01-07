function f = ApplyFluxFunction(q)
% compute the flux for this problem

Globals2D;

% reshape q so we can apply the flux function column by column

        if p==0
            n_local_dofs = 1;
        else
            n_local_dofs = 3*p;
        end
        calN = K*n_local_dofs;
n_total_dofs = n_vars*calN;

q_reordered = cell(n_vars,1);

for i = 1:n_vars
%     q_reordered(i,:) = q((i-1)*calN+1:i*calN)';
    q_reordered{i} = q((i-1)*calN+1:i*calN);
end

% q{i} is a column vector???

f_reordered = FluxFunction(q_reordered);%for advection will output [f_x1 , f_x_2]

% finally put f into the standard ordering


%     f((i-1)*calN+1:i*calN) = f_reordered(i,:)';

switch bndry_condition_string
    case 'Problem1a'
        f = zeros(n_total_dofs,n_dim);
        i=1;
        f((i-1)*calN+1:i*calN,1:n_dim) = f_reordered{i};
    case 'Problem1d'
%         f = zeros(n_total_dofs,n_dim);
%         i=1;
%         f((i-1)*calN+1:i*calN,1:n_dim) = f_reordered{i};
f = f_reordered;
    case 'Problem2'
        for i = 1:n_vars
            f((i-1)*calN+1:i*calN,1:2) = f_reordered{i};
        end
    case 'Problem2b'
        for i = 1:n_vars
            f((i-1)*calN+1:i*calN,1:2) = f_reordered{i};
        end
    case 'Problem3'
        for i = 1:n_vars
            f((i-1)*calN+1:i*calN,1:2) = f_reordered{i};
        end
otherwise
        for i = 1:n_vars
            f((i-1)*calN+1:i*calN,1:2) = f_reordered{i};
        end
end



