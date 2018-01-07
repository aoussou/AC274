function f_boundary = ApplyBoundaryFluxFunction(q,time)
% compute the flux for various problems

Globals1D;

% reshape q so we can apply the flux function column by column
calN = K*(p+1);
q_reordered = zeros(n_vars,calN);

for var = 1:n_vars
    first_index = calN*(var-1) + 1;
    last_index = var*calN;
    q_reordered(var,:) = q(first_index:last_index)';
end

% initialize flux to 0
n_total_dofs = n_vars*calN;
f_boundary = zeros(n_total_dofs,1);
f_reordered = 0*q_reordered;


% calculate boudary flux for various problems
switch bndry_condition_string
    
    case 'None'
        % do nothing
        
    case 'Outflow2'
        
        % since we have a BC on q1 on both boundaries, we need to set
        % the boundary values of q1 to zero here
        q_reordered(1,1) = 0;
        
        % the right-most node index depends on p
        q_reordered(1,end) = 0;
        
        % the wave equation matrix
        A = [0 wave_speed^2; 1 0];
        
        f_reordered = A*q_reordered;
        
    case 'Problem2a'
        
        % since we have a BC on q1 on both boundaries, we need to set
        % the boundary values of q1 to zero here
        q_reordered(1,1) = 0;
        
        q_reordered(2,end) = 0;
        
        % the wave equation matrix
        A = [0 wave_speed^2; 1 0];
        
        f_reordered = A*q_reordered;
        
    case 'Problem2b'
        
        % since we have a BC on q1 on both boundaries, we need to set
        % the boundary values of q1 to zero here
        q_reordered(1,1) = 0;
        
        q_reordered(1,end) = 0;
        
        % the wave equation matrix
        A = [0 wave_speed^2; 1 0];
        
        f_reordered = A*q_reordered;
        
        
    case 'Problem2c'
        
        % since we have a BC on q1 on both boundaries, we need to set
        % the boundary values of q1 to zero here
        a=0;b=2;
        q_reordered(1,1) = -0.6*wave_speed*pi*cos(0.6*pi*(a+0.5))*sin(0.6*wave_speed*pi*time);

        q_reordered(1,end) = -0.6*wave_speed*pi*cos(0.6*pi*(b+0.5))*sin(0.6*wave_speed*pi*time);
        
        % the wave equation matrix
        A = [0 wave_speed^2; 1 0];
        
        f_reordered = A*q_reordered;
        
    case 'Problem3a'
        q_reordered(2,1) = 0;%no flux BC imposed by zeroing out the velocity at a = -5 and b = 5
        q_reordered(2,end) =0;
        
        g = 1;%normalized gravitational constant
        
        f_reordered = [q_reordered(2,:);q_reordered(2,:).^2./q_reordered(1,:) + 0.5*g*q_reordered(1,:).^2];
        
        if problem_3a_flag
            C = 3;
        end
        
    case 'Problem3b'
        q_reordered(2,1) = 0;%no flux BC imposed by zeroing out the velocity at a = -5 and b = 5
        q_reordered(2,end) =0;
        
        g = 1;%normalized gravitational constant
        
        f_reordered = [q_reordered(2,:);q_reordered(2,:).^2./q_reordered(1,:) + 0.5*g*q_reordered(1,:).^2];
        
        if problem_3b_flag
            C = 3;
        end
        
end


% finally put f into the standard ordering
n_total_dofs = n_vars*calN;
f_boundary = zeros(n_total_dofs,1);
for var = 1:n_vars
    first_index = calN*(var-1) + 1;
    last_index = var*calN;
    f_boundary(first_index:last_index) = f_reordered(var,:)';
end