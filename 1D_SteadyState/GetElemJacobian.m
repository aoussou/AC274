function ElemJacobian = GetElemJacobian(elem_id)
% Get the Jacobian of the mapping from reference to physical element

% Globals1D;
global EToN
global VX
% get the node IDs of the nodes that belong to elem
node_ids = EToN(elem_id,:);

% get the corresponding node locations
x1 = VX(node_ids(1));
x2 = VX(node_ids(2));
x2_minus_x1 = x2 - x1;


% the mapping from [-1,1] to [x1,x2] is
% x = x1 + 0.5*(1+x_hat)*(x2-x1)
ElemJacobian = 0.5*x2_minus_x1;

