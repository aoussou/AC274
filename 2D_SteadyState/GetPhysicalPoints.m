function phys_pts = GetPhysicalPoints(elem_id, reference_pts)
% Map reference_pts to physical space

Globals1D;

% get the node IDs of the nodes that belong to elem
node_ids = EToN(elem_id,:);

% get the corresponding node locations
x1 = VX(node_ids(1));
x2 = VX(node_ids(2));
x2_minus_x1 = x2 - x1;
if periodic_BCs_flag
    x2_minus_x1 = x2_minus_x1 - domain_length.*(x2_minus_x1>domain_length/2) + domain_length.*(x2_minus_x1<-domain_length/2);
end

% the mapping from [-1,1] to [x1,x2] is
phys_pts = x1 + 0.5*(1+reference_pts)*x2_minus_x1;

