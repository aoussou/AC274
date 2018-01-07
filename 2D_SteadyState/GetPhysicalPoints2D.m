function phys_pts = GetPhysicalPoints2D(elem_id, reference_pts)
% Map reference_pts to physical space

global EToN
global VX

v1 = VX(EToN(elem_id,1),:);
v2 = VX(EToN(elem_id,2),:);
v3 = VX(EToN(elem_id,3),:);

A = [-v1(1)+v2(1) -v1(1)+v3(1);...
     -v1(2)+v2(2) -v1(2)+v3(2)];
B = [v1(1);v1(2)];


phys_pts = A*reference_pts'+ B*ones(1,size(reference_pts,1));
