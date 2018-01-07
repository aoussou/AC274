function ElemJacobian = GetElemJacobian2D(elem_id)
% Get the Jacobian of the mapping from reference to physical element

global EToN
global VX

v1 = VX(EToN(elem_id,1),:);
v2 = VX(EToN(elem_id,2),:);
v3 = VX(EToN(elem_id,3),:);


ElemJacobian = [-v1(1)+v2(1) -v1(1)+v3(1);...
    -v1(2)+v2(2) -v1(2)+v3(2)];
