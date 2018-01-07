function [VX, EToN] = MeshGen1D(a,b,K)
% Generate uniform mesh with K elements in the interval [a,b]
% VX:   The x-positions of the nodes
% EToN: K x 2 array, where EToN(k,:) gives the
%       nodes associated with the k^th element

% Generate node coordinates
N_nodes = K+1;
VX = zeros(1,N_nodes);
for i = 1:N_nodes
  VX(i) = (b-a)*(i-1)/K + a;
end

% generate element to node connectivity
EToN = zeros(K, 2);
for k = 1:K
  EToN(k,1) = k; EToN(k,2) = k+1;
end
