function phi = GetPhi2D(xi_eta)
% Build an array storing the basis functions at the
% points on the reference element, reference_pts
% (matlab implementation seems fastest if I write out the basis explicitly)

%xi_eta is the unit element coordinate matrix

Globals2D;

xi = xi_eta(:,1);
eta = xi_eta(:,2);

npts = size(xi_eta,1);



% [XI,ETA] = meshgrid(xi,eta);

% phi0 = zeros(npts,npts);

switch p
    
    case 0
        phi(:,1) = ones(size(xi,1),1);
    case 1
        
        phi(:,1) = 1 - xi - eta;
        phi(:,2) = xi;
        phi(:,3) = eta;
        
    case 2
        phi = zeros(npts,6);
        phi(:,1) = 2*(0.5 - xi - eta).*(1 - xi - eta);
        phi(:,2) = 2*xi.*(xi - 0.5);
        phi(:,3) = 2*eta.*(eta - 0.5);
        phi(:,4) = 4*xi.*(1 - xi - eta);
        phi(:,5) = 4*xi.*eta;
        phi(:,6) = 4*eta.*(1 - xi - eta);
        
end