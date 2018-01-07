function dphi = GetDPhi2D(xi_eta)
% Build an array storing the basis functions at the
% points on the reference element, reference_pts
% (matlab implementation seems fastest if I write out the basis explicitly)

%xi_eta is the unit element coordinate matrix

Globals2D;

xi = xi_eta(:,1);
eta = xi_eta(:,2);

% [xi,eta] = meshgrid(xi,eta);

% phi0 = zeros(npts,npts);

switch p
    
    case 0
        dphi = [0 0];
    case 1
        
        % phi(:,1)   = 1 - xi - eta;
        % phi(:,:,2) = xi;
        % phi(:,:,3) = eta;
        
        % partial wrt xi in first column and wrt eta in second column, phi
        % on diff cards
        
        dphi(:,:,1) = -1*ones(length(xi),2);
        dphi(:,:,2) = [ones(length(xi),1)  zeros(length(xi),1) ];
        dphi(:,:,3) = [zeros(length(xi),1) ones(length(xi),1)  ];
        
        
    case 2
        
        
        dphi(:,1,1) = 4*eta + 4*xi - 3;
        dphi(:,2,1) = 4*eta + 4*xi - 3;
        
        dphi(:,1,2) = 4*xi - 1;
        dphi(:,2,2) = zeros(length(eta),1);
        
        dphi(:,1,3) = zeros(length(eta),1);
        dphi(:,2,3) = 4*eta - 1;
        
        dphi(:,1,4) = 4 - 8*xi - 4*eta;
        dphi(:,2,4) = -4*xi;
        
        dphi(:,1,5) = 4*eta;
        dphi(:,2,5) = 4*xi;
        
        dphi(:,1,6) = -4*eta;
        dphi(:,2,6) = 4 - 4*xi - 8*eta;
        
        
        
        % phi = [phi1; phi2; phi3; phi4; phi5; phi6];
end