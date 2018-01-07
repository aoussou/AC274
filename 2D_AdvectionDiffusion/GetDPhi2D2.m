function [dphidxi, dphideta] = GetDPhi2D2(xi_eta)
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
        
    case 1
        
% phi(:,1)   = 1 - xi - eta;
% phi(:,:,2) = xi;
% phi(:,:,3) = eta;
        
% partial wrt xi in first column and wrt eta in second column

dphidxi(:,1) = -1;
dphideta(:,1) = -1;

dphidxi(:,2) = 1;
dphideta(:,2) = 0;

dphidxi(:,3) = 0;
dphideta(:,3) = 1;

    case 2

%         phi1 = 2*(0.5 - xi - eta).*(1 - xi - eta);
%         phi2 = 2*xi.*(xi - 0.5);
%         phi3 = 2*eta.*(eta - 0.5);
%         phi4 = 4*xi.*(1 - xi - eta);
%         phi5 = 4*xi.*eta;
%         phi6 = 4*eta.*(1 - xi - eta);
        
        dphidxi(:,1) = 4*eta + 4*xi - 3;
        dphideta(:,1) = 4*eta + 4*xi - 3;

        dphidxi(:,2) = 4*xi - 1;
        dphideta(:,2) = zeros(length(eta),1);

        dphidxi(:,3) = zeros(length(eta),1);
        dphideta(:,3) = 4*eta - 1;
        
        dphidxi(:,4) = 4 - 8*xi - 4*eta;
        dphideta(:,4) = -4*xi;
        
        dphidxi(:,5) = 4*eta;
        dphideta(:,5) = 4*xi;
        
        dphidxi(:,6) = -4*eta;
        dphideta(:,6) = 4 - 4*xi - 8*eta;
        
        
        
% phi = [phi1; phi2; phi3; phi4; phi5; phi6];        
 end 