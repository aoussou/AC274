function dphi = GetDPhi(p, x_hat)
% Build an array storing the derivative of the basis functions at the
% quadrature points x_quad

NumPts = length(x_hat);
dphi = zeros(NumPts,p+1);

switch p
    
    case 0
        % don't need to do anything, dphi is already zero
        
    case 1
        dphi(:,1) = -0.5;
        dphi(:,2) = 0.5;
        
    case 2
        dphi(:,1) = x_hat - 0.5*ones(size(x_hat));
        dphi(:,2) = -2*x_hat;
        dphi(:,3) = x_hat + 0.5*ones(size(x_hat));
        
    case 3
        dphi(:,1) = -9/16*( (x_hat+1/3).*(x_hat-1/3) + (x_hat+1/3).*(x_hat-1) + (x_hat-1/3).*(x_hat-1) );
        dphi(:,2) = 27/16*( (x_hat+1).*(x_hat-1/3) + (x_hat+1).*(x_hat-1) + (x_hat-1/3).*(x_hat-1) );
        dphi(:,3) = -27/16*( (x_hat+1).*(x_hat+1/3) + (x_hat+1).*(x_hat-1) + (x_hat+1/3).*(x_hat-1) );  
        dphi(:,4) = 9/16*( (x_hat+1).*(x_hat+1/3) + (x_hat+1).*(x_hat-1/3) + (x_hat+1/3).*(x_hat-1/3) ); 

    case 4
        dphi(:,1) =  2/3*( (x_hat+.5).*(x_hat-0).*(x_hat-.5) + (x_hat+.5).*(x_hat-0).*(x_hat-1) + (x_hat+.5).*(x_hat-.5).*(x_hat-1) + (x_hat-0).*(x_hat-.5).*(x_hat-1) );
        dphi(:,2) = -8/3*( (x_hat+1).*(x_hat-0).*(x_hat-.5) + (x_hat+1).*(x_hat-0).*(x_hat-1) + (x_hat+1).*(x_hat-.5).*(x_hat-1) + (x_hat-0).*(x_hat-.5).*(x_hat-1) );
        dphi(:,3) =    4*( (x_hat+1).*(x_hat+.5).*(x_hat-.5) + (x_hat+1).*(x_hat+.5).*(x_hat-1) + (x_hat+1).*(x_hat-.5).*(x_hat-1) + (x_hat+.5).*(x_hat-.5).*(x_hat-1) );  
        dphi(:,4) = -8/3*( (x_hat+1).*(x_hat+.5).*(x_hat-0) + (x_hat+1).*(x_hat+.5).*(x_hat-1) + (x_hat+1).*(x_hat-0).*(x_hat-1) + (x_hat+.5).*(x_hat-0).*(x_hat-1) );
        dphi(:,5) =  2/3*( (x_hat+1).*(x_hat+.5).*(x_hat-0) + (x_hat+1).*(x_hat+.5).*(x_hat-.5) + (x_hat+1).*(x_hat-0).*(x_hat-.5) + (x_hat+.5).*(x_hat-0).*(x_hat-.5) );        
        
    otherwise
       error('Invalid polynomial order in GetDPhi')
        
end
