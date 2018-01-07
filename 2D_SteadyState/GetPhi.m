function phi = GetPhi( x_hat)
% Build an array storing the basis functions at the
% points on the reference element, reference_pts
Globals2D;

NumPts = length(x_hat);
phi = zeros(NumPts,p+1);

switch p
    
    case 0
        phi(:,1) = ones(size(x_hat));
        
    case 1
        phi(:,1) = 0.5*(1-x_hat);
        phi(:,2) = 0.5*(1+x_hat);
        
    case 2
        phi(:,1) = 0.5*x_hat.*(x_hat-1);
        phi(:,2) = (1-x_hat).*(1+x_hat);
        phi(:,3) = 0.5*x_hat.*(1+x_hat);

    case 3
        phi(:,1) = -9/16 *(x_hat+1/3).*(x_hat-1/3).*(x_hat-1);
        phi(:,2) = 27/16 *(x_hat+1).*(x_hat-1/3).*(x_hat-1);
        phi(:,3) = -27/16*(x_hat+1).*(x_hat+1/3).*(x_hat-1);  
        phi(:,4) = 9/16  *(x_hat+1).*(x_hat+1/3).*(x_hat-1/3); 

    case 4
        phi(:,1) =  2/3*(x_hat+.5).*(x_hat-0).*(x_hat-.5).*(x_hat-1);
        phi(:,2) = -8/3*(x_hat+1) .*(x_hat-0).*(x_hat-.5).*(x_hat-1);
        phi(:,3) =    4*(x_hat+1).*(x_hat+.5).*(x_hat-.5).*(x_hat-1);  
        phi(:,4) = -8/3*(x_hat+1).*(x_hat+.5).*(x_hat-0) .*(x_hat-1);
        phi(:,5) =  2/3*(x_hat+1).*(x_hat+.5).*(x_hat-0) .*(x_hat-.5);
        
    otherwise
       error('Invalid polynomial order in GetPhi')     
        
end
