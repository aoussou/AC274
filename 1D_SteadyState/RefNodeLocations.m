function ref_nodes = RefNodeLocations(p)
% Build an array storing the locations of the nodes
% on the reference element

ref_nodes = zeros(p+1,1);

switch p
    
    case 0
        ref_nodes(1) = 0;
        
    case 1
        ref_nodes(1) = -1;
        ref_nodes(2) = 1;
        
    case 2
        ref_nodes(1) = -1;
        ref_nodes(2) = 0;
        ref_nodes(3) = 1;

    case 3
        ref_nodes(1) = -1;
        ref_nodes(2) = -1/3;
        ref_nodes(3) = 1/3;  
        ref_nodes(4) = 1; 

    case 4
        ref_nodes(1) = -1;
        ref_nodes(2) = -0.5;
        ref_nodes(3) = 0;  
        ref_nodes(4) = 0.5; 
        ref_nodes(5) = 1;
        
    otherwise
       error('Invalid polynomial order in RefNodeLocations')     
        
end
