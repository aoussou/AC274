%% plot the solution coefficients
%plot_q

figure(2)

%% compute the values of the interpolation function
if n_dim == 2
xplotref = [0 0.5 1 0.5 0 0]';
yplotref = [0 0 0 0.5 1 0.5]';
end

Pplotref = [xplotref, yplotref];
phi = GetPhi2D(Pplotref);



for elem_id = 1:size(EToN,1)

hold on
    
    phys_points = GetPhysicalPoints2D(elem_id, Pplotref);%takes the points in
    %             phys_pts_check((k-1)*(numQpts)+1:k*(numQpts))= phys_points;
    ElemJacobian = GetElemJacobian2D(elem_id);
    xpoly = phys_points(1,:);
    ypoly = phys_points(2,:);
    switch p
        case 0
            mappedValues = [solution_coeffs(DofMap{1}(elem_id)) solution_coeffs(DofMap{1}(elem_id)) solution_coeffs(DofMap{1}(elem_id))...
                   solution_coeffs(DofMap{1}(elem_id)) solution_coeffs(DofMap{1}(elem_id)) solution_coeffs(DofMap{1}(elem_id))];
        case 1
            
            mappedValues = [solution_coeffs(DofMap{1}(elem_id,1)) solution_coeffs(DofMap{1}(elem_id,2)) solution_coeffs(DofMap{1}(elem_id,2))...
                   solution_coeffs(DofMap{1}(elem_id,3)) solution_coeffs(DofMap{1}(elem_id,3)) solution_coeffs(DofMap{1}(elem_id,1))];
        case 2
    mappedValues =[solution_coeffs(DofMap{1}(elem_id,1)) solution_coeffs(DofMap{1}(elem_id,4)) solution_coeffs(DofMap{1}(elem_id,2))...
                   solution_coeffs(DofMap{1}(elem_id,5)) solution_coeffs(DofMap{1}(elem_id,3)) solution_coeffs(DofMap{1}(elem_id,6))];
    end
    
    patch(xpoly,ypoly,mappedValues');
    
    
    hold on

    
end
        
hold off
axis equal
colorbar
switch bndry_condition_string
case 'Problem2'
caxis([-.02 0.02])
end