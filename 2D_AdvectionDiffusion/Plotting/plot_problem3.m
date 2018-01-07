
%% plot the solution coefficients
%plot_problem3

addpath('..');
load('solutionCoeffs.mat')

titles = cell(4,1);
titles{1} = 'rho';
titles{2} = 'rho * u';
titles{3} = 'rho * v';
titles{4} = 'E';


%% compute the values of the interpolation function
if p==2
xplotref = [0 0.5 1 0.5 0 0]';
yplotref = [0 0 0 0.5 1 0.5]';
elseif p ==0
    xplotref = [0 1 0 ]';
    yplotref = [0 0 1 ]';
end

Pplotref = [xplotref, yplotref];
phi = GetPhi2D(Pplotref);


for i = 2
    
figure(i);
for elem_id = 1:size(EToN,1)

hold on
    
    phys_points = GetPhysicalPoints2D(elem_id, Pplotref);%takes the points in
    %             phys_pts_check((k-1)*(numQpts)+1:k*(numQpts))= phys_points;
    ElemJacobian = GetElemJacobian2D(elem_id);
    xpoly = phys_points(1,:);
    ypoly = phys_points(2,:);
    
    switch p
        case 0
            mappedValues = [solution_coeffs(DofMap{i}(elem_id)) solution_coeffs(DofMap{i}(elem_id)) solution_coeffs(DofMap{i}(elem_id))];
        case 1
            mappedValues = [solution_coeffs(DofMap{1}(elem_id,1)) solution_coeffs(DofMap{1}(elem_id,2)) solution_coeffs(DofMap{1}(elem_id,2))...
                   solution_coeffs(DofMap{1}(elem_id,3)) solution_coeffs(DofMap{1}(elem_id,3)) solution_coeffs(DofMap{1}(elem_id,1))];
        case 2
    mappedValues =[solution_coeffs(DofMap{1}(elem_id,1)) solution_coeffs(DofMap{1}(elem_id,4)) solution_coeffs(DofMap{1}(elem_id,2))...
                   solution_coeffs(DofMap{1}(elem_id,5)) solution_coeffs(DofMap{1}(elem_id,3)) solution_coeffs(DofMap{1}(elem_id,6))];
    end
    patch(xpoly,ypoly,mappedValues');

end
hold off
title(sprintf('%s  at t = 4',titles{i}))
end
        



% 
% 
% if plotMach
% 
% 
% %% compute the values of the interpolation function
% if n_dim == 2
% xplotref = [0 0.5 1 0.5 0 0]';
% yplotref = [0 0 0 0.5 1 0.5]';
% end
% 
% Pplotref = [xplotref, yplotref];
% phi = GetPhi2D(Pplotref);
% 
% figure(2)
% 
% for elem_id = 1:K
% 
% hold on
%     
%     phys_points = GetPhysicalPoints2D(elem_id, Pplotref);%takes the points in
%     %             phys_pts_check((k-1)*(numQpts)+1:k*(numQpts))= phys_points;
%     ElemJacobian = GetElemJacobian2D(elem_id);
%     xpoly = phys_points(1,:);
%     ypoly = phys_points(2,:);
%     
% for i = 1:4
% subplot(2,2,i)
%     mappedValues =[Mach(DofMap{i}(elem_id,1)) Mach(DofMap{i}(elem_id,4)) Mach(DofMap{i}(elem_id,2))...
%                    Mach(DofMap{i}(elem_id,5)) Mach(DofMap{i}(elem_id,3)) Mach(DofMap{i}(elem_id,6))];
% end
%     patch(xpoly,ypoly,mappedValues');
%     
%     hold on
% 
% end
%         
% hold off
% axis equal
% colorbar
% 
% end
% 
% 
% 



% caxis([-.02 0.02])