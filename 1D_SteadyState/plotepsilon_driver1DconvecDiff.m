PLot epsilon

if wpsilon == 1
%     c= 'b';
%     if norm(double(time==[1.5]))
%         h = figure(2);
%         reference_pts = linspace(-1,1,10)';
%         n_reference_pts = length(reference_pts);
%         phi = GetPhi(p, reference_pts);
%         
%         n_local_dofs = p+1;
%         
%         %account for ghost elements
%         if ghost_cell_flag
%             elem1 = 2;     elemK = K-1;
%         else
%             elem1 = 1;     elemK = K;
%         end
%         
%         for elem_id=elem1:elemK
%             local_dofs = DofMap(elem_id,:);
%             
%             local_qh = zeros(n_reference_pts,1);
%             for j=1:n_local_dofs
%                 local_dof_j = local_dofs(j);
%                 local_qh = local_qh + ...
%                     phi(:,j)*solution_coeffs(local_dof_j);
%             end
%             
%             local_x = GetPhysicalPoints(elem_id, reference_pts);
%             
%             
%             h3 = plot(local_x, local_qh,c,'linewidth',2);
%             
%             
%             if(elem_id==elem1)
%                 ylim([-0.5,1.5]);
%                 hold on
%             end
%         end
%         
%         h_legend =legend([h1 h2 h3],'epsilon = 0.01','epsilon = 0.1','epsilon = 1');
%         set(h_legend,'FontSize',14)
%         title(sprintf('Diffusion of a Gaussian t=%g (p=0,K=100,dt = 0.0001)',time))
%         drawnow
%         hf = gcf;
%         saveas(hf,sprintf('gaussian_t15_p0_3eps_a1_bneg1.png'))
%     end
% elseif epsilon ==0.1
%     c = 'r';
%     if norm(double(time==[1.5]))
%         h = figure(2);
%         reference_pts = linspace(-1,1,10)';
%         n_reference_pts = length(reference_pts);
%         phi = GetPhi(p, reference_pts);
%         
%         n_local_dofs = p+1;
%         
%         %account for ghost elements
%         if ghost_cell_flag
%             elem1 = 2;     elemK = K-1;
%         else
%             elem1 = 1;     elemK = K;
%         end
%         
%         for elem_id=elem1:elemK
%             local_dofs = DofMap(elem_id,:);
%             
%             local_qh = zeros(n_reference_pts,1);
%             for j=1:n_local_dofs
%                 local_dof_j = local_dofs(j);
%                 local_qh = local_qh + ...
%                     phi(:,j)*solution_coeffs(local_dof_j);
%             end
%             
%             local_x = GetPhysicalPoints(elem_id, reference_pts);
%             
%             
%             h2 = plot(local_x, local_qh,c,'linewidth',2);
%             
%             
%             if(elem_id==elem1)
%                 ylim([-0.5,1.5]);
%                 hold on
%             end
%         end
%         
% %         h_legend =legend(sprintf('epsilon = %g',epsilon));
% %         set(h_legend,'FontSize',14)
% %         title('Diffusion of a Square Pulse  (p=0)')
%         drawnow
%     end
% elseif epsilon == 0.01
%     c= 'm';
%     if norm(double(time==[1.5]))
%         h = figure(2);
%         reference_pts = linspace(-1,1,10)';
%         n_reference_pts = length(reference_pts);
%         phi = GetPhi(p, reference_pts);
%         
%         n_local_dofs = p+1;
%         
%         %account for ghost elements
%         if ghost_cell_flag
%             elem1 = 2;     elemK = K-1;
%         else
%             elem1 = 1;     elemK = K;
%         end
%         
%         for elem_id=elem1:elemK
%             local_dofs = DofMap(elem_id,:);
%             
%             local_qh = zeros(n_reference_pts,1);
%             for j=1:n_local_dofs
%                 local_dof_j = local_dofs(j);
%                 local_qh = local_qh + ...
%                     phi(:,j)*solution_coeffs(local_dof_j);
%             end
%             
%             local_x = GetPhysicalPoints(elem_id, reference_pts);
%             
%             
%             h1 = plot(local_x, local_qh,c,'linewidth',2);
%             
%             
%             if(elem_id==elem1)
%                 ylim([-0.5,1.5]);
%                 hold on
%             end
%         end
%         
% %         h_legend =legend(sprintf('epsilon = %g',epsilon));
% %         set(h_legend,'FontSize',14)
% %         title('Diffusion of a Square Pulse  (p=0)')
%         drawnow
% %         hold all
%     end
% end



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Compare p for constant epsilon 
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if p == 0
            c= 'sb';
            if norm(double(time==[1.5]))
                h = figure(2);
                reference_pts = linspace(-1,1,10)';
                n_reference_pts = length(reference_pts);
                phi = GetPhi(p, reference_pts);
                
                n_local_dofs = p+1;
                
                %account for ghost elements
                if ghost_cell_flag
                    elem1 = 2;     elemK = K-1;
                else
                    elem1 = 1;     elemK = K;
                end
                
                for elem_id=elem1:elemK
                    local_dofs = DofMap(elem_id,:);
                    
                    local_qh = zeros(n_reference_pts,1);
                    for j=1:n_local_dofs
                        local_dof_j = local_dofs(j);
                        local_qh = local_qh + ...
                            phi(:,j)*solution_coeffs(local_dof_j);
                    end
                    
                    local_x = GetPhysicalPoints(elem_id, reference_pts);
                    
                    
                    h1 = plot(local_x, local_qh,c,'linewidth',1);
                    
                    
                    if(elem_id==elem1)
                        ylim([-0.5,1.5]);
                        hold on
                    end
                end
                
                drawnow
            end
        elseif p==1
            c = '^g';
            if norm(double(time==[1.5]))
                h = figure(2);
                reference_pts = linspace(-1,1,10)';
                n_reference_pts = length(reference_pts);
                phi = GetPhi(p, reference_pts);
                
                n_local_dofs = p+1;
                
                %account for ghost elements
                if ghost_cell_flag
                    elem1 = 2;     elemK = K-1;
                else
                    elem1 = 1;     elemK = K;
                end
                
                for elem_id=elem1:elemK
                    local_dofs = DofMap(elem_id,:);
                    
                    local_qh = zeros(n_reference_pts,1);
                    for j=1:n_local_dofs
                        local_dof_j = local_dofs(j);
                        local_qh = local_qh + ...
                            phi(:,j)*solution_coeffs(local_dof_j);
                    end
                    
                    local_x = GetPhysicalPoints(elem_id, reference_pts);
                    
                    
                    h2 = plot(local_x, local_qh,c,'linewidth',1);
                    
                    
                    if(elem_id==elem1)
                        ylim([-0.5,1.5]);
                        hold on
                    end
                end
                
                drawnow
            end
        elseif p ==2
            c= 'r';
            if norm(double(time==[1.5]))
                h = figure(2);
                reference_pts = linspace(-1,1,10)';
                n_reference_pts = length(reference_pts);
                phi = GetPhi(p, reference_pts);
                
                n_local_dofs = p+1;
                
                %account for ghost elements
                if ghost_cell_flag
                    elem1 = 2;     elemK = K-1;
                else
                    elem1 = 1;     elemK = K;
                end
                for elem_id=elem1:elemK
                    local_dofs = DofMap(elem_id,:);
                    
                    local_qh = zeros(n_reference_pts,1);
                    for j=1:n_local_dofs
                        local_dof_j = local_dofs(j);
                        local_qh = local_qh + ...
                            phi(:,j)*solution_coeffs(local_dof_j);
                    end
                    
                    local_x = GetPhysicalPoints(elem_id, reference_pts);
                    
                    h3 = plot(local_x, local_qh,c,'linewidth',1);
                    
                    if(elem_id==elem1)
                        ylim([-0.5,1.5]);
                        hold on
                    end
                end
                
                h_legend =legend([h1 h2 h3],'p=0','p=1','p=2');
                set(h_legend,'FontSize',14)
                title(sprintf('Diffusion of a Gaussian t=%g (epsilon = 0.01,K=100,dt = 0.0001)',time))
                drawnow
                hf = gcf;
                saveas(hf,sprintf('gaussian_eps01_3p_a1_bneg1.png'))
                saveas(hf,sprintf('gaussian_eps01_3p_a1_bneg1.fig'))            
            end
        end