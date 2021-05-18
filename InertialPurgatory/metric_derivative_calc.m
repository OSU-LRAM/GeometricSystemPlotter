function [dM,Gamma2] = metric_derivative_calc(r,i,Mp_raw,alpha_field,shv,M)%-----Removed method

%------------------------ Metric Derivation ------------------------------%
gi = inv(M);
%-----Pretty sure Hossein said we'll only use the numerical case-----%
% switch method
%     
%     case 'Analytical'   % If you have your metric equation you can use this option
% 
%         dM{i,1} = dMx_func;
%         dM{i,2} = dMy_func(r(i,1),r(i,2),1,1,1,1);
% %         dM{i,3} = dMz_func(r(i,1),r(i,2),r(i,3),1,1,1,1);
% 
% %         dM_newis1x = dMx_func;
% %         dM_newis1y = dMy_func(r(i-1,1),r(i-1,2),r(i-1,3),1,1,1,1);
% %         dMis1z = dMz_func(r(i-1,1),r(i-1,2),r(i-1,3),1,1,1,1);
% % 
% %         dM_newip1x = dMx_func;
% %         dM_newip1y = dMy_func(r(i+1,1),r(i+1,2),r(i+1,3),1,1,1,1);
% %         dMip1z = dMz_func(r(i+1,1),r(i+1,2),r(i+1,3),1,1,1,1);


%     case 'Numerical'    % When the metric is given for all the configuration you can use this option
        % Find the points around the current point
        
        for j = 1:shv
            h_L_original = [];
            h_R_original = [];
        
            L = r(i,j)-0.001;
            R = r(i,j)+0.001;

            for k = 1:shv
                if k ~=j
                    p1 = r(i,k);
                    p2 = r(i,k);
                else
                    p1 = L;
                    p2 = R;
                end
                h_L_original{k} = [p1];
                h_R_original{k} = [p2];
            end
            
            M_L = Granular_metric_calc2(h_L_original,Mp_raw,alpha_field{:});
            M_R = Granular_metric_calc2(h_R_original,Mp_raw,alpha_field{:});
            
            dg{j} = M_R-M_L;
            
            dU(j) = (h_R_original{j} - h_L_original{j});
            
            dM{i,j} = dg{j}/dU(j);
            
            if isnan(dM{i,j})
                dM{i,j} = zeros(shv,shv);
            end
            
        end
        
        n = size(r,2);
        Gamma2 = zeros (n, n, n);
        for l = 1:n % delta spans the variable labels over which symbols will be calculated
          for ii = 1:n
            for j = 1:n
              for k = 1:n
                Gamma2 (ii, j, l) = Gamma2 (ii, j, l) + ...
                  0.5 * gi (l, k) * ...
                    ((dg{ii}(j,k))/dU(ii) + ...
                     (dg{j}(ii,k))/dU(j)  - ...
                     (dg{k}(ii,j))/dU(k));
              end
            end
          end
          Gamma2 (:, :, l); % display for verification
        end
        
        
        
        
%         h_L_original = [r(i,1)+0.001 r(i,2)];
%         h_R_original = [r(i,1)-0.001 r(i,2)];
%         h_U_original = [r(i,1) r(i,2)+0.001 ];
%         h_D_original = [r(i,1) r(i,2)-0.001];
% %         h_F_original = [r(i,1) r(i,2) r(i,3)+0.001];
% %         h_B_original = [r(i,1) r(i,2) r(i,3)-0.001];
% 
% 
%         % Find the metric coorespoding to the orginal coordinate 
%         M_L_temp = Granular_metric_calc(h_L_original(1),h_L_original(2),Mp_raw,a,b);
%         M_R_temp = Granular_metric_calc(h_R_original(1),h_R_original(2),Mp_raw,a,b);
%         M_U_temp = Granular_metric_calc(h_U_original(1),h_U_original(2),Mp_raw,a,b);
%         M_D_temp = Granular_metric_calc(h_D_original(1),h_D_original(2),Mp_raw,a,b);
% %         M_F_temp = Granular_metric_calc3(h_F_original(1),h_F_original(2),h_F_original(3),Mp_raw,a,b,c);
% %         M_B_temp = Granular_metric_calc3(h_B_original(1),h_B_original(2),h_B_original(3),Mp_raw,a,b,c);
% 
% 
%         % Transform the metric to the local coordinate
%         M_L = M_L_temp;
%         M_R = M_R_temp;
%         M_U = M_U_temp;
%         M_D = M_D_temp;
% %         M_F = M_F_temp;
% %         M_B = M_B_temp;
% 
% 
%         % Find the difference of metrics in x and y direction
%         dg{1} = M_R-M_L;
%         dg{2} = M_U-M_D;
% %         dg{3} = M_F-M_B;
%         
%         % Find the x and y distance in local coordinate
%         dU(1) = (h_R_original(1) - h_L_original(1));
%         dU(2) = (h_U_original(2) - h_D_original(2));
% %         dU(3) = (h_F_original(3) - h_B_original(3));
%         
%         
%         dM{i,1} = dg{1}/dU(1);
% 
%         dM{i,2} = dg{2}/dU(2);
% 
% %         dM{i,3} = dg{3}/dU(3);
% 
%         if isnan(dM{i,1})
%             dM{i,1} = zeros(3,3);
%         end
% 
%         if isnan(dM{i,2})
%             dM{i,2} = zeros(3,3);
%         end
% 
% %         if isnan(dM{i,3})
% %             dM{i,3} = zeros(3,3);
% %         end

        
end
