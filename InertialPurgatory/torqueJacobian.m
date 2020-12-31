function [newJ,jacob] = torqueJacobian(M,Metric,r,i,T,sign_kappa)

    dr = (cdiff( r' ))'/T;
    
    if dr(i,1) == 0
        dr(i,:) = [1 0];
        
        if i == 1
        M1 = (M{i}+M{i+1})/2;
        elseif i == size(r,1)
            M1 = (M{i}+M{i-1})/2;
        else
            M1 = M{i};
        end
    else
        M1 = (2*M{i}+M{i-1}+M{i+1})/4;
    end

    % Make the jacobian
    norm_u1 = sqrt([dr(i,1) dr(i,2)]*M1*[dr(i,1) dr(i,2)]');
    u1 = [dr(i,1) dr(i,2) 0]/norm_u1;
    
    % Find the normal vector to the u1
    u2_temp = [-u1(2)*sign_kappa u1(1)*sign_kappa 0];%cross(u1,z1);%/norm(cross(u1,z1));
    
    % Find projection of u2 over u1
    proj_u2 = ((u2_temp(1:2)*M1*u1(1:2)')/(u1(1:2)*M1*u1(1:2)'))*u1(1:2);
    
    u2 = u2_temp(1:2) - proj_u2;
    u2 = u2/sqrt(u2*M1*u2');
    
    J = inv([u1(1:2)' u2']);
    newJ = J*J';
    jacob = J;
    
    % Find the derivative of metrix
    shv = 2;
    
%     for j = 1:2
%         h_L_original = [];
%         h_R_original = [];
% 
%         L = r(i,j)-0.001;
%         R = r(i,j)+0.001;
% 
%         for k = 1:shv
%             if k ~=j
%                 p1 = r(i,k);
%                 p2 = r(i,k);
%             else
%                 p1 = L;
%                 p2 = R;
%             end
%             h_L_original{k} = [p1];
%             h_R_original{k} = [p2];
%         end
% 
%         
%         var = [cell2mat(h_L_original); cell2mat(h_R_original)];
%         
%         for ite = 1:shv
%             
%             M2 = Metric(var(ite,1),var(ite,2));
%             % Make the jacobian
%             norm_u1 = sqrt(var(ite,:)*M2*var(ite,:)');
%             u1 = [var(ite,:) 0]/norm_u1;
% 
%             % Find the normal vector to the u1
%             u2_temp = [-u1(2) u1(1) 0];%cross(u1,z1);%/norm(cross(u1,z1));
% 
%             % Find projection of u2 over u1
%             proj_u2 = ((u2_temp(1:2)*M2*u1(1:2)')/(u1(1:2)*M2*u1(1:2)'))*u1(1:2);
% 
%             u2 = u2_temp(1:2) - proj_u2;
%             u2 = u2/sqrt(u2*M2*u2');
% 
%             J = inv([u1(1:2)' u2']);
%             newJ_dir{ite} = J*J';
% 
%         end
% 
%         dg{j} = newJ_dir{2}-newJ_dir{1};
% 
%         dU(j) = (h_R_original{j} - h_L_original{j});
% 
%         dnewJ{1,j} = dg{j}/dU(j);
% 
%         if isnan(dnewJ{1,j})
%             dnewJ{1,j} = zeros(shv,shv);
%         end
% 
%     end
    
    
    