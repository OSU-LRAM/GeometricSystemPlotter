function [grad,delta_vee,acc] = gradientstroke1(i,P,r,M,dM,cons,shv,Gamma,acc,Metric) %-----Removed method as fifth input
    

    T = P.dt;
    % Calculate the Kappa and Gradient of Kappa
    [kappa1,kappa2,F, dkappa1,dkappa2] = curvature_calc4(i,r,M,P,cons,Gamma,dM); %-----Removed method
    
    % Calculate the velocity and Gradient of velocity
    [L_i,L_is1,L_is2,L_ip1, dL_i_x,dL_is1_x,dL_i_y,dL_is1_y] = velocity_calc1(i,r,M,dM,shv,cons);

    % Calculate the tangential acceleation
    a_i = (L_i - L_is1);
    acc(i) = a_i;
    a_s1 = L_is1 - L_is2;
    a_ip1 = (L_ip1 - L_i);
    
    % Calculate the gradient of pathlength
    grad_L = (dL_i_y + dL_is1_y + dL_i_x + dL_is1_x);

    % Calculate the Gradient for acceleration 
    grad_ax1 = (dL_i_x+0*dL_i_y - dL_is1_x-0*dL_is1_y);
    grad_axp1 = -(dL_i_x+0*dL_i_y);
    grad_axs1 = (dL_is1_x+0*dL_is1_y);


    grada = grad_ax1;%(a_i*grad_ax1 + a_ip1*grad_axp1 + a_s1*grad_axs1); %(grad_ax1 + grad_axp1 + grad_axs1);% grad_ax1
    
    % Check type of the cost function whether the velcity is contant or
    % variable
    %-----Pretty sure Hossein said to only use noncte case-----%
%     if strcmp(P.initial_velocity,'cte')
%     
%         T = 1;
%         
%         L = (L_is1+L_i)/2;
%         
%         a_n = kappa2*(L)^2/(T^2);
%         a_t = a_i;
%     
%         delta_vee = ([a_t a_n]*[a_t;a_n])^0.25;
%         
%         delta_vee4 = ([a_t a_n]*[a_t;a_n])*T;
%         
%         grad_an = (dkappa2*L^2)/(T^2) + 2*grad_L*kappa2*L/(T^2);
%         grad_at = grada;
%     
%         grad1 = 0.25*2*[grad_at(1) grad_an(1)]*[a_t;a_n]*delta_vee4^(-3/4);
%         grad2 = 0.25*2*[grad_at(2) grad_an(2)]*[a_t;a_n]*delta_vee4^(-3/4);
%         
%         
%         grad = [grad1 grad2];
%         
%         
%     else    % variable velocity
        
        L = (L_is1+L_i)/2;
        
        a_n = kappa2*(L)^2;
        a_t = a_i;
        
        % calculate the JJ^T
        [newJ,~] = torqueJacobian(M,P.M_t,r,i,1,-1);
    
        delta_vee = ([a_t a_n]*newJ*[a_t;a_n])^0.25;
        
        delta_vee4 = ([a_t a_n]*newJ*[a_t;a_n]);
        
        grad_an = dkappa2*L^2 + 2*grad_L*kappa2*L;
        grad_at = grada;
    
        grad1 = 0.25*2*[grad_at(1) grad_an(1)]*newJ*[a_t;a_n]*delta_vee4^(-3/4);
        grad2 = 0.25*2*[grad_at(2) grad_an(2)]*newJ*[a_t;a_n]*delta_vee4^(-3/4);
        
        grad = [grad1 grad2];
        
%     end
      
    
end