function del_cost = powerquality_gradient_helper(t,X,s,gait,grad_alpha,grad_alphadot,grad_alphaddot,metric,dM,ddM)
% Helper function to calculate the gradient of cost for inertia systems.
% Designed to work with ode45; solves for the gradient of cost at an
% instant of time t.
% Inputs:
%   t: Time period at which gradient of cost is being evaluated.
%   X: Unused, required by ode45
%   s: system struct used by sysplotter
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alphadot: Gradient of shape velocity with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alpha: Gradient of shape position with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

	% Evaluate gradient of shape variables
    grad_alpha_eval = cellfun(@(C) C(t), grad_alpha, 'UniformOutput', false);
    grad_alphadot_eval = cellfun(@(C) C(t), grad_alphadot, 'UniformOutput', false);
    grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
    del_cost = zeros(size(grad_alpha_eval));
    % Get the shape and shape derivative at the current time
	shape = readGait(gait.phi_def,t);
	shapelist = num2cell(shape);
	dshape = readGait(gait.dphi_def,t);
    ddshape = readGait(gait.ddphi_def,t);
    
    % Get mass and partial mass matrices
    M = metric;
    
    % Regular torque calculation
    C = calc_coriolis_matrix(dM,shape,dshape);
    tau = M*ddshape(:) + C;
    
    for i = 1:numel(grad_alpha_eval)
        % Partial of shape variables with respect to fourier coefficient i
        del_shape = grad_alpha_eval{i};
        del_dshape = grad_alphadot_eval{i};
        del_ddshape = grad_alphaddot_eval{i};
    
        % Gradient of torque calculation
        % Start with effect of gradient on M_alpha*alphaddot
        M_temp = zeros(length(shapelist));
        for j = 1:length(shapelist)
            M_temp = M_temp + dM{j}*del_shape(j);
        end
        % Catching for debugging
        try
            M_grad = M_temp*ddshape(:) + M*del_ddshape(:);
        catch
            M_temp
        end
        % Effect of gradient on dM_alphadalpha*alphadot*alphadot
        C1_partialgrad = zeros(length(shapelist));
        C1_shapegrad = zeros(length(shapelist));
        C1_outergrad = zeros(length(shapelist));
        del_dM_alphadalpha = cell(size(shapelist));
        for j = 1:length(shapelist)
            Cj_temp = zeros(length(shapelist));
            for k = 1:length(shapelist)
                Cj_temp = Cj_temp + ddM{j,k}*del_shape(k);
            end
            del_dM_alphadalpha{j} = Cj_temp;
            C1_partialgrad = C1_partialgrad + Cj_temp*dshape(j);
            C1_shapegrad = C1_shapegrad + dM{j}*del_dshape(j);
            C1_outergrad = C1_outergrad + dM{j}*dshape(j);
        end
        C1_grad = (C1_partialgrad + C1_shapegrad)*dshape(:) + ...
            C1_outergrad*del_dshape(:);

        % Effect of gradient on -(1/2)*alphadot'*dM_alphadalpha*alphadot
        C2_grad = zeros(size(shapelist(:)));
        for j = 1:length(shapelist)
            C2_grad(j) = del_dshape(:)'*dM{j}*dshape(:) + ...
                dshape(:)'*del_dM_alphadalpha{j}*dshape(:) + ...
                dshape(:)'*dM{j}*del_dshape(:);
        end
        
        % Gradient of torque
        del_tau = M_grad + C1_grad - (1/2)*C2_grad;
        % Gradient of (P1+P2)^2, aka gradient of v'*tau*v'*tau where v
        % is the shape velocity
        qual_1 = del_dshape(:)'*tau*dshape(:)'*tau + ...
                 dshape(:)'*del_tau*dshape(:)'*tau + ...
                 dshape(:)'*tau*del_dshape(:)'*tau + ...
                 dshape(:)'*tau*dshape(:)'*del_tau;
        % Gradient of P1^2 + P2^2
        qual_2 = (2*del_dshape(:).*dshape(:))'*(tau.*tau) + ...
                 (dshape(:).*dshape(:))'*(2*del_tau.*tau);
        % Gradient of power quality cost wrt this fourier coeff
        del_cost(i) = qual_1 - qual_2;
    end
    del_cost = del_cost(:);
end