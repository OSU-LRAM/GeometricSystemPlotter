function del_cost = accelerationcoord_gradient_helper(t,X,s,gait,grad_alphaddot,metric,dM,ddM)
% Helper function to calculate the gradient of shape-space acceleration
% cost
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

	% Evaluate gradient of shape variables
    grad_alphaddot_eval = cellfun(@(C) C(t), grad_alphaddot, 'UniformOutput', false);
    del_cost = zeros(size(grad_alphaddot_eval));
    % Get the shape derivative at the current time
    
    ddshape_def = readGait(gait.ddphi_def,t);
    actual_size = min(numel(s.grid.eval),numel(ddshape_def));
    ddshape = zeros(size(s.grid.eval));
    
    ddshape(1:actual_size) = ddshape_def(1:actual_size);
   
    % Regular cost calculation
    cost = sqrt(ddshape(:)'*ddshape(:));
    
    for i = 1:numel(grad_alphaddot_eval)
        % Partial of shape acceleration with respect to fourier coefficient i
        del_ddshape = grad_alphaddot_eval{i};
    
        % Gradient of shape-space acceleration cost
        del_cost(i) = del_ddshape(:)'*ddshape(:)+ddshape(:)'*del_ddshape(:);
    end
    del_cost = del_cost(:);
end