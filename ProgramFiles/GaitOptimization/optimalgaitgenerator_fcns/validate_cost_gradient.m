function validate_cost_gradient(s,n,y,g,p)
% Function to help verify that the cost gradient calculated by
% inertia_gradient_helper is valid. Values returned by function and
% individually calculated costs are printed to the terminal along with the
% difference between the two methods. Should be relatively small for full
% range of fourier coefficients.
% Inputs:
%   s: System struct used by sysplotter.
%   n: Number of points at which gait should be evaluated in the shape
%       space.
%   y: Fourier coefficients that parametrize the gait.
%   g: Time period over which gait is executed.
%   p: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t

% Set the delta in the fourier coefficients between individual cost
% evaluations
fourier_delta = 0.001;
[grad_alphaddot,grad_alphadot,grad_alpha] = shape_grad(n,y,g);
% Go through each of the fourier coefficients and test how applying
% fourier_delta to each of them matches with the cost gradient calculated
% by inertia_gradient_helper
for fourier_test = 1:numel(y)
    % Perturb the fourier coefficient to be tested by fourier_delta and
    % create a parametrization for calculating the point in the gait at a
    % time t
    y2 = y;
    y2(fourier_test) = y2(fourier_test) + fourier_delta;
    w1 = y2(end,1); % Frequency of Fourier transform
    w2 = y2(end,2);
    p2 = makeGait(y2);
            
    % Perturb again in the opposite direction
    y2 = y;
    y2(fourier_test) = y2(fourier_test) - fourier_delta;
    w1 = y2(end,1); % Frequency of Fourier transform
    w2 = y2(end,2);
    p3 = makeGait(y2);
     
    % Get the shape and shape derivative at a random time for each gait
    t = g*rand(1);
	shape = readGait(p3.phi_def,t);
	shapelist = num2cell(shape);
	dshape = readGait(p3.dphi_def,t);
    ddshape = readGait(p3.ddphi_def,t);
    shape_delta = readGait(p2.phi_def,t);
	shapelist_delta = num2cell(shape_delta);
	dshape_delta = readGait(p2.dphi_def,t);
    ddshape_delta = readGait(p2.ddphi_def,t);
    
    % Get mass matrices at both locations
    M = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    M_delta = cellfun(@(C) interpn(s.grid.mass_eval{:},C,...
        shapelist_delta{:},'spline'),s.massfield.mass_eval.content.M_alpha);
    
    % Get partial mass matrices at both locations
    dM_alphadalpha = calc_partial_mass(s,shapelist);
    dM_alphadalpha_delta = calc_partial_mass(s,shapelist_delta);
    
    % Calculate the cost using the two different gaits
    cost = torque_cost(M,dM_alphadalpha,shape,dshape,ddshape);
    cost_delta = torque_cost(M_delta,dM_alphadalpha_delta,shape_delta,dshape_delta,ddshape_delta);
    % Calculate what the gradient of cost is for this particular point
    cost_grad = inertia_gradient_helper(t,[],s,p,grad_alpha,grad_alphadot,grad_alphaddot);
    cost_grad = reshape(cost_grad,size(y));
    cost_grad_rel = cost_grad(fourier_test)
    cost_grad_calc = (cost_delta-cost)/(2*fourier_delta)
    
    % Find what the difference is between cost_grad and the costs evaluated
    % at distance fourier_delta
    err = cost_grad_rel - cost_grad_calc
end

end