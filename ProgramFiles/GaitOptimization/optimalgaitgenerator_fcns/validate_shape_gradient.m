function validate_shape_gradient(n,y,g,grad_alphaddot,grad_alphadot,grad_alpha) %#ok<DEFNU>
% Function that helps validate that the gradient of shape position,
% velocity, and acceleration are correctly calculated. Difference between
% the input gradients and calculation-verified gradients are printed to the
% terminal. Should be very close to zero.
% Inputs:
%   n: Number of points at which gait should be evaluated in shape space
%   y: Fourier coefficients that parametrize the gait
%   g: Time period over which gait is executed
%   grad_alphaddot: Gradient of shape acceleration with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alphadot: Gradient of shape velocity with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time
%   grad_alpha: Gradient of shape position with respect to the
%       fourier coefficients, cell array of same dimension as y where each 
%       entry is an array function of time

    % Select 10 random times at which to evaluate the gradient of shape
    t = g*rand(1,10);
    w1 = y(end,1); % Frequency of Fourier transform
    w2 = y(end,2);
    
    % Evaluate the gradient at each random time against a hard-coded
    % gradient calculation
    for i = 1:length(t)
        grad_alpha_eval = cellfun(@(C) C(t(i)), grad_alpha, 'UniformOutput', false);
        grad_alpha1_eval = cell2mat(grad_alpha_eval(:,1));
        grad_alpha2_eval = cell2mat(grad_alpha_eval(:,2));
        grad_alpha1_calc = [1,0;cos(w1*t(i)),0;sin(w1*t(i)),0;cos(2*w1*t(i)),0;sin(2*w1*t(i)),0;cos(3*w1*t(i)),0;sin(3*w1*t(i)),0;cos(4*w1*t(i)),0;sin(4*w1*t(i)),0;0,0];
        grad_alpha2_calc = [0,1;0,cos(w2*t(i));0,sin(w2*t(i));0,cos(2*w2*t(i));0,sin(2*w2*t(i));0,cos(3*w2*t(i));0,sin(3*w2*t(i));0,cos(4*w2*t(i));0,sin(4*w2*t(i));0,0];
        grad_alpha1_err = grad_alpha1_eval - grad_alpha1_calc
        grad_alpha2_err = grad_alpha2_eval - grad_alpha2_calc
        
        grad_alphadot_eval = cellfun(@(C) C(t(i)), grad_alphadot, 'UniformOutput', false);
        grad_alphadot1_eval = cell2mat(grad_alphadot_eval(:,1));
        grad_alphadot2_eval = cell2mat(grad_alphadot_eval(:,2));
        grad_alphadot1_calc = [0,0;-w1*sin(w1*t(i)),0;w1*cos(w1*t(i)),0;-2*w1*sin(2*w1*t(i)),0;2*w1*cos(2*w1*t(i)),0;-3*w1*sin(3*w1*t(i)),0;3*w1*cos(3*w1*t(i)),0;-4*w1*sin(4*w1*t(i)),0;4*w1*cos(4*w1*t(i)),0;0,0];
        grad_alphadot2_calc = [0,0;0,-w2*sin(w2*t(i));0,w2*cos(w2*t(i));0,-2*w2*sin(2*w2*t(i));0,2*w2*cos(2*w2*t(i));0,-3*w2*sin(3*w2*t(i));0,3*w2*cos(3*w2*t(i));0,-4*w2*sin(4*w2*t(i));0,4*w2*cos(4*w2*t(i));0,0];
        grad_alphadot1_err = grad_alphadot1_eval - grad_alphadot1_calc;
        grad_alphadot2_err = grad_alphadot2_eval - grad_alphadot2_calc;
        
        grad_alphaddot_eval = cellfun(@(C) C(t(i)), grad_alphaddot, 'UniformOutput', false);
        grad_alphaddot1_eval = cell2mat(grad_alphaddot_eval(:,1));
        grad_alphaddot2_eval = cell2mat(grad_alphaddot_eval(:,2));
        grad_alphaddot1_calc = [0,0;-w1^2*cos(w1*t(i)),0;-w1^2*sin(w1*t(i)),0;-4*w1^2*cos(2*w1*t(i)),0;-4*w1^2*sin(2*w1*t(i)),0;-9*w1^2*cos(3*w1*t(i)),0;-9*w1^2*sin(3*w1*t(i)),0;-16*w1^2*cos(4*w1*t(i)),0;-16*w1^2*sin(4*w1*t(i)),0;0,0];
        grad_alphaddot2_calc = [0,0;0,-w2^2*cos(w2*t(i));0,-w2^2*sin(w2*t(i));0,-4*w2^2*cos(2*w2*t(i));0,-4*w2^2*sin(2*w2*t(i));0,-9*w2^2*cos(3*w2*t(i));0,-9*w2^2*sin(3*w2*t(i));0,-16*w2^2*cos(4*w2*t(i));0,-16*w2^2*sin(4*w2*t(i));0,0];
        grad_alphaddot1_err = grad_alphaddot1_eval - grad_alphaddot1_calc;
        grad_alphaddot2_err = grad_alphaddot2_eval - grad_alphaddot2_calc;
    end
end