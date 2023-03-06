function [grad_alphaddot,grad_alphadot,grad_alpha] = shape_grad(n,y,g)
% Calculates the gradient of the shape position, velocity, and acceleration
% with respect to the fourier coefficients.
% Inputs:
%   n: Number of points composing the gait when using shape-space
%       parametrization
%   y: Set of fourier coefficients parametrizing the gait
%   g: Time period of the gait

% Get the fourier frequency, number of shape variables, and number of
% fourier coefficientsinertia_cost_gradient
w = y(end,:);
dim = size(y,2);
num_coeffs = size(y,1);
% Initialize cell array of function handles to hold the partials of the
% shape variables with respect to the fourier coefficients
grad_alpha = cell(num_coeffs,dim);
grad_alphadot = cell(num_coeffs,dim);
grad_alphaddot = cell(num_coeffs,dim);

for i = 1:num_coeffs 
    for j = 1:dim
        % Coefficient a_0 is a lone scalar, so partials with respect to it
        % are zero for the dotted terms and 1 for alpha
        if i == 1
            grad_alpha{i,j} = @(t) [0*(1:j-1), 1, 0*(j+1:dim)];
            grad_alphadot{i,j} = @(t) zeros(1,dim);
            grad_alphaddot{i,j} = @(t) zeros(1,dim);
            continue
        elseif i == num_coeffs % Partial w.r.t. frequency
            grad_alpha{i,j} = @(t) [0*(1:j-1), t*(-y(2,j)*sin(w(j)*t) + y(3,j)*cos(w(j)*t) - ...
                                               2*y(4,j)*sin(2*w(j)*t) + 2*y(5,j)*cos(2*w(j)*t) - ...
                                               3*y(6,j)*sin(3*w(j)*t) + 3*y(7,j)*cos(3*w(j)*t) - ...
                                               4*y(8,j)*sin(4*w(j)*t) + 4*y(9,j)*cos(4*w(j)*t)), ...
                                    0*(j+1:dim)];
            
            grad_alphadot{i,j} = @(t) [0*(1:j-1), -y(2,j)*(w(j)*t*cos(w(j)*t)+sin(w(j)*t)) + y(3,j)*(-w(j)*t*sin(w(j)*t)+cos(w(j)*t)) + ...
                                                  -2*y(4,j)*(2*w(j)*t*cos(2*w(j)*t)+sin(2*w(j)*t)) + 2*y(5,j)*(-2*w(j)*t*sin(2*w(j)*t)+cos(2*w(j)*t)) + ...
                                                  -3*y(6,j)*(3*w(j)*t*cos(3*w(j)*t)+sin(3*w(j)*t)) + 3*y(7,j)*(-3*w(j)*t*sin(3*w(j)*t)+cos(3*w(j)*t)) + ...
                                                  -4*y(8,j)*(4*w(j)*t*cos(4*w(j)*t)+sin(4*w(j)*t)) + 4*y(9,j)*(-4*w(j)*t*sin(4*w(j)*t)+cos(4*w(j)*t)), ...
                                      0*(j+1:dim)];
            
            grad_alphaddot{i,j} = @(t) [0*(1:j-1), -y(2,j)*w(j)*(-t*w(j)*sin(w(j)*t)+2*cos(w(j)*t)) - y(3,j)*w(j)*(t*w(j)*cos(w(j)*t)+2*sin(w(j)*t)) + ...
                                                   -4*y(4,j)*w(j)*(-2*t*w(j)*sin(2*w(j)*t)+2*cos(2*w(j)*t)) - 4*y(5,j)*w(j)*(2*t*w(j)*cos(2*w(j)*t)+2*sin(2*w(j)*t)) + ...
                                                   -9*y(6,j)*w(j)*(-3*t*w(j)*sin(3*w(j)*t)+2*cos(3*w(j)*t)) - 9*y(7,j)*w(j)*(3*t*w(j)*cos(3*w(j)*t)+2*sin(3*w(j)*t)) + ...
                                                   -16*y(8,j)*w(j)*(-4*t*w(j)*sin(4*w(j)*t)+2*cos(4*w(j)*t)) - 16*y(9,j)*w(j)*(4*t*w(j)*cos(4*w(j)*t)+2*sin(4*w(j)*t)), ...
                                       0*(j+1:dim)];
            continue
        end
        % For partial alpha, a_n is associated with cosine and b_n is
        % associated with sine; a_n terms are every second row entry in y
        % with the b_n terms in between
        if mod(i,2) == 0
            trig = @cos;
        else
            trig = @sin;
        end
        % mult comes from the multiplier of the natural frequency for
        % increasing fourier coefficients
        mult = floor(i/2);
        
        grad_alpha{i,j} = @(t) [0*(1:j-1), trig(mult*w(j)*t), 0*(j+1:dim)];
        % For partial alphadot, a_n is associated with sine and b_n is
        % associated with cosine; the a_n terms are every second row entry 
        % in y with the b_n terms in between
        if mod(i,2) == 0
            trig = @sin;
        else
            trig = @cos;
        end
        
        grad_alphadot{i,j} = @(t) [0*(1:j-1), (-1)^(i-1)*mult*w(j)*trig(mult*w(j)*t), 0*(j+1:dim)];
        
        % For partial alphaddot, a_n is associated with cosine and b_n is
        % associated with sine
        if mod(i,2) == 0
            trig = @cos;
        else
            trig = @sin;
        end
        
        grad_alphaddot{i,j} = @(t) [0*(1:j-1), -mult^2*w(j)^2*trig(mult*w(j)*t), 0*(j+1:dim)];
    end
end
end