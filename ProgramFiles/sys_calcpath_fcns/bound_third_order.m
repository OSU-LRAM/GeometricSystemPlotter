% estimates an upper bound gait diameter based on a desired third-order bound
% inputs:
    % s:
        % system info
    % shape:
        % vector of "alphas" defining center of gait(s) of interest
    % A_est, cBVI_est:
        % functions taking (s, shape), that define reasonable estimates for
        % the local connection and cBVI within the context of third-order
        % error bound
    % P:
        % upper bound for third order error, expressed as a proportion of
        % the cBVI
    % phi:
        % starting gait phase
% outputs:
    % l:
        % upper bound gait diameter (characteristic length) that will not
        % exceed specified third order error
    % cBVI_fun:
        % function, returning estimate of the cBVI with respect to length
    % third_order_fun:
        % function, returning worst-case third-order estimate with respect
        % to length
function [l, cBVI_fun, third_order_fun] = bound_third_order(s, shape, A_est, cBVI_est, P, phi)
    % sanitize
    if ~exist('phi', 'var')
        phi = 0;
    end
    % anonymous helper functions
    R = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];
    lb_wc = @(x,y) [abs(x(2)*y(3)) + abs(y(2)*x(3));...
                    abs(y(1)*x(3)) + abs(x(1)*y(3)); 0];
                
    % use phi to construct unit tangent vector
    init_tan = [cos(phi); sin(phi)];
    % construct worst-case alpha+beta
    l_sym = sym('l', {'real' 'positive'});
    A_fun = A_est(s, shape);
    a_b_wc = pi/4 * l_sym * (abs(A_fun(l_sym) * R(pi/8) * init_tan) +...
                             abs(A_fun(l_sym) * R(5*pi/8) * init_tan));
    % estimate integral of DA inside circle of diameter l
    cBVI_fun = cBVI_est(s, shape);
    % do worst-case lie bracket
    third_order_wc = lb_wc(a_b_wc, cBVI_fun(l_sym))/2;
    
    % solve for ceiling characteristic length
    l = double(abs(vpasolve(norm(third_order_wc)/norm(cBVI_fun(l_sym)) == P)));
    % get functions for sanity
    third_order_fun = matlabFunction(third_order_wc);
end