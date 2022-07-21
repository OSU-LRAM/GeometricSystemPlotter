function [f,g]=solvedifffmincon(y,s,n,dimension,direction,~,~,~)%,lb,ub,writerObj)
%%%%%%%%%%%%%
% This function calculates efficiency (or displacement, if
% that is the objective function) and its gradient with respect to the coefficients obtained
% by the fourier series parametrization
%
% Inputs:
%
% y: Matrix containing the Fourier series coefficients
% s: System file which contains the connection vector field, CCF's and
%   metric data
% dimension: Indicates the number of shape variables of the system
% n: The number of points desired in a direct transcription parametrization
%   of the gaits
% direction: direction in which to optimize motion: 1-x, 2-y, 3-theta
% costfunction: costfunction type to optimize for
% lb: Lower bound of shape variables for each point which is obtained from the grid inside
%   which an optimal gait is desired
% ub: Upper bound of shape variables for each point which is obtained from the grid inside
%   which an optimal gait is desired
%
% Outputs:
%
% f: Objective function value (This is negative of efficiency by default, can be
%   changed to displacement)
% g: Gradient of the objective function
%%%%%%%%%%%%%
global bestCost bestDisp bestEff;

if(nargout > 1)
    [jacobfourier,net_disp_opt,totalstroke] = evaluate_jacobian_fourier(y,s,n,dimension,direction);
else
    p = makeGait(y);
    [~, net_disp_opt, totalstroke] = evaluate_displacement_and_cost1(s,p,[0, 2*pi/y(end,1)],'interpolated','fixed_step');
end
lineint = net_disp_opt(direction);

%% minimizing negative of efficiency(or displacement)
f=-lineint/(totalstroke); % Optimizing for displacement over cost
%     f = -lineint; % Optimizing for displacement only

if abs(f) > bestEff
    bestEff = abs(f);
    bestDisp = abs(lineint);
    bestCost = abs(totalstroke);
end

if nargout>1
    if strcmpi(s.costfunction,'torque') || strcmpi(s.costfunction,'covariant acceleration') || strcmpi(s.costfunction,'acceleration coord') || strcmpi(s.costfunction,'power quality')
        % Return the gradient of efficiency as previously calculated for
        % inertia systems
        g = -jacobfourier.eff; % Optimizing for displacement over cost
        % g = -jacobfourier.disp; % Optimizing for displacement only
    else
        % Return the gradient of efficiency plus row of zeros for frequency
        % terms for drag systems
        g=[-jacobfourier.eff;zeros(1,dimension)]; % Optimizing for displacement over cost
        % g = [-jacobfourier.disp;zeros(1,dimension)]; % Optimizing for displacement only
    end
end
    
    
end