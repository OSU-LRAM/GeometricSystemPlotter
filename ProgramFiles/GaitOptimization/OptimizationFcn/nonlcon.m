function [A,Aeq,gA,gAeq]=nonlcon(y,s,n,dimension,lb,ub,direction,constraint)
%%%%%%%%% 
%
%This function imposes the nonlinear constraint that all the points forming the gait stay within bounds
%
%Inputs:
%
%y: Fourier series coefficients that describe the gait
%s: System file which contains the connection vector field, CCF's and
%   metric data
%n: Number of points used to parametrize the gaits in a direct
%   transcription method
%dimension: Indicates the number of shape variables of the system
%lb: Lower bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
%ub: Upper bound of shape variables for each point which is obtained from the grid inside which an optimal gait is desired
% 
%%%%%%%%%

%% Getting a discrete waypoint series of Fourier Gait.

% The first step is to obtain a direct transciption parametrization of 
% the gait from the fourier series parametrization
y1 = path_from_fourier(y,n,dimension);
y2=y1(:);

%b=length(y2);

%% Constraints for upper bound and lower bound.

% A1 and A2 together impose the constraint that all the points forming the gait stay in bounds
A1=y2+lb;
A2=-y2-ub;

A = [A1;A2];

% the gradient of waypoint with respect to Fourier series.
if nargout > 2
    chy=chy_generator(y,n+1,dimension);
    for i = 1:length(chy)
        chy{i} = [chy{i}; zeros(1,n+1)];
    end
    gA1 = blkdiag(chy{:});
    gA2 = -gA1;
    gA = [gA1 gA2];
end

%% Make sure the frequency doesn't get changed from 2*pi
Aeq = y(end,:)' - 2*pi;
if nargout > 2
    gAeq = zeros(numel(y),2);
end


%% Make sure that the displacement in the other two directions is zero
w1 = y(end,1); % Frequency of Fourier transform
w2 = y(end,2);

% Assign a time period for executing the gait
if (constraint(1) == 1)||(constraint(2) == 1)
    [~, net_disp_opt,~] = evaluate_displacement_and_cost1(s,y);
    net_disp_opt = se2log(net_disp_opt);
end

% Constrain solutions to only displace in the desired direction
if constraint(1) == 1
    for idx=1:3
        if idx ~= direction
            Aeq(end+1) = net_disp_opt(idx);
            if nargout > 2
                Jacob = evaluate_jacobian_fourier(y,s,n,dimension,direction);
                gdisp = [Jacob.disp; 2*pi*ones(1,dimension)];
                gdisp = reshape(gdisp,[size(y,1)*dimension,1]);
                gAeq(:,end+1) = gdisp.';
            end
        end
    end
end

% If optimizing for translation, restrict to zero rotation.
if constraint(2) == 1
    if direction ~=3
        Aeq(end+1) = net_disp_opt(3);
    end
end

% if(nargout > 3)
%     gradA = [];    
%     
%     gradAeq = ones(size(A,1),dimension);
%     
%     jacobiandisp = zeros(n,dimension);
%     for i=2:1:n-1
%         jacobiandisp(i,:)=jacobiandispcalculator3(y(i-1,:),y(i,:),y(i+1,:),ccf(i,:),dimension);
%     end
%     jacobiandisp(1,:)=jacobiandispcalculator3(y(n,:),y(1,:),y(2,:),ccf(1,:),dimension);
%     jacobiandisp(n,:)=jacobiandispcalculator3(y(n-1,:),y(n,:),y(1,:),ccf(n,:),dimension);
% end
% end
