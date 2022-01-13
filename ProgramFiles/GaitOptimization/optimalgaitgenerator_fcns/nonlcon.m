function [A,Aeq]=nonlcon(y,s,n,dimension,lb,ub,direction)
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

% % The first step is to obtain a direct transciption parametrization of the gait from the 
% % fourier series parametrization
y1 = path_from_fourier(y,n,dimension);
y2=y1(:);

%b=length(y2);

% A1 and A2 together impose the constraint that all the points forming the gait stay in bounds
A1=y2+lb;
A2=-y2-ub;

A = [A1;A2];


% Make sure the frequency doesn't get changed from 2*pi
Aeq = y(end,:)' - 2*pi;


% Make sure that the displacement in the other two directions is zero
w1 = y(end,1); % Frequency of Fourier transform
w2 = y(end,2);
% Assign a time period for executing the gait
T = 2*pi/w1;
p = makeGait(y);
[~, net_disp_opt] = evaluate_displacement_and_cost1(s,p,[0, T],'interpolated','fixed_step');

% % If optimizing for translation, restrict to zero rotation.
% if direction ~=3
%     Aeq(end+1) = net_disp_opt(3);
% end

% Constrain solutions to only displace in the desired direction
for idx=1:3
    if idx ~= direction
        Aeq(end+1) = net_disp_opt(idx);
    end
end

%
% gradA = [];
% 
% 
% gradAeq = ones(size(A,1),dimension);
% 
% jacobiandisp = zeros(n,dimension);
% for i=2:1:n-1
%     jacobiandisp(i,:)=jacobiandispcalculator3(y(i-1,:),y(i,:),y(i+1,:),ccf(i,:),dimension);
% end
% jacobiandisp(1,:)=jacobiandispcalculator3(y(n,:),y(1,:),y(2,:),ccf(1,:),dimension);
% jacobiandisp(n,:)=jacobiandispcalculator3(y(n-1,:),y(n,:),y(1,:),ccf(n,:),dimension);

end
