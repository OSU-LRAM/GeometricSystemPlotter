function [A,Aeq,C,Ceq]=nonlcon_passive(y,s,n,dimension,lb,ub,direction)
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

global lastCoeffSet lastTransferGrads;
y = lastCoeffSet;
tf = lastTransferGrads;

Aeq = [];
A = getConstraintVector(s,y,n,dimension,lb,ub);

if nargout == 4
   
    C = zeros(size(y,1),size(A,2));
    Ceq = [];
    for i = 1:size(y,1)
        C(i,:) = getGrad(s,y,n,dimension,lb,ub,i,tf{i},A);
    end
    
end

end

function thisGrad = getGrad(s,coeffs,n,dimension,lb,ub,fourierIndex,transferGrads,A0)

    stepSize = 0.01;
    
    newCoeffs = coeffs;
    newCoeffs = newCoeffs + transferGrads*stepSize;
    
    A_new = getConstraintVector(s,coeffs,n,dimension,lb,ub);
    
    thisGrad = (A_new - A0)/stepSize;

end

function A = getConstraintVector(s,coeffs,n,dimension,lb,ub)

    y1 = path_from_fourier(coeffs,n,dimension);
    y2=y1(:);

    w = coeffs(end,2); % Frequency of Fourier transform
    T = 2*pi/w;
    p = makeGait(coeffs);
    ts = linspace(0,T,n)';

    A = [];
    Aeq = [];

    %Imposes the constraint of maximum motor acceleration
    if isfield(s.physics,'maxAcc')
        motorAccelerations = p.ddphi_def{2}(ts);
        worstAcc = max(abs(motorAccelerations));
        %disp(['Max Acc: ',num2str(s.physics.maxAcc),', This Acc: ',num2str(worstAcc)]);
        A = [A;motorAccelerations-s.physics.maxAcc;-s.physics.maxAcc-motorAccelerations];
    end 

    %Imposes the constraint that the gait must be (at least almost) closed.
    %Because the simulation only reaches approximate steady-state, 
    %the system is unlikely to close exactly
    closeTolerance = 0.1;
    closeGap = norm(y1(1,:)-y1(end,:));

    A = [A;closeGap-closeTolerance];

    % A1 and A2 together impose the constraint that all the points forming the gait stay in bounds
    A1=lb-y2;
    A2=y2-ub;

    A = [A;A1;A2]';
    
end
